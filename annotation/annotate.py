import json
import os
import sqlite3
import sys
import textwrap

import disease
import argparse

db = './hic.db'
const_loop_range = 1e5


def get_enhancer_by_gene(gene):
    """
    Param:
        Gene Name

    Returns:
        [] - Corresponding Enhancer
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select chrom, txStart, txEnd from enhancer where gene like '%{}%' ".format(
        gene)
    cursor.execute(sql)
    enhancer_logs = cursor.fetchall()
    enhancer_info = ''
    for enhancer_log in enhancer_logs:
        enhancer_info += enhancer_log[0] + ':' + \
            str(enhancer_log[1]) + '-' + str(enhancer_log[2]) + ','
    conn.commit()
    conn.close()
    return enhancer_info[:-1]


def get_super_enhancer_by_gene(gene):
    """
    Param:
        Gene Name

    Returns:
        [] - Corresponding Super Enhancer
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select  chrom, txStart, txEnd from superEnhancer where cloestactivegene like '%{}%' ".format(
        gene)
    cursor.execute(sql)
    super_enhancer_logs = cursor.fetchall()
    super_enhancer_info = ''
    for super_enhancer_log in super_enhancer_logs:
        super_enhancer_info += super_enhancer_log[0] + ':' + str(
            super_enhancer_log[1]) + '-' + str(super_enhancer_log[2]) + ','
    conn.commit()
    conn.close()
    return super_enhancer_info[:-1]


def get_genes_by_name(name):
    """
    Param:
        Gene Name

    Returns:
        [] - Gene Info
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    # get gene info
    sql = "select name1, chrom, strand, txStart, txEnd, name2 ,count(distinct name2) from refGene where name2 like '%{}%' group by  name2".format(
        name)
    cursor.execute(sql)
    # ('NM_001136131', 'chr21', '-', 26174731, 26465317, 'APP', 1)
    gene_logs = cursor.fetchall()
    conn.commit()
    conn.close()
    return gene_logs


def get_disease_by_genes(q_genes):
    """
    Param:
        Gene Names

    Returns:
        [] - Corresponding Disease
    """
    names = ','.join(map(str, q_genes.keys()))
    diseases = disease.query_disease_by_gene(names)
    if diseases is None:
        return None
    d_infos = {}
    for d in diseases:
        if not isinstance(d, dict):
            continue
        g_name = d['gene_symbol']
        idx = q_genes[g_name]
        if idx in d_infos:
            d_infos[idx] += d['disease_name'] + ','
        else:
            d_infos[idx] = d['disease_name'] + ','
    return d_infos


def get_promoter_by_gene(gene):
    """
    Param:
        Gene Name

    Returns:
        [] - Corresponding Promoter
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select  chrom, txStart, txEnd, mode from epd where gene like '%{}%' ".format(
        gene)
    cursor.execute(sql)
    promoter_logs = cursor.fetchall()
    promoter_info = ''
    for promoter_log in promoter_logs:
        promoter_info += promoter_log[0] + ':' + \
            str(promoter_log[1]) + '-' + str(promoter_log[2]) + ','
    conn.commit()
    conn.close()
    return promoter_info[:-1]


def get_tf_target_by_gene(name):
    """
    Param:
        Gene Name

    Returns:
        [] - Corresponding TF and Target
    """
    tf_target = {}
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select *, count(distinct Target) from TRRUST where TF='{}' group by Target".format(name)
    cursor.execute(sql)
    targets = cursor.fetchall()
    target_info = ''
    for target in targets:
        if target[2] == 'Activation':
            mode = '+'
        elif target[2] == 'Repression':
            mode = '-'
        else:
            mode = '?'
        target_info += target[1] + "({})".format(mode) + ','
    target_info = target_info[:-1]
    tf_target['target'] = target_info
    sql = "select *, count(distinct TF) from TRRUST where Target='{}' group by TF".format(name)
    cursor.execute(sql)
    tfs = cursor.fetchall()
    tf_info = ''
    for tf in tfs:
        if tf[2] == 'Activation':
            mode = '+'
        elif tf[2] == 'Repression':
            mode = '-'
        else:
            mode = '?'
        tf_info += tf[0] + "({})".format(mode) + ','
    tf_info = tf_info[:-1]
    tf_target['TF'] = tf_info
    conn.commit()
    conn.close()
    return tf_target


def get_gene_locus_by_name(name):
    """
    Param:
        Gene Name

    Returns:
        [] - Corresponding Gene Locus
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select chrom, txStart, txEnd from refGene where name2='{}'".format(
        name)
    cursor.execute(sql)
    locus = ''
    logs = cursor.fetchall()
    if len(logs) >= 1:
        log = logs[0]
        locus = log[0] + ':' + format(log[1], ',') + '-' + format(log[2], ',')
    conn.commit()
    conn.close()
    return locus


def generate_res(gene_logs):
    res = []
    genes = {}
    for idx, gene_log in enumerate(gene_logs):
        card = dict()
        card['id'] = gene_log[0]
        card['locus'] = gene_log[1] + ':' + \
            format(gene_log[3], ',') + '-' + format(gene_log[4], ',')
        card['strand'] = gene_log[2]
        card['name'] = gene_log[5]
        # get enhancer info
        card['enhancer'] = get_enhancer_by_gene(gene_log[5])
        card['promoter'] = get_promoter_by_gene(gene_log[5])
        card['super_enhancer'] = get_super_enhancer_by_gene(gene_log[5])
        genes[gene_log[5]] = idx
        card['diseases'] = ''
        tf_target = get_tf_target_by_gene(gene_log[5])
        card['TF'] = tf_target['TF']
        card['target'] = tf_target['target']
        res.append(card)
    disease_infos = get_disease_by_genes(genes)
    if disease_infos is not None:
        for k in disease_infos.keys():
            res[k]['diseases'] = disease_infos[k][:-1]
    return res


def get_left_close_gene(chr_x, x_s, x_e):
    """
    Param:
        chr_x  - chromosome index  chr16
        x_s - start
        x_e - end

    Returns:
        [] - upstream corresponding gene info
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select name1, chrom, strand, txStart, txEnd, name2, count(distinct name2) from refGene where chrom='{}' and txEnd < {} and txEnd > {} group by name2".format(
        chr_x, x_s, x_s - const_loop_range)
    cursor.execute(sql)
    logs = cursor.fetchall()
    conn.commit()
    conn.close()
    return logs


def get_right_close_gene(chr_x, x_s, x_e):
    """
        Param:
            chr_x  - chromosome index  chr16
            x_s - start
            x_e - end

        Returns:
            [] - downstream corresponding gene info
        """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "select name1, chrom, strand, txStart, txEnd, name2, count(distinct name2) from refGene where chrom='{}' and txStart > {} and refGene.txStart < {} group by name2".format(
        chr_x, x_e, x_e + const_loop_range)
    cursor.execute(sql)
    logs = cursor.fetchall()
    conn.commit()
    conn.close()
    return logs


def get_cross_gene_by_locus(chr_x, x_s, x_e):
    """
        Param:
            chr_x  - chromosome index  chr16
            x_s - start
            x_e - end

        Returns:
            [] - corresponding gene info that located in TAD
        """
    conn = sqlite3.connect(db)

    cursor = conn.cursor()
    sql = "select name1, chrom, strand, txStart, txEnd, name2, count(distinct name2) from refGene where chrom='{}' and txStart <= {} and txEnd >= {} group by name2".format(
        chr_x,
        x_e,
        x_s)
    cursor.execute(sql)
    logs = cursor.fetchall()
    conn.commit()
    conn.close()
    return logs


def query_by_locus(chr_x, x_s, x_e):
    genes = []
    c_logs = get_cross_gene_by_locus(chr_x, x_s, x_e)
    genes += c_logs
    # cha zhao zuo he you ji yin
    if len(genes) < 1:
        l_logs = get_left_close_gene(chr_x, x_s, x_e)
        r_logs = get_right_close_gene(chr_x, x_s, x_e)

        genes += l_logs
        genes += r_logs
    if len(genes) > 0:
        res = generate_res(genes)
    else:
        res = None
    return res


def annotation(_c, _s, _e, _r):
    _s *= _r
    _e *= _r
    _c = "chr" + str(_c)
    res = query_by_locus(_c, _s, _e)
    return res


def process_bulk(_arr, _c, _r, _o):
    with open(_o, 'w') as f:
        for i in range(0, len(_arr) - 1, 2):
            _s, _e = int(_arr[i]), int(_arr[i + 1])
            res = annotation(_c, _s, _e, _r)
            f.write(
                "# Annotation the TAD of {}:{}-{}, resolution:{} \n".format(_c, _s, _e, _r))
            res = json.dumps(res, indent=4)
            f.write(res)
            f.write('\n')
    return


def process_file(_path, _c, _r, _o):
    fo = open(_o, 'w')
    with open(_path, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            line = line.split('\t')
            if len(line) < 2:
                print("file error")
                sys.exit()
            _s = int(line[0])
            _e = int(line[1])
            res = annotation(_c, _s, _e, _r)
            fo.write(
                "# Annotation the TAD of {}:{}-{}, resolution:{} \n".format(_c, _s, _e, _r))
            res = json.dumps(res, indent=4)
            fo.write(res)
            fo.write('\n')
    fo.close()
    return


def process_one(_o, _c, _r, _s, _e):
    with open(_o, 'w') as f:
        res = annotation(_c, _s, _e, _r)
        f.write(
            "# Annotation the TAD of {}:{}-{}, resolution:{} \n".format(_c, _s, _e, _r))
        res = json.dumps(res, indent=4)
        f.write(res)
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawDescriptionHelpFormatter,
        description= textwrap.dedent('''\
        A simple tools to annotate Human TAD!
        -------------------------------------
            Use example: python annotate.py -c 16 -r 40000 -b 1170,1181,677,684
            Use example: python annotate.py -c 16 -r 40000 -s 1170 -e 181
            Use example: python annotate.py -c 16 -r 40000 -f './test_TAD'
        ''')
    )
    parser.add_argument(
        '-c',
        type=int,
        help="index of Chromosome",
        required=True
    )

    parser.add_argument(
        '-s',
        type=int,
        help="index of starting bin of TAD")

    parser.add_argument(
        '-e',
        type=int,
        help="index of ending bin of TAD")

    parser.add_argument(
        '-r',
        type=int,
        required=True,
        help="resolution(unit b) of Hi-C contact matrix")

    parser.add_argument(
        '-o',
        type=str,
        default="./annotation_out",
        help="output file to store result,the default is ./annotation.out")

    parser.add_argument(
        '-b',
        type=str,
        help="a list of TAD. For example, 1,5,9,13 contains TAD (1,5) and TAD (9, 13)")

    parser.add_argument(
        '-f',
        type=str,
        help="a file record TAD. One TAD per line tab separation")
    args = parser.parse_args()
    c, o, s, e, r, b, f = args.c, args.o, args.s, args.e, args.r, args.b, args.f
    if b is not None:
        print("bulk processing")
        b = b.split(',')
        process_bulk(_arr=b, _c=c, _o=o, _r=r)
        print("bulk result already write to {}".format(o))
        sys.exit()
    if f is not None:
        print("file processing")

        process_file(_path=f, _o=o, _c=c, _r=r)
        print("file process already write to {}".format(o))
        sys.exit()
    if s is None:
        print("No Start Input")
        sys.exit()
    if e is None:
        print("No End Input")
        sys.exit()
    print("process {} {} {}".format(c, s, e))
    _res = process_one(_o=o, _r=r, _c=c, _e=e, _s=s)
    print(_res)
