# LPAD
Detection of Topological Associated Domains from Hi-C Data using Network Construction and Label propagation

## 1. Hi-C Data 

In our experiment, we used 3 three normalized Hi-C matrix pre-processed by Bing Ren's Lab. Data is available at : http://chromosome.sdsc.edu/mouse/hi-c/download.html.
> As of January 1, 2023, this website is dead temporarily, you can use our backup at: https://mailnankaieducn-my.sharepoint.com/:f:/g/personal/2120210438_mail_nankai_edu_cn/ElUG7nVOdzpKlakCfYZPPSQB355VjGguq_t93bVLvlY-wg?e=acg2YC

Besides,we also used a simulated data extract from chr6 on mouse cortex cell line, normalized Matrix is in the folder  `data`.

## 2.Usage:

run the PLA_2.py

```sh
python PLA_2.py -h
```

Parameters are as follows:

 * f: the path of a `intra-chromosomal` Hi-C matrix seperated by Tab with N by N shaped
 * w: window-size,the default is 5.
 * c: the restart probability, the default is 0.9
 * o: the storage path and filename of result, the default is ./out
 * p (<b>optional</b>): the Penalty coefficient, the default is 1.
 * k (<b>optional</b>): the top k, the default is 0.6.

### Requirements
 * numpy
 * matplotlib
 * seaborn
 * colormap
 * pyrwr: we use this package download from https://github.com/jinhongjung/pyrwr and modify the source code. <b>you do not need to install it.</b>
 * scipy
 * tdqm
 * pandas

 ## 3.Output

LPAD produces a files. Note that,The output file will be named by the user in which list the TADs extracted from the input data. The first column is the start bin, and the second column is the stop bin.

# Annotation TAD
```
cd ./annotation
python annotation -h
```

Use example: python annotate.py -c 16 -r 40000 -b 1170,1181,677,684

Use example: python annotate.py -c 16 -r 40000 -s 1170 -e 181

Use example: python annotate.py -c 16 -r 40000 -f './test_TAD'


optional arguments:
  -h, --help  show this help message and exit
  -c C        index of Chromosome
  -s S        index of starting bin of TAD
  -e E        index of ending bin of TAD
  -r R        resolution(unit b) of Hi-C contact matrix
  -o O        output file to store result,the default is ./annotation.out
  -b B        a list of TAD. For example, 1,5,9,13 contains TAD (1,5) and TAD (9, 13)
  -f F        a file record TAD. One TAD per line tab separation
