# varViewer

## Brief description

 varViewer is designed for displaying the reads mapping information of the variant (including SNP/InDel, Exome CNV and Fusion gene).  Specially,  the varViewer can help clinical researcher to review the mutation sites to filter false positive mutations. 

​    varViewer depends on hitlib and boost library.

## How to download and install
On a linux or mac system with git and g++ installed varViewer can be downloaded and installed as follows:

```
 git clone  git@github.com:huangfei304/varViewer.git
 cd varViewer/src
 mkdir build && cd build
 cmake ..
 make
```

## result
#### 1. SNP/InDel plot

For SNP/InDel plot, at least 3 infiles is need.

```
1. the bam format mapping file with index
2. the genome sequence file with index
3. the variant information file(format: chr<tab>start<tab>end<tab>varType)

example:
varViewer -r hg19.fa -b test.bam -i test.var.infor -t SNV -p outprefix
```

 ![snp](/data/snp.png)

Note:  for Pair-End reads (+:lightpink, -: lightslateblue), x54: reads 54 copies, Lowercase base: low quaility base (<20).

#### 2. CNV plot

For CNV plot, at least 6 files

```
1. 1. the bam format mapping file with index
2. the genome sequence file with index
3. the reference sample bam in this batch (one sample per line)
4. the SNP/InDel variant information file(format: chr<tab>pos<tab>ref<tab>allale<tab>infor)
5. the CNV information file(format: chr<tab>varType<tab>gene<tab>exon_start<tab>exon_end)
6. Mapability file


example:
varViewer -r hg19.fa -b test.bam -n refer.bam.list -i test.CNV.infor -v test.vcf.infor -t CNV -a target.bed -k Mapability.bigWig.gz -p outprefix
```

![cnv](/data/cnv.png)

Note: The normailized reads depth and SNP reads rate both do not supported the BRCA1 23 exon Duplication. The 23th exon Duplication maybe false-positive.

#### 3. Fusion gene plot

For Fusin gene plot, at least 2 files

```
1. the bam format mapping file with index
2. the genome sequence file with index
3. the variant information file(format: gene1<tab>gene2<tab>gene1_bp<tab>gene2_bp)
e.g.: MUC16   AATF    chr19:9021043   chr17:35327512

example:
varViewer -r hg19.fa -b test.bam -i test.sv.infor -t SV -p outprefix
```

![SV_ALK_EML4](/data/SV_ALK_EML4.svg)




## Citing and references
1. coming soon...

