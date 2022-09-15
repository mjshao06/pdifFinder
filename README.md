# PdifFinder
Author:     Shao Mengjie; Liang Qian

Email:      1437819081@qq.com; norah-liang@dmicrobe.com

Institute:  Hangzhou Digital-Micro Biotech Co.,Ltd

This program is designed for annotation of antimicrobal resistance(AMR), pdif site and pdif-ARGs module in bacteria.

### Install:
PdifFinder is a python3.X script, running on linux. 
You should install BLAST and add it in environment variable, you can download from `https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/`. BLAST version is 2.10.1 in pdifFinder.

* One:
  You can download from github by `git clone https://github.com/mjshao06/pdifFinder.git`. Then execute `cd pdifFinder`. Last execute `pip install .`.
* Two:
  You can install PdifFinder from [PyPI](https://pypi.org/project/PdifFinder) by `pip install PdifFinder`.


### Run:
PdifFinder can accept FASTA and GENBANK format file(single or multi sequences in one file). Attention on GENBANK format file, it should follow standard format.
There are three input parameter, "-i" means FASTA, "-g" means GENBANK, "-d" means input dir contains FASTA or GENBANK.
* Simply, you can just run:
```
pdifFinder -n FASTA -o outdir
pdifFinder -g GENBANK -o outdir
pdifFinder -d inputdir -o outdir
```
* For more parameter, you can run:
```
pdifFinder -h
```
* Here are some import parameter:

parameter  | description
---- | -----
--nucleotide(-n) | FASTA file
--genbank(-g) | GENBANK file
--indir(-d) | input dirname
--resultdir(-o) | output dirname
--map(-c) | output graph format,default is circle

### Databases:
Here are databases structure:
<pre>
  .
  ├── AMRDB
  │   ├── sequence.fasta         Resistance gene reference sequences in FASTA format
  │   │                     sequence id must be database name~~~gene~~~accession~~~description,
  │   │                     eg:  ncbi~~~1567214_ble~~~NG_047553.1~~~BLEOMYCIN BLMA family bleomycin binding protein
  │   ├── Res.nhr
  │   ├── Res.nin
  │   └── Res.nsq
  │   └── Res.ndb
  │   └── Res.not
  │   └── Res.nto
  │   └── Res.ntf
  └── data
      ├── redundant.seed.fa  Pdif site reference sequences in FASTA format
                             sequence id must be database >plasmid accession number in NCBI
                             eg: >KY984047.1_1 ACTGCGCATAAGAGATTTTATGTTAAAT
      ├── pdifdatabase.fasta ALL pdif sites from 481 plasmids
      └── genecolor.txt
</pre>      
### Output:

filename  | description
---- | -----
AMRgene.txt | resistance gene annotation
pdif_site.txt | pdif site annotation
pdifmodule_list.txt | pdif-ARGs module annotation
pdifmoduleseq.fasta | pdif-ARGs module sequence
pdifmodule.svg | pdif-ARGs figure
plasmid.html | circular graph for above features
