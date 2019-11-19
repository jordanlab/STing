# STing: Ultrafast sequence typing and gene detection from NGS raw reads

Rapid sequence typing.  Identification of bacterial sequence type (ST) is of utmost importance for molecular epidemiology and outbreak control.  The most common approach to typing is multilocus sequence typing (MLST) that defines an ST based on seven to nine housekeeping loci.  Availability of whole-genome sequencing has allowed increased resolution in MLST by incorporating more loci into the scheme.  Such “super-MLST” schemes include ribosomal MLST (rMLST, ~50 loci), core-genome MLST (cgMLST, >1,500 loci) and whole-genome MLST (>2,000 loci).  While these super-MLST schemes provide substantially higher resolution, the classical analysis of whole genome sequence entails substantial downstream bioinformatics analyses before the ST can be determined.  IHRC-GIT ABiL’s tool STing (Sequence Typing) is an ultrafast utility that that identifies an ST directly from raw sequencing reads using a novel, alignment-free algorithm.  STing correctly assigns sequence types with 100% accuracy in traditional MLST and near 100% accuracy in rMLST and cgMLST.    

## Requirements


- Linux OS
- gcc >= 4.8 and gcc <= 5.3 </br>
    You may use gcc >= 5.3 to compile STing but you will see multiple warning messages. Those warning messages will not interfere with the binary generation.
- autotols (Ubuntu: autotools)
    - Autoconf >= 2.69
    - Automake >= 1.15.1
    - Libtool >= 2.4.6
- autotools developer headers (Ubuntu: autotools-dev)
- zlib >= 1.2.8 (Ubuntu: zlib)
- zlib developer headers  (Ubuntu: zlib1g-dev)

## Installation

### Using the pre-built static binaries

```bash
wget https://github.com/jordanlab/STing/releases/download/1.0.0/sting_v1.0.0.tar.gz
tar xfv sting_v1.0.0.tar.gz
export PATH=$PWD/sting:$PATH
```

### From source: 
```bash
./autogen.sh
./configure
make
make install
```


By default, `make install` will install all the files in ```/usr/local/bin```, ```/usr/local/lib``` etc.  You can specify an installation prefix other than ```/usr/local``` using the ```--prefix``` options from ```./configure```, for instance ```./configure --prefix=$HOME```.  Please check all the available options of ```./configure``` by executing ```./configure --help```.

## STing applications

STing has three applications:

* **```indexer```**
: Creates the indexes (databases) required to execute a locus-based typing analysis and detect genes.

* **```typer```**
: Predicts STs of a read sample based on an index built from a species-specific locus-based typing scheme (allelic profile table and sequences of observed alleles).

* **```detector```**
: Detects genes in a read sample based on an index built from sequences of a set of genes of interest.

To explore the usage and available options for each tool, run the corresponding application using the option ```-h``` or ```--help```:

```
indexer -h
typer -h
detector -h
```
