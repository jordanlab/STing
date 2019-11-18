# STing

#### Table of contents:
- [Documentation](#documention)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quickstart guide](#quickstart-guide)
    - [Typing](#typing)
    - [Gene detection](#gene-detection)
    - [Preparing the database:](#preparing-the-database)
- [STing applications](#sting-applications)


## Documentation

See the full documentation at: https://sting.readthedocs.io/ 


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

```
./autogen.sh
./configure
make
make install
```

By default, `make install` will install all the files in ```/usr/local/bin```, ```/usr/local/lib``` etc.  You can specify an installation prefix other than ```/usr/local``` using the ```--prefix``` options from ```./configure```, for instance ```./configure --prefix=$HOME```.  Please check all the available options of ```./configure``` by executing ```./configure --help```.

## Quickstart guide

### Typing 

Preparing directory:

```
mkdir STing_demo
cd STing_demo
```

Downloading the *Neisseria spp.* MLST database from PubMLST.org and build a STing index from it:

```
db_util.py fetch --query "Neisseria spp." --out_dir my_dbs --build_index
```

Downloading a WGS sample to analyze:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_2.fastq.gz
```

Run STing typer:

```
typer -x my_dbs/neisseria_spp/db/index -1 ERR026529_1.fastq.gz -2 ERR026529_2.fastq.gz -s ERR026529 -c -a -d -t ERR026529.depth.tsv -y -o ERR026529.results.tsv --sensitive
```


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

