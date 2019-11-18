# STing

#### Table of contents:
<!-- MarkdownTOC levels="1,2,3" -->

- [Requirements](#requirements)
- [Installation](#installation)
- [Quickstart guide](#quickstart-guide)
    - [Typing](#typing)
    - [Gene detection](#gene-detection)
    - [Preparing the database:](#preparing-the-database)
- [STing applications](#sting-applications)
- [STing companion scripts](#sting-companion-scripts)
    - [**`db_util.py`**](#db_utilpy)
    - [**`plot_kmer_depth.R`**](#plot_kmer_depthr)
- [Using STing](#using-sting)
    - [Typing \(MLST\)](#typing-mlst)
    - [Gene Detection](#gene-detection-1)
- [Application extended help](#application-extended-help)
    - [indexer](#indexer)
    - [typer](#typer)
- [**SeqAn version:** 2.1.1](#seqan-version-211)
    - [detector](#detector)

<!-- /MarkdownTOC -->


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

### Gene detection

Preparing directory:

```
mkdir -p STing_demo/my_dbs/card
cd STing_demo/my_dbs/card
```

### Preparing the database:

Downloading the antimicrobial resistance gene database from [CARD](https://card.mcmaster.ca/)

```
wget https://card.mcmaster.ca/latest/data -O card.tgz
tar xvf card.tgz
```

Changing IUPAC extended DNA characters to N on the CARD nucleotide file:

```
sed -e '/^[^>]/s/[^ATGCatgc]/N/g' nucleotide_fasta_protein_homolog_model.fasta > card.fasta
```

Creating STing db config file:

```
echo -e "[loci]\ncard\tcard.fasta\n" > config.txt
```

The new config file should look like this:
```
[loci]
card    card.fasta
```

Creating the STing database index:

```
indexer -c config.txt -p db/index -m GDETECT
```

Downloading a WGS sample to analyze:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_2.fastq.gz
```

Run STing detector:

```
cd ../../
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

## STing companion scripts

STing includes two companion scripts:

### **`db_util.py`**
This script provides a set of utilities to download databases from PubMLST and build STing indices from them.

```

This script provides a set of utilities to download databases from PubMLST and build STing indices from them.

optional arguments:
  -h, --help          show this help message and exit
  --version           show program's version number and exit

subcommands:
  {list,query,fetch}
    list              List all the available databases at PubMLST
    query             Search a database in PubMLST
    fetch             Fetch a database from PubMLST
```

To list the available MLST schemes and their last update time, use the `list` sub-command:

```bash
./scripts/db_util.py list
# #       Database        #Profiles       Retrieved       DB_URL
# 1       Achromobacter spp.      476     2019-11-18      https://pubmlst.org/achromobacter
# 2       Acinetobacter baumannii#1       2058    2019-11-18      https://pubmlst.org/abaumannii/

```

To search the available PubMLST schemes by search term (e.g. 'cholera'), use the `query` sub-command:

```bash
./scripts/db_util.py query "cholera"
# 1        -b 984     2019-11-18      https://pubmlst.org/vcholerae/
# 2       Vibrio cholerae#2       422     2019-11-18      http://pubmlst.org/vcholerae
```

Finally, to download and build a PubMLST database, use the `fetch` sub-command:

```bash
./scripts/db_util.py fetch -q "Vibrio cholerae" -b -o testdb
# Database: "Vibrio cholerae"
#  Fetching allele sequences:
#  - https://pubmlst.org/data/alleles/vcholerae/adk.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/adk.fa
#  - https://pubmlst.org/data/alleles/vcholerae/gyrB.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/gyrB.fa
#  - https://pubmlst.org/data/alleles/vcholerae/mdh.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/mdh.fa
#  - https://pubmlst.org/data/alleles/vcholerae/metE.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/metE.fa
#  - https://pubmlst.org/data/alleles/vcholerae/pntA.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/pntA.fa
#  - https://pubmlst.org/data/alleles/vcholerae/purM.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/purM.fa
#  - https://pubmlst.org/data/alleles/vcholerae/pyrC.tfa -> /storage/aroon/sting/databases/vibrio_cholerae/pyrC.fa
#  Fetching profiles:
#  - https://pubmlst.org/data/profiles/vcholerae.txt -> /storage/aroon/sting/databases/vibrio_cholerae/profile.txt
# Building STing index:
# /data/home/achande3/bin/indexer -c /storage/aroon/sting/databases/vibrio_cholerae/config.txt -p /storage/aroon/sting/databases/vibrio_cholerae/db/index
# Loading sequences from sequences files:

# #       Seqs.   File
# 1       130     /storage/aroon/sting/databases/vibrio_cholerae/adk.fa
# 2       151     /storage/aroon/sting/databases/vibrio_cholerae/gyrB.fa
# 3       164     /storage/aroon/sting/databases/vibrio_cholerae/mdh.fa
# 4       329     /storage/aroon/sting/databases/vibrio_cholerae/metE.fa
# 5       185     /storage/aroon/sting/databases/vibrio_cholerae/pntA.fa
# 6       128     /storage/aroon/sting/databases/vibrio_cholerae/purM.fa
# 7       259     /storage/aroon/sting/databases/vibrio_cholerae/pyrC.fa

# Total loaded sequences: 1346

# Creating and saving ESA index from loaded sequences...
# Index successfuly created!

```

### **`plot_kmer_depth.R`** 

This script Generates __k__-mer depth plots from output depth files generated by the `typer` and `detector` applications.

#### Requirements:

`plot_kmer_depth.R` requires the following R packages:

- argparser
- ggsci
- gridExtra
- RColorBrewer
- stringr
- svglite
- tidyverse

By default, `plot_kmer_depth.R` will try to install automatically the required packages to the personal R library directory (usually something like `~/R/x86_64-pc-linux-gnu-library/3.4`).  

## Using STing

Next subsections present a brief description of the MLST analysis and gene detection workflows of STing. For the detailed help and usage of each application, please check the last section of this document ([Application extended help](#extended-help-of-the-applications)).
In order to test the different STing tools, the folder ```test``` is included at the STing root directory. It contains the example files described below (configuration, profiles, and sequences) in addition to some read samples.

### Typing (MLST)

Two main steps are required to execute an MLST analysis with STing: *Database indexing* and *ST prediction*

#### Database indexing

1. Create a config file that contains the path to loci and allelic profiles files. The format details for each file are the following:

    ##### Config file
    
    A tab separated file that follows this format:
    
    ```
    [loci]
    locus1   relative/path/to/locusFile1.fa
    locus2   relative/path/to/locusFile2.fa
    [profile]
    profile  relative/path/to/profileFile.txt
    ```
    
    Blank lines and comments (lines starting with ```#```) in this file, will be ignored. Paths are relative to the config file itself.
    This is an example of a configuration file for the species *Neisseria spp.* (```test/mlst/pubmlst_db_files/Neisseria_spp/config.txt```):
    
    ```
    [loci]
    abcZ  abcZ.fa
    adk adk.fa
    aroE  aroE.fa
    fumC  fumC.fa
    gdh gdh.fa
    pdhC  pdhC.fa
    pgm pgm.fa
    [profile]
    Neisseria_spp profile.txt
    ```
    ##### Allele sequence file
    
    A standard multi-FASTA file in which the id of each sequence consists of a locus name and an allele number separated by ```_```:
    
    ```
    >abcZ_1
    TTTGATACTGTTGCCGA...
    >abcZ_2
    TTTGATACTGTTGCCGA...
    ```
    
    ##### Profile file
    
    A tab separated file which contains each ST and its corresponding allelic profile (```test/mlst/pubmlst_db_files/Neisseria_spp/profile.txt```):
    
    ```
    ST  abcZ  adk  aroE  fumC  gdh  pdhC  pgm  clonal_complex
    1   1     3    1     1     1    1     3    ST-1 complex/subgroup I/II
    2   1     3    4     7     1    1     3    ST-1 complex/subgroup I/II
    3   1     3    1     1     1    23    13   ST-1 complex/subgroup I/II
    4   1     3    3     1     4    2     3    ST-4 complex/subgroup IV
    ...
    ```
    
    ST and loci columns (columns 1 to 7) are required in this file.
    
2. Build the database using the ```indexer``` tool:

    ```
    indexer -m <MODE> -c <CONFIG_FILE> -p <PREFIX_FILE>
    ```
    
    Example:
    
    ```
    indexer -m MLST -c test/mlst/pubmlst_db_files/Neisseria_spp/config.txt -p test/mlst/dbs/Neisseria_spp/db
    ```
    
    The command above will create an index for MLST analysis (mode MLST specified by ```-m MLST```) using the MLST database files defined in the config file ```test/mlst/pubmlst_db_files/Neisseria_spp/config.txt```. As a result, the indexer will generate a series of files named with the prefix ```db``` inside the directory ```test/mlst/dbs/Neisseria_spp```.
    
    The output of the command looks like this:
    
    ```
    Loading sequences from sequences files:
    N  Loci  #Seqs.  File
    1  abcZ  1019    Neisseria_spp/abcZ.fa
    2  adk   788     Neisseria_spp/adk.fa
    3  aroE  1048    Neisseria_spp/aroE.fa
    4  fumC  1125    Neisseria_spp/fumC.fa
    5  gdh   1080    Neisseria_spp/gdh.fa
    6  pdhC  1025    Neisseria_spp/pdhC.fa
    7  pgm   1097    Neisseria_spp/pgm.fa
    Total sequences loaded: 7182
    Loading the profiles file... Done!
    Creating and saving ESA index from loaded sequences... Done!
    Index created successfully!
    ```
         
#### ST prediction

Use the ```typer```tool to predict the ST of a read set:

```
typer -x <INDEX_PREFIX_FILENAME> -1 <FASTQ1> -2 <FASTQ2> -k <KMER_LENGTH> -s <SAMPLE_NAME>
```

Example:

```
typer -x test/mlst/dbs/Neisseria_spp/db -1 test/mlst/samples/ERR017011_1.fastq.gz -2 test/mlst/samples/ERR017011_2.fastq.gz -s ERR017011
```

The command above will predict the ST of the read sample called ```ERR017011``` (```-s```), which comprises the input files ```test/mlst/samples/ERR017011_1.fastq.gz``` and ```test/mlst/samples/ERR017011_2.fastq.gz``` (specified by ```-1``` and ```-2```), using the index located at the directory ```test/mlst/dbs/Neisseria_spp``` with the prefix ```db``` (```-x test/mlst/dbs/Neisseria_spp/db```).

The output of the previous command looks like this:

```
Sample  Line_type       Status  ST      abcZ    adk     aroE    fumC    gdh     pdhC    pgm     Total_k-mers  Total_reads      Input_files
ERR017011       allelic_profile st_predicted    5       1       1       2       1       3       2       3     76657    3597    ERR017011_1.fastq.gz,ERR017011_2.fastq.gz
```

By default, the typer application will send the header to ```stderr```, and the prediction result to ```stdout```. You can use the option ```-o``` to save the whole results (header and prediction) to a file, e.g., ```-o my_results.tsv```.

### Gene Detection

Two main steps are required to detect genes using STing: *Database building* and *Detection*.

#### Database building

1. Create a config file that contains the path to gene files. The format details are the following:

    ##### Config file
    A tab separated file with the following format:
    
    ```
    [loci]
    gene1   /path/to/geneFile1.fa
    gene2   /path/to/geneFile2.fa
    ```
    
    Blank lines and comments (lines starting with ```#```) in this file, will be ignored.
    Note that there are no ```[profile]``` section for a configuration file for gene detection. If the file contains this section, the ```indexer``` tool will show an error message.
    This is an example of a configuration file of AMR genes (```test/amr/amr_db_files/set_01/config.txt```):
    
    ```
    [loci]
    erm erm.fasta
    ksg ksg.fasta
    pen pen.fasta
    qac qac.fasta
    aac2ic  aac2ic.fasta
    aph6id  aph6id.fasta
    bl2a_1  bl2a_1.fasta
    ermb  ermb.fasta
    mepa  mepa.fasta
    pbp2b pbp2b.fasta
    pbp2x pbp2x.fasta
    tetpa tetpa.fasta
    ```
    
    ##### Gene sequence file
    
    A standard multi-FASTA file in which the id is the name of the gene. In case of having genes with the same name, you should add a number to the name separated by ```_```:
    
    ```
    >pen_1
    TTTGATACTGTTGCCGA...
    >pen_2
    TTTGATACTGTTGCCGA...
    ```
    
2. Build the database using the ```indexer``` tool:

    ```
    indexer -m <MODE> -c <CONFIG_FILE> -p <PREFIX_FILE>
    ```
    
    Example:
    
    ```
    indexer -m GDETECT -p databases/amr -c test/amr/amr_db_files/set_01/config.txt
    ```
    
    The command above will create an index for gene detection (mode GDETECT in the config file ```test/amr/amr_db_files/set_01/config.txt```. As a result, the indexer will create a series of files named with the prefix ```amr``` inside the directory ```databases```.
    
    The output of the command looks like this:
    
    ```
    Loading sequences from sequences files:
    N     Loci    #Seqs.  File
    1     aac2ic  1       set_01/aac2ic.fasta
    2     aph6id  1       set_01/aph6id.fasta
    3     bl2a_1  1       set_01/bl2a_1.fasta
    4     erm     1       set_01/erm.fasta
    5     ermb    1       set_01/ermb.fasta
    6     ksg     2       set_01/ksg.fasta
    7     mepa    1       set_01/mepa.fasta
    8     pbp2b   1       set_01/pbp2b.fasta
    9     pbp2x   1       set_01/pbp2x.fasta
    10    pen     4       set_01/pen.fasta
    11    qac     1       set_01/qac.fasta
    12    tetpa   1       set_01/tetpa.fasta
    Total loaded sequences: 16
    Creating and saving ESA index from loaded sequences...
    Index created successfuly!
    ```
    
#### Detection

You must use the ```detector``` tool to detect genes in a read set:

```
detector -x <INDEX_PREFIX_FILENAME> -1 <FASTQ1> -2 <FASTQ2> -k <KMER_LENGTH>
```

Example:

```
detector -x databases/amr -1 test/amr/samples/GCF_000008805.fasta.40.1.fq.gz -2 test/amr/samples/GCF_000008805.fasta.40.2.fq.gz -k 30 -s GCF_000008805
```

The command above will detect presence/absence (1/0) of the genes from the database, in the read sample called ```GCF_000008805``` (```-s```) which comprises the input files ```test/amr/samples/GCF_000008805.fasta.40.1.fq.gz``` and ```test/amr/samples/GCF_000008805.fasta.40.2.fq.gz``` (specified by ```-1``` and ```-2```), using the index located at the directory ```databases``` with the prefix ```amr``` (```-x databases/amr```). Additionally, the tool will use *k*-mers of size 30 (```-k 30```) to process the input reads.

The output of the previous command looks like this:

```
Sample  Line_type ermC  ksga1 ksga2 pbp2b pen1  pen2  pen3  pen4  qacE1 Total_hits  Total_kmers Total_reads Input_files 
GCF_000008805 presence  1 1 0 1 1 1 1 1 1 288106170395  1396  GCF_000008805.fasta.40.1.fq.gz,GCF_000008805.fasta.40.2.fq.gz
```

By default, the detector application will send the header to ```stderr```, and the prediction result to ```stdout```. 

## Application extended help

### indexer

#### Synopsis

**indexer** [OPTIONS] -c <CONFIG_FILE>  

#### Description

STing **indexer** creates indexes (DBs) required by the STing **typer** and **detector** programs for loci-based typing analysis and detecting genes, respectively, from NGS raw sequence reads.

<dl>
  <dt>-h, --help</dt>
    <dd>Display the help message.</dd>
  <dt>--version</dt>
    <dd>Display version information.</dd>
  <dt>-c, --config-file CONFIG_FILE</dt>
    <dd>A tab delimited file whith names and paths to the typing scheme files (see the FILE FORMAT DETAILS section below).</dd>
  <dt>-p, --db-prefix PREFIX</dt>
    <dd>Filename prefix for the DB files to be created. You can specify a folder structure here to store your DB at a particular location, e.g., path/to/my/db/prefix. Default: name of the config file CONFIG_FILE</dd>
  <dt>-m, --mode MODE</dt>
    <dd>Indexing mode. Valid options: MLST, GDETECT. Select MLST to create a database for MLST analysis or GDETECT to create a database for gene detection. Default: _MLST_.</dd>
</dl>

#### File format details

##### CONFIG_FILE

A tab separated file with the name and location of files for creating a DB.
Config files for MLST DBs (MLST mode) must have two sections: **[loci]** that describes names and paths to alleles sequence files for each locus, and **[profile]** that describes the name and path to the profile file.
Config files for gene detection DBs (GDETECT mode), only require the [loci] section.
An example of a config file for a MLST DB is as follows:

  ```
  [loci]
  locus1 relative/path/to/locusFile1
  locus2 relative/path/to/locusFile2
  locusN relative/path/to/locusFileN
  [profile]
  profile relative/path/to/profileFile
  ```

Paths are relative to the config file itself. Blank lines and comments (lines starting with '#') will be ignored.

##### Allele sequence file

A standard multi-FASTA file (.fa or .fasta) in which each sequence description must be the locus name and the allele number separated by '_':

  ```
  >abcZ_1
  TTTGATACTGTTGCCGA...
  >abcZ_2
  TTTGATACTGTTGCCGA...
  ```
  
##### Profile file

A tab separated file that contains the ST and the corresponding allelic profile:

  ```
  ST abcZ adk aroE fumC gdh pdhC pgm clonal_complex
  1 1 3 1 1 1 1 3 ST-1 complex/subgroup I/II
  2 1 3 4 7 1 1 3 ST-1 complex/subgroup I/II
  3 1 3 1 1 1 23 13 ST-1 complex/subgroup I/II
  4 1 3 3 1 4 2 3 ST-4 complex/subgroup IV
  ```
  
#### Version

**Last update:** 09/11/2018  
**indexer version:** 0.23.2  
**SeqAn version:** 2.1.1

--------------------------------------------------------------------------------

### typer

#### Synopsis

**typer** -x <INDEX_PREFIX_FILENAME> -1 <FASTQ1> [OPTIONS]  

#### Description

STing **typer** is an ultrafast assembly- and alignment-free program for sequence typing directly from whole-genome raw sequence reads. STing **typer** is based on k-mer frequencies, and works with locus-based typing schemes like those defined in the traditional MLST method and its derivatives (e.g, rMLST, cgMLST or wgMLST). STing **typer** requires an index (DB) created with the STing **indexer** program.

<dl>
    <dt>-h, --help</dt>
    <dd>Display the help message.</dd>
    <dt>--version</dt>
    <dd>Display version information.</dd>
</dl>

##### Required input parameters:

<dl>
    <dt>-x, --index-prefix INDEX_PREFIX_FILENAME</dt>
    <dd>Index prefix filename.</dd>
    <dt>-1, --fastq-1-files FASTQ1</dt>
    <dd>Files with #1 mates, paired with files in _FASTQ2_. Valid file extensions are _.fq_, _.fastq_ (uncompressed fastq), and _.gz_ (gzip'ed fastq).</dd>
</dl>

##### Input options:

<dl>
    <dt>-2, --fastq-2-files FASTQ2</dt>
    <dd>Files with #2 mates, paired with files in _FASTQ1_. Valid file extensions are _.fq_, _.fastq_ (uncompressed fastq), and _.gz_ (gzip'ed fastq).</dd>
    <dt>-s, --sample-name SAMPLE_NAME</dt>
    <dd>Name of the sample to be analized.</dd>
    <dt>-k, --kmer-length KMER_LENGTH</dt>
    <dd>Length of the k-mers to process the input reads. Default: _30_.</dd>
    <dt>-n, --n-top-alleles N_TOP_ALLELES</dt>
    <dd>Number of top alleles by k-mer frequencies, from which best alleles will be called. Default: _2_.</dd>
</dl>

##### Output options:

<dl>
    <dt>-c, --kmer-counts</dt>
    <dd>Select and print the normalized k-mer counts at each locus (k-mer hits divided by allele length).</dd>
    <dt>-a, --allele-cov</dt>
    <dd>Select to calculate and print the percent of the length of each allele that is covered by k-mer hits. Note that this calculation will require more execution time.</dd>
    <dt>-d, --kmer-depth</dt>
    <dd>Select to calculate and print the average k-mer depth of each allele. Note that this calculation will require more execution time.</dd>
    <dt>-t, --k-depth-file KMER_DEPTH_FILENAME</dt>
    <dd>Output filename to save the detailed per-base k-mer depth data.</dd>
    <dt>-e, --ext-k-depth-file EXT_KMER_DEPTH_FILENAME</dt>
    <dd>Output filename to save the extended per-base k-mer depth data.</dd>
    <dt>-o, --output-file OUTPUT_FILENAME</dt>
    <dd>Output filename to save the typing results (instead of stdout).</dd>
    <dt>-y, --print-tidy</dt>
    <dd>Select to print results in a tidy format.</dd>
    <dt>-v, --verbose</dt>
    <dd>Select to print informative messages (to stderr).</dd>
</dl>

##### Presets:

<dl>
    <dt>--fast</dt>
    <dd>(Default). Select to set the _fast_ mode to call alleles based solely on k-mer frequencies. The best allele of each locus is that with the highest k-mer hit frequency.</dd>
    <dt>--sensitive</dt>
    <dd>Select to set the _sensitive_ mode to call alleles based on k-mer frequencies and coverage information. The best allele of each locus is that with the highest k-mer hit frequency and the highest allele coverage. Allele ties are solved selecting the allele with the minimum k-mer depth standard deviation.</dd>
</dl>

#### Note on FASTQ1 and FASTQ2

_FASTQ1_ and _FASTQ2_ can be comma-separated lists (no whitespace) and can be specified many times, e.g., -1 file1.fq,file2.fq -1 file3.fq

#### Output conventions

##### Status values:

<dl>
    <dt>st_predicted</dt>
    <dd>ST inferred from k-mer hits and profiles table.</dd>
    <dt>no_kmer_hits</dt>
    <dd>There are no k-mer hits for one or more of the loci to infer a profile and its associated ST. Probable causes:</dd>
    <dd>1) k-mer length not adequate (too long),</dd>
    <dd>2) low quality data/too many N's in the data.</dd>
    <dt>no_st_in_table</dt>
    <dd>There is no ST associated to the inferred allelic profile. Probable causes:</dd>
    <dd>1) k-mer length not adequate (too short),</dd>
    <dd>2) this could be a new allelic profile.</dd>
</dl>

##### Other values/symbols:

<dl>
    <dt>NA</dt>
    <dd>No ST associated or no k-mer hits for any allele of a locus.</dd>
    <dt>*</dt>
    <dd>Low confidence. Allele predicted with a length coverage below 100%. Probable causes:</dd>
    <dd>1) No enough k-mer hits to cover the whole length of the allele</dd>
    <dd>2) No enough k-mer depth in some part along the range [10, length-10] bp of the allele to consider it as covered.</dd>
</dl>

#### Version

**Last update:** 02/28/2019  
**typer version:** 0.24.2  
**SeqAn version:** 2.1.1
--------------------------------------------------------------------------------

### detector

#### Synopsis
**detector** -x <INDEX_PREFIX_FILENAME> -1 <FASTQ1> [OPTIONS]  

#### Description
STing **detector** is an ultrafast assembly- and alignment-free program for detecting genes directly from NGS raw sequence reads. STing **detector** is based on k-mer frequencies. STing **detector** requires an index (DB) created with the STing **indexer** program (using the GDETECT mode).
<dl>
  <dt>-h, --help</dt>
    <dd>Display the help message.</dd>
  <dt>--version</dt>
    <dd>Display version information.</dd>
</dl>
#### Required input parameters:
<dl>
  <dt>-x, --index-prefix INDEX_PREFIX_FILENAME</dt>
    <dd>Database prefix filename.</dd>
  <dt>-1, --fastq-1-files FASTQ1</dt>
    <dd>Files with #1 mates, paired with files in _FASTQ1_.</dd>
</dl>
#### Input options:
<dl>
  <dt>-2, --fastq-2-files FASTQ2</dt>
    <dd>Files with #2 mates, paired with files in _FASTQ2_.</dd>
  <dt>-s, --sample-name SAMPLE_NAME</dt>
    <dd>Name of the sample to be analized.</dd>
  <dt>-k, --kmer-length KMER_LENGTH</dt>
    <dd>Length of the k-mers to process the input reads. Default: _30_.</dd>
  <dt>-r, --threshold THRESHOLD</dt>
    <dd>Minimum length coverage (%) required to consider a gene as present in a sample. In range [1.0..100.0]. Default: _75_.</dd>
</dl>

##### Output options:

<dl>
  <dt>-c, --kmer-counts</dt>
    <dd>Select to print the number of k-mer matches at each gene.</dd>
  <dt>-p, --kmer-perc</dt>
    <dd>Select to print the percentage of k-mer matches from the total of processed k-mers.</dd>
  <dt>-g, --gene-cov</dt>
    <dd>Select to print the percent of the gene length that is covered by the corresponding k-mer matches.</dd>
  <dt>-d, --kmer-depth</dt>
    <dd>Select to print the mean k-mer depth of each gene.</dd>
  <dt>-y, --print-tidy</dt>
    <dd>Select to print results in a tidy format.</dd>
  <dt>-t, --k-depth-file KMER_DEPTH_FILENAME</dt>
    <dd>Output filename to save the detailed per-base k-mer depth data.</dd>
  <dt>-v, --verbose</dt>
    <dd>Select to print informative messages (to stderr). Default non verbose.</dd>
</dl>
_FASTQ1_ and _FASTQ2_ can be comma-separated lists (no whitespace) and can be specified many times. E.g. -1 file1.fq,file2.fq -1 file3.fq

#### Version

**Last update:** 10/13/2018  
**detector version:** 0.23.3  
**SeqAn version:** 2.1.1
