# Gene Detection

Two main steps are required to detect genes using STing: *Database building* and *Detection*.

## Database building

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
    
## Running Gene Detection

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

# detector

## Synopsis
**detector** -x <INDEX_PREFIX_FILENAME> -1 <FASTQ1> [OPTIONS]  

## Description
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

## Input options:
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

### Output options:

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
