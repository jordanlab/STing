# Typing (MLST)

Two main steps are required to execute an MLST analysis with STing: *Database indexing* and *ST prediction*

## Database indexing

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
         
## ST prediction

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


# typer

## Synopsis

**typer** -x <INDEX_PREFIX_FILENAME> -1 <FASTQ1> [OPTIONS]  

## Description

STing **typer** is an ultrafast assembly- and alignment-free program for sequence typing directly from whole-genome raw sequence reads. STing **typer** is based on k-mer frequencies, and works with locus-based typing schemes like those defined in the traditional MLST method and its derivatives (e.g, rMLST, cgMLST or wgMLST). STing **typer** requires an index (DB) created with the STing **indexer** program.

<dl>
    <dt>-h, --help</dt>
    <dd>Display the help message.</dd>
    <dt>--version</dt>
    <dd>Display version information.</dd>
</dl>

### Required input parameters:

<dl>
    <dt>-x, --index-prefix INDEX_PREFIX_FILENAME</dt>
    <dd>Index prefix filename.</dd>
    <dt>-1, --fastq-1-files FASTQ1</dt>
    <dd>Files with #1 mates, paired with files in _FASTQ2_. Valid file extensions are _.fq_, _.fastq_ (uncompressed fastq), and _.gz_ (gzip'ed fastq).</dd>
</dl>

### Input options:

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

### Output options:

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

### Presets:

<dl>
    <dt>--fast</dt>
    <dd>(Default). Select to set the _fast_ mode to call alleles based solely on k-mer frequencies. The best allele of each locus is that with the highest k-mer hit frequency.</dd>
    <dt>--sensitive</dt>
    <dd>Select to set the _sensitive_ mode to call alleles based on k-mer frequencies and coverage information. The best allele of each locus is that with the highest k-mer hit frequency and the highest allele coverage. Allele ties are solved selecting the allele with the minimum k-mer depth standard deviation.</dd>
</dl>

## Note on FASTQ1 and FASTQ2

_FASTQ1_ and _FASTQ2_ can be comma-separated lists (no whitespace) and can be specified many times, e.g., -1 file1.fq,file2.fq -1 file3.fq

## Output conventions

### Status values:

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

### Other values/symbols:

<dl>
    <dt>NA</dt>
    <dd>No ST associated or no k-mer hits for any allele of a locus.</dd>
    <dt>*</dt>
    <dd>Low confidence. Allele predicted with a length coverage below 100%. Probable causes:</dd>
    <dd>1) No enough k-mer hits to cover the whole length of the allele</dd>
    <dd>2) No enough k-mer depth in some part along the range [10, length-10] bp of the allele to consider it as covered.</dd>
</dl>
