# indexer

## Synopsis

**indexer** [OPTIONS] -c <CONFIG_FILE>  

## Description

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

## File format details

### CONFIG_FILE

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

### Allele sequence file

A standard multi-FASTA file (.fa or .fasta) in which each sequence description must be the locus name and the allele number separated by '_':

  ```
  >abcZ_1
  TTTGATACTGTTGCCGA...
  >abcZ_2
  TTTGATACTGTTGCCGA...
  ```
  
### Profile file

A tab separated file that contains the ST and the corresponding allelic profile:

  ```
  ST abcZ adk aroE fumC gdh pdhC pgm clonal_complex
  1 1 3 1 1 1 1 3 ST-1 complex/subgroup I/II
  2 1 3 4 7 1 1 3 ST-1 complex/subgroup I/II
  3 1 3 1 1 1 23 13 ST-1 complex/subgroup I/II
  4 1 3 3 1 4 2 3 ST-4 complex/subgroup IV
  ```
  
