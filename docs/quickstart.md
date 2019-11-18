# Typing 

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

# Gene detection

Preparing directory:

```
mkdir -p STing_demo/my_dbs/card
cd STing_demo/my_dbs/card
```

## Preparing the database:

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

## Run STing detector:

```
cd ../../
typer -x my_dbs/neisseria_spp/db/index -1 ERR026529_1.fastq.gz -2 ERR026529_2.fastq.gz -s ERR026529 -c -a -d -t ERR026529.depth.tsv -y -o ERR026529.results.tsv --sensitive
```