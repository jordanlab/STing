#!/usr/bin/env bash

# process fasta files
# parameters: {db id} {file #1} {file #2}

# location of svg images will be a sub directory under file #1

# configuration variables

DB_HOME=${DB_HOME:-$HOME/pubmlst_dbs}


db_id="$1"
file1="$2"
file2="$3"

if ! [[ "$db_id" =~ ^[0-9]+$ ]]; then
	echo "ERROR: db id must be an integer, was: $db_id"
	exit 1
fi

if [ ! -f "$file1" ]; then
	echo "ERROR: file ${file1} not found"
	exit 1
fi

if [ ! -f "$file2" ]; then
	echo "ERROR: file ${file2} not found"
	exit 1
fi

file1=`realpath "$file1"`
file2=`realpath "$file2"`

tmp_dir=`echo "$file1" |sed 's/.fastq.*//i'`
tmp_dir=`mktemp -d "${tmp_dir}-XXXXX"`

bin_base=`realpath "$0"`
bin_base=`dirname "$bin_base"`


# step 1: Run STing
cd $tmp_dir
sting_command=`"$bin_base/get_sting_command.py" \
	-f "$bin_base/../data/organism-typing_schemes.tsv" \
	-i "$db_id" \
	-l "${file1},${file2}" \
	-b "$HOME/.local/STing/bin" \
	-d "$DB_HOME"`

$sting_command


# step 2: create report

STing_Report \
	--input `ls -tr *.k-depth.log|tail -n 1` \
	--genes `ls -tr *.log|grep -v k-depth.log|tail -n 1` \
	--base `ls -tr *.k-depth.log|tail -n 1|sed 's/.k-depth.log//'` \
	--plot \
	--no-depths \
	--parallel | \
parallel



# step: 3 output JSON

cat *json

