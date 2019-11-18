#!/bin/bash

binary="ntyper";

flags="-g -O1";

# Compiling
command="g++ -DHAVE_CONFIG_H -I. -I./src  -I include -W -Wall -pedantic -std=c++11 ${flags} -c src/KmerMatcherOptions.cpp -o src/KmerMatcherOptions.o";
echo $command;
$command;

command="g++ -DHAVE_CONFIG_H -I. -I./src  -I include -W -Wall -pedantic -std=c++11 ${flags} -c src/KmerMatcher.cpp -o src/KmerMatcher.o"
echo $command;
$command;

command="g++ -DHAVE_CONFIG_H -I. -I./src  -I include -W -Wall -pedantic -std=c++11 ${flags} -c src/Typer.cpp -o src/Typer.o"
echo $command;
$command;

# Linking
command="g++ -I include -W -Wall -pedantic -std=c++11 ${flags} -static-libstdc++ -static-libgcc -o ${binary} src/Typer.o src/KmerMatcherOptions.o src/KmerMatcher.o -lz -lrt"
echo $command;
$command;

# command="./${binary}"
# command="./${binary} --help"
# echo $command;
# $command;

# command="/usr/bin/time -v ./styper -d dbs/N_meningitidis_db -1 test/fastq_files/sample_11.1.fq -2 test/fastq_files/sample_11.2.fq"
# echo $command;
# $command;

# command="/usr/bin/time -v ./styper -d dbs/N_meningitidis_db -1 test/fastq_files/ERR028666_1.fq -2 test/fastq_files/ERR028666_2.fq"
# echo $command;
# $command;