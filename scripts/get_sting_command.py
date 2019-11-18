#!/usr/bin/env python3

# ##############################################################################
#  Author:   Hector Fabio Espitia-Navarro
#            Georgia Institute of Technology
#  Version:  1.1
#  Date:     04/15/2019
# ##############################################################################

import argparse
import os
# import pprint as pp
import re
import sys
import textwrap

def setup_argument_parser():
    class CustomHelpFormatter(
                argparse.ArgumentDefaultsHelpFormatter,
                argparse.RawTextHelpFormatter,
                argparse.RawDescriptionHelpFormatter):
            pass
    
    exampleText = textwrap.dedent('''\
        examples:  
          Two input read files:
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 1 -l "sample_1.fq,sample_2.fq" -b "path/to/bin" -d "path/to/dbs" 
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 142 -l "sample_123_1.fastq.gz,sample_123_2.fastq.gz" -b "path/to/bin" -d "path/to/dbs" 
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 154 -l "some/very/long/path/to/fastq/sample/files/sample.N431_1.fastq.gz,some/very/long/path/to/fastq/sample/files/sample.N431_2.fastq.gz" -b "path/to/bin" -d "path/to/dbs" 
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 155 -l "path/to/samples/SampleName.N342_S1_L001_R1_001.fastq.gz,path/to/samples/SampleName.N342_S1_L001_R2_001.fastq.gz" -b "path/to/bin" -d "path/to/dbs"
          One input read file:
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 1 -l "sample_1.fq" -b "path/to/bin" -d "path/to/dbs" 
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 142 -l "sample.T562_1.fastq.gz" -b "path/to/bin" -d "path/to/dbs" 
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 154 -l "some/very/long/path/to/fastq/sample/files/sample#1N23_1.fastq.gz" -b "path/to/bin" -d "path/to/dbs"
            ./get_sting_command.py -f organism-typing_schemes.tsv -i 155 -l "path/to/samples/SampleName.N342_S1_L001_R1_001.fastq.gz" -b "path/to/bin" -d "path/to/dbs"
            ''')
    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''
            This script returns a STing command line, given: a file with the corresponding organism dbs and schemes, an 
            organism db id from the file, a list of input read files, and the paths to the STing binaries and main db 
            directories.'''),
        epilog=exampleText, 
        formatter_class=CustomHelpFormatter)
    
    parser.add_argument('--version', action='version', 
        version='%(prog)s\n  version:\t{}\n  last update:\t{}'.format(
            '1.1',
            '04/15/2019'))
    
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    
    dataFileHelpText = textwrap.dedent('''\
        Tab separated file with 4 columns: 
        - id: id of the organism database
        - display_name: organism db to be shown in GUI 
        - scheme: name of the scheme (e.g., MLST, cgMLST, AMR)
        - original_name: original organism name in the origin website (e.g. in PubMLST)''')
    required.add_argument('-f', '--data_file', required=1, help=dataFileHelpText)
    required.add_argument('-i', '--db_id', required=1, help='Id of the organism db (column 1 in data file)')
    required.add_argument('-l', '--file_list', required=1, help='List of input files (comma separated, no spaces)')
    required.add_argument('-b', '--bin_path', required=1, help='Path to the STing binaries location')
    required.add_argument('-d', '--dbs_path', required=1, help='Path to the STing main database directory location')
    
    parser._action_groups.append(optional)
    return parser

def loadDbFile(filename):
    dbFile = open(filename, 'r')
    dbDict = {}
    with dbFile:
        for line in dbFile:
            fields     = line.rstrip().split('\t')
            id         = fields[0]
            if not id.isnumeric():  # Skip header
                continue
            dbDict[id] = {
                'display_name'  : fields[1],
                'scheme'        : fields[2],
                'original_name' : fields[3]#,
                # 'norm_name'     : '',
                # 'norm_scheme'   : '',
                # 'index_path'    : ''
            }
    return dbDict

def normalizeName(name):
    normName = name.lower()
    normName = re.sub(r' |#|-', '_', normName)
    normName = re.sub(r'([^0-9])\.', r'\1', normName)
    normName = re.sub(r'/','-', normName)
    normName = re.sub(r'\(|\)','', normName)
    normName = re.sub(r'_{2,}','_', normName)
    return normName
    
def getIndexPath(dbsPath, normDbName, normScheme):
    indexPath = os.path.join(dbsPath, normScheme, normDbName, 'db', 'index')
    return indexPath

def getSampleName(inputFiles):
    illuminaRegex     = re.compile("([A-Za-z0-9-_.#]+)_S[0-9]+_L[0-9]+_R([1-2])_[0-9]+\.(fastq|fq)(\.gz)*") # https://regex101.com/r/3aEF1g/7
    normalRegex       = re.compile("([A-Za-z0-9-_.#]+)_R*([0-9]{1})\.(fastq|fq)(\.gz)*")                    # https://regex101.com/r/3aEF1g/4
    inFile1           = inputFiles.split(',')[0].rstrip()
    dirname, filename = os.path.split(os.path.abspath(inFile1))
    illuminaMatchObj  = illuminaRegex.fullmatch(filename)
    normalMatchObj    = normalRegex.fullmatch(filename)
    sampleName        = ''
    
    if illuminaMatchObj:
        sampleName = illuminaMatchObj.group(1)
    elif normalMatchObj:
        sampleName = normalMatchObj.group(1)
    
    return sampleName

def main():
    parser = setup_argument_parser()
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    
    dataFile      = args.data_file
    inputFileList = args.file_list
    dbId          = args.db_id
    binPath       = args.bin_path
    dbsPath       = args.dbs_path
    
    dbDict        = loadDbFile(dataFile)
    # pp.pprint(dbDict)
    dbName        = dbDict[dbId]['original_name']
    scheme        = dbDict[dbId]['scheme']
    normDbName    = normalizeName(dbName)
    normScheme    = normalizeName(scheme)
    inFile1       = ''
    inFile2       = ''
    
    if inputFileList != "":
        fields  = inputFileList.split(',')
        inFile1 = fields[0].rstrip()
        if len(fields) == 2:
            inFile2 = fields[1].rstrip()
    
    inputFiles = '-1 {}'.format(inFile1)
    if inFile2 != '':
        inputFiles = '{} -2 {}'.format(inputFiles, inFile2)
    
    sampleName = getSampleName(inputFileList)
    
    typerOptions    = '-k 30 -s {} -o {}.log -t {}.k-depth.log'.format(sampleName, sampleName, sampleName)
    detectorOptions = '-k 30 -s {} -o {}.log -t {}.k-depth.log'.format(sampleName, sampleName, sampleName)
    
    tool    = os.path.join(binPath, 'typer')
    options = typerOptions
    if normScheme == 'amr':
        tool    = os.path.join(binPath, 'detector')
        options = detectorOptions
    
    indexPath = getIndexPath(dbsPath, normDbName, normScheme)
    
    command = '{} -x {} {} {}'.format(tool, indexPath, inputFiles, options)
    
    print(command)
    
if __name__ == '__main__':
    main()
