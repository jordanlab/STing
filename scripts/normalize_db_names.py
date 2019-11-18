#!/usr/bin/env python3

import argparse
import os
import re
import sys
from get_sting_command import normalizeName, getIndexPath

def setup_argument_parser():
    parser = argparse.ArgumentParser(description="""This script produces 
        normalizaed names of dbs and schemes, and paths to the corresponding 
        STing index.""")
    parser.add_argument('-i', '--input_file', required=1, help='A tab separated text file with the pairs organism-scheme')
    parser.add_argument('-d', '--dbs_path', required=1, help='Path to STing main database dir location')
    return parser

# def normalizeName(name):
#     normName = name.lower()
#     normName = re.sub(r' |#|-', '_', normName)
#     normName = re.sub(r'([^0-9])\.', r'\1', normName)
#     normName = re.sub(r'/','-', normName)
#     normName = re.sub(r'\(|\)','', normName)
#     normName = re.sub(r'_{2,}','_', normName)
#     return normName

# def getIndexPath(dbsPath, normOrganism, normScheme):
#     indexPath = os.path.join(dbsPath, normScheme, normOrganism, 'db', 'index')
#     return indexPath

def main():
    parser = setup_argument_parser()
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    
    inputFile = args.input_file
    dbsPath   = args.dbs_path
    inFile    = open(inputFile, 'r')
    
    print('n\tdisplay_name\tscheme\toriginal_name\tnormalized_name\tnorm_scheme\tindex_path')
    with inFile:
        inFile.readline()   # Skip header
        for line in inFile:
            fields     = line.rstrip().split('\t')
            name       = fields[1].rstrip()
            scheme     = fields[2].rstrip()
            normName   = normalizeName(name)
            normScheme = normalizeName(scheme)
            indexPath  = getIndexPath(dbsPath, normName, normScheme)
            fields.extend([normName, normScheme, indexPath])
            print('{}'.format('\t'.join(fields)))

if __name__ == '__main__':
    main()
