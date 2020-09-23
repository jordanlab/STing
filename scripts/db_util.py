#!/usr/bin/env python3

# ##############################################################################
#  Author:   Hector Fabio Espitia-Navarro
#            Georgia Institute of Technology
#  Version:  1.1
#  Date:     04/15/2019
# ##############################################################################

import argparse
import os
import re
import sys
import subprocess
import urllib.request
from xml.etree import ElementTree

# ==============================================================================
class Error(Exception):
   """Base class for other exceptions"""
   pass
class BinaryNotFoundError(Error):
   """Raised when the input value is too small"""
   pass

# ==============================================================================
def setup_argument_parser():
    class CustomHelpFormatter(
                argparse.ArgumentDefaultsHelpFormatter,
                argparse.RawTextHelpFormatter,
                argparse.RawDescriptionHelpFormatter):
            pass
            
    parser = argparse.ArgumentParser(description="This script provides a set of utilities to download databases from PubMLST and build STing indices from them.",
        formatter_class=CustomHelpFormatter)
    subparsers = parser.add_subparsers(title="subcommands",
                                       # description="Valid subcommands:",
                                       # help='valid subcommands',
                                       dest="subparser_name")
    
    # list task parser
    parserList = subparsers.add_parser('list', 
                                        help = 'List all the available databases at PubMLST')
    parserList.set_defaults(func = list_pubmlst_dbs)
    
    # query task parser
    parserQuery = subparsers.add_parser('query', 
                                        help = 'Search a database in PubMLST')
    parserQuery.set_defaults(func = search_db)
    parserQuery.add_argument('query', help = "Query term")
    
    # fetch task parser
    parserFetch = subparsers.add_parser('fetch', 
                                        help = 'Fetch a database from PubMLST')
    parserFetch.set_defaults(func = fetch)
    parserFetch.add_argument('-q', '--query', help = "Query term")
    parserFetch.add_argument('-a', '--all', 
                             action ='store_true', 
                             help   = "Fetch all databases available in PubMLST")
    parserFetch.add_argument('-b', '--build_index', 
                             action ='store_true', 
                             help   = "Build a STing MLST index from DB fetched")
    parserFetch.add_argument('-o', '--out_dir', 
                             default ='pubmlst_dbs',
                             help    = "Output directory for saving retrieved DBs. [Default: pubmlst_dbs]")
    
    parser.add_argument('--version', action='version', version='%(prog)s\n  version:\t{}\n  last update:\t{}'.format(
                '1.2.1', 
                '07/09/2020'))
    return parser

# ==============================================================================
def get_db_file_content(fileUrl):
    response = urllib.request.urlopen(fileUrl)
    data     = response.read()      # a `bytes` object
    content  = data.decode('utf-8')
    return content

# ==============================================================================
def list_pubmlst_dbs(fileContent):
    root    = ElementTree.fromstring(fileContent)
    counter = 0
    
    print("{}\t{}\t{}\t{}\t{}".format("#", "Database", "#Profiles", "Retrieved", "DB_URL"))
    # for species in root.findall('species'):
    for species in root.iter('species'):
        counter  += 1
        spName    = species.text.strip("\n")
        nProfiles = 0
        retDate   = ""
        dbUrl     = ""
        for count in species.iterfind("./mlst/database/profiles/count"):
            nProfiles = count.text
        for retrieved in species.iterfind("./mlst/database/retrieved"):
            retDate = retrieved.text
        for url in species.iterfind("./mlst/database/url"):
            dbUrl = url.text
        
        print("{}\t{}\t{}\t{}\t{}".format(counter, spName, nProfiles, retDate, dbUrl))

# ==============================================================================
def search_db(fileContent, query):
    # [tag='text']
    root    = ElementTree.fromstring(fileContent)
    counter = 0
    for species in root.iter("species"):
        spName    = species.text.strip("\n")
        nProfiles = 0
        retDate   = ""
        dbUrl     = ""
        
        if spName.lower().find(query.lower()) != -1 :
            counter += 1
            for count in species.iterfind("./mlst/database/profiles/count"):
                nProfiles = count.text
            for retrieved in species.iterfind("./mlst/database/retrieved"):
                retDate = retrieved.text
            for url in species.iterfind("./mlst/database/url"):
                dbUrl = url.text
        
            print("{}\t{}\t{}\t{}\t{}".format(counter, spName,nProfiles, retDate, dbUrl))

# ==============================================================================
def getOutDbDir(dbName):
    normName = dbName.lower()
    normName = re.sub(r' |#|-', '_', normName)
    normName = re.sub(r'([^0-9])\.', r'\1', normName)
    normName = re.sub(r'/','-', normName)
    normName = re.sub(r'\(|\)','', normName)
    normName = re.sub(r'_{2,}','_', normName)
    return normName

# ==============================================================================
def fetch_db(species, outDbDir):
    spName   = species.text.strip("\n")
    nRecords = sum(1 for _ in species.iterfind("./mlst/database/loci/locus/url"))
    
    # Fetch files and make config file
    if nRecords > 0 :
        if not os.path.exists(outDbDir):
            os.makedirs(outDbDir)
        
        # Open config file for writing
        configFileName = os.path.abspath(os.path.join(outDbDir, "config.txt"))
        configFile     = open(configFileName, 'w')
        
        print("Database: \"{}\"".format(spName))
        
        configFile.write("[loci]\n")
        
        print(" Fetching allele sequences: ")
        
        # Fetch loci sequence files
        for locus in species.iterfind("./mlst/database/loci/locus"):
            locusName = locus.text.strip()
            locusUrl = locus.find("./url")
            fileName  = os.path.join(outDbDir,locusName+".fa")
            urllib.request.urlretrieve(locusUrl.text, fileName)
            print(" - {} -> {}".format(locusUrl.text, fileName))
            configFile.write("{}\t{}.fa\n".format(locusName, locusName))
        
        configFile.write("\n[profile]\n")
        
        # Fetch profile table
        for profileUrl in species.iterfind("./mlst/database/profiles/url") :
            fileName = os.path.join(outDbDir, "profile.txt")
            print(" Fetching profiles: ")
            urllib.request.urlretrieve(profileUrl.text, fileName)
            print(" - {} -> {}".format(profileUrl.text, fileName, ""))
            configFile.write("{}\tprofile.txt\n".format(os.path.basename(outDbDir)))
        
        configFile.close()
        
# ==============================================================================
def fetch(fileContent, args, parser):
    root = ElementTree.fromstring(fileContent)
    
    if args.all or (args.query != None) :
        outDir = os.path.join(os.getcwd(), args.out_dir)
            
        for species in root.iter("species"):
            spName   = species.text.strip("\n")
            outDbDir = os.path.join(outDir, getOutDbDir(spName))
            if args.query != None:
                if spName.lower().find(args.query.lower()) == -1 :
                    continue
            
            fetch_db(species, outDbDir)
            
            if args.build_index:
                print("Building STing index:")
                outPrefix = os.path.abspath(os.path.join(outDbDir, 'db'))
                if not os.path.exists(outPrefix):
                    os.makedirs(outPrefix)
                outPrefix = os.path.join(outPrefix, 'index')
                configFileName = os.path.abspath(os.path.join(outDbDir, "config.txt"))
                try:
                    build_index(configFileName, outPrefix)
                except BinaryNotFoundError:
                    print("{} ERROR - fetch: STing indexer is not installed. The STing index was not built.".format(os.path.basename(sys.argv[0])))
    else:
        print("{} ERROR - fetch: One of --all or --query options are required".format(os.path.basename(sys.argv[0])))
        
        subparsers_actions = [
            action for action in parser._actions 
            if isinstance(action, argparse._SubParsersAction)]
        # there will probably only be one subparser_action,
        # but better save than sorry
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                if choice == "fetch":
                    print(subparser.format_help())

# ==============================================================================
def which(program):
    ''' Look for a binary and return its path.
    From https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    '''
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
    
# ==============================================================================
def build_index(configFileName, outPrefix):
    indexerPath = which("indexer")
    if not indexerPath:
        raise BinaryNotFoundError
    else:
        commandString = "{} -c {} -p {}".format(indexerPath, configFileName,
                                                  outPrefix)
        print(commandString)
        outputString = subprocess.check_output(commandString, shell = True, 
                                               stderr=subprocess.STDOUT).decode("utf-8")
        print(outputString)

# ==============================================================================
def main():
    parser = setup_argument_parser()
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    
    # URL to PubMLST databases file
    dbFileUrl     = "https://pubmlst.org/data/dbases.xml"
    dbFileContent = get_db_file_content(dbFileUrl)
    
    if args.subparser_name == "list" :
        args.func(dbFileContent)
    elif args.subparser_name == "query" :
        args.func(dbFileContent, args.query)
    elif args.subparser_name == "fetch" :
        args.func(dbFileContent, args, parser)

# ==============================================================================
if __name__ == '__main__':
    main()
