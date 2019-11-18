#!/usr/bin/env python3

# ##############################################################################
#  Author:   Hector Fabio Espitia-Navarro
#            Georgia Institute of Technology
#  Version:  1.2
#  Date:     04/21/2019
# ##############################################################################

# ToDo: Check if binaries required exist
# ToDo: Check for return status from requests and display error messages
# ToDo: Implement retry with tenacity

import argparse
import csv
import json
import os
import pprint as pp
import re
import requests
import shutil
import subprocess
import sys
import tarfile
import tempfile
import textwrap
import urllib.request
from collections import defaultdict
from tenacity import *
from xml.etree import ElementTree

# Global variables
failed_resources = defaultdict(list)

def setup_argument_parser():
    class CustomHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawTextHelpFormatter,
        argparse.RawDescriptionHelpFormatter):
        pass

    example_text = textwrap.dedent('''\
        example:  
            ./update_dbs.py -f organism-typing_schemes.tsv -o /path/to/sting/main/db/dir -e /path/to/edirect/binaries/dir/ -b path/to/sting/binaries/dir
            ''')
    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''
            This script download all the databases for the MLST and cgMLST typing schemes, and the AMR gene detection, 
            given a file with the corresponding organism dbs and schemes.
            
            Requirements: 
              - Python tenacity module (https://github.com/jd/tenacity)
              - NCBI Entrez Direct command line tools (available at https://www.ncbi.nlm.nih.gov/books/NBK179288/, or via anaconda's bioconda channel)'''),
        epilog=example_text,
        formatter_class=CustomHelpFormatter)
    
    parser.add_argument('--version', action='version', 
        version='%(prog)s\n  version:\t{}\n  last update:\t{}'.format(
                '1.2', 
                '04/21/2019'))
    
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    data_file_help_text = textwrap.dedent('''\
        Tab separated file with 4 columns: 
        - id: id of the organism database
        - display_name: organism db to be shown in GUI 
        - scheme: name of the scheme (e.g., MLST, cgMLST, AMR)
        - original_name: original organism name in the origin website (e.g. in PubMLST)
        - resource_uri: URI of the db''')
    required.add_argument('-f', '--data_file', required=1, help=data_file_help_text)
    required.add_argument('-o', '--out_dir', required=0, help = "Output directory for saving retrieved DBs.", default='pubmlst_dbs', )
    required.add_argument('-e', '--edirect_path', required=1, help='Path to the EDirect (NCBI Entrez Direct) binaries location.')
    required.add_argument('-s', '--sting_path', required=1, help='Path to the STing binaries location')
    # required.add_argument('-d', '--dbs_path', required=1, help='Path to the STing main database directory location')

    parser._action_groups.append(optional)
    return parser


def load_db_file(filename):
    db_file = open(filename, 'r')
    db_dict = defaultdict(list)
    with db_file:
        for line in db_file:
            row = line.rstrip().split('\t')
            id = row[0]
            if not id.isnumeric():  # Skip header
                continue
            display_name  = row[1]
            scheme        = row[2]
            original_name = row[3]
            resource_uri  = row[4]
            db_dict[scheme].append({
                'id'            : id,
                'display_name'  : display_name,
                'original_name' : original_name,
                'resource_uri'  : resource_uri
            })
    return db_dict


def normalize_name(name):
    norm_name = name.lower()
    norm_name = re.sub(r' |#|-', '_', norm_name)
    norm_name = re.sub(r'([^0-9])\.', r'\1', norm_name)
    norm_name = re.sub(r'/', '-', norm_name)
    norm_name = re.sub(r'\(|\)', '', norm_name)
    norm_name = re.sub(r'_{2,}', '_', norm_name)
    return norm_name


def get_sample_name(input_files):
    illumina_regex = re.compile(
        "([A-Za-z0-9-_.#]+)_S[0-9]+_L[0-9]+_R([1-2])_[0-9]+\.(fastq|fq)(\.gz)*"
    )  # https://regex101.com/r/3aEF1g/7
    normal_regex   = re.compile(
        "([A-Za-z0-9-_.#]+)_R*([0-9]{1})\.(fastq|fq)(\.gz)*"
    )  # https://regex101.com/r/3aEF1g/4
    in_file1       = input_files.split(',')[0].rstrip()

    illumina_match_obj = illumina_regex.fullmatch(in_file1)
    normal_match_obj   = normal_regex.fullmatch(in_file1)
    sample_name        = ''

    if illumina_match_obj:
        sample_name = illumina_match_obj.group(1)
    elif normal_match_obj:
        sample_name = normal_match_obj.group(1)

    return sample_name


def get_uri_json_content(reource_uri):
    response = urllib.request.urlopen(reource_uri)
    data     = response.read()      # a `bytes` object
    content  = json.loads(data)
    return content


def fetch_mlst_db(species, out_db_dir):
    sp_name   = species.text.strip("\n")
    n_records = sum(1 for _ in species.iterfind("./mlst/database/loci/locus/url"))
    
    # Fetch files and make config file
    if n_records > 0 :
        if not os.path.exists(out_db_dir):
            os.makedirs(out_db_dir)
        
        # Open config file for writing
        config_file_name = os.path.abspath(os.path.join(out_db_dir, "config.txt"))
        config_file      = open(config_file_name, 'w')
        
        print("Database: \"{}\"".format(sp_name))
        config_file.write("[loci]\n")
        print(" Fetching allele sequences: ")
        
        # Fetch loci sequence files
        for locus_url in species.iterfind("./mlst/database/loci/locus/url") :
            locus_name = os.path.splitext(os.path.basename(locus_url.text))[0]
            file_name  = os.path.join(out_db_dir,locus_name+".fa")
            urllib.request.urlretrieve(locus_url.text, file_name)
            print(" - {} -> {}".format(locus_url.text, file_name))
            config_file.write("{}\t{}\n".format(locus_name, os.path.basename(file_name)))
        
        config_file.write("\n[profile]\n")
        
        # Fetch profile table
        for profile_url in species.iterfind("./mlst/database/profiles/url") :
            file_name = os.path.join(out_db_dir, "profile.txt")
            print(" Fetching profiles: ")
            print(" - {} -> {}".format(profile_url.text, file_name, ""))
            urllib.request.urlretrieve(profile_url.text, file_name)
            config_file.write("{}\t{}\n".format(os.path.basename(out_db_dir), os.path.basename(file_name)))
        
        config_file.close()


def build_index(config_file_name, out_prefix, sting_bin_dir, db_detection_mode = False):
    # (parent, child) = os.path.split(os.path.abspath(os.path.dirname(sys.argv[0])))
    indexer_path    = os.path.join(sting_bin_dir, "indexer")
   
    db_mode_options = ""
    if db_detection_mode:
        db_mode_options = "-m GDETECT"
        
    command_string = "{} -c {} -p {} {}".format(indexer_path, config_file_name,
                                                out_prefix, db_mode_options)
    print(command_string)
    
    output_string = subprocess.check_output(command_string, shell=True, 
                                            stderr=subprocess.STDOUT).decode("utf-8")
    print(output_string)


def rollback_db(out_db_dir):
    if os.path.exists(out_db_dir):
        try:
            shutil.rmtree(out_db_dir)
            return 0
        except OSError as e:
            print ("  Error while trying to remove the directory \"{}\": {} - {}.".format(out_db_dir, e.filename, e.strerror))
            return 1


def download_mlst_dbs(resource_uri, out_dir, sting_bin_dir):
    response = urllib.request.urlopen(resource_uri)
    data     = response.read()      # a `bytes` object
    content  = data.decode('utf-8')
    root     = ElementTree.fromstring(content)
    
    for species in root.iter("species"):
        sp_name    = species.text.strip("\n").rstrip()
        out_db_dir = os.path.join(out_dir, normalize_name(sp_name))
        
        fetch_mlst_db(species, out_db_dir)
        
        print("Building STing index:")
        out_prefix       = os.path.abspath(os.path.join(out_db_dir, 'db', "index"))
        config_file_name = os.path.abspath(os.path.join(out_db_dir, "config.txt"))
        build_index(config_file_name, out_prefix, sting_bin_dir)


def fecth_cgmlst_db(resource_uri, out_db_dir):
    scheme       = get_uri_json_content(resource_uri)
    profiles_uri = scheme['profiles_csv']
    loci_uri     = scheme['loci']
    
    if not os.path.exists(out_db_dir):
        os.makedirs(out_db_dir)
   
    # Open config file for writing
    config_file_name = os.path.abspath(os.path.join(out_db_dir, "config.txt"))
    
    with open(config_file_name, 'w') as config_file:
        config_file.write("[loci]\n")
        print(" Fetching allele sequences: ")
        
        for locus_uri in loci_uri:
            locus       = get_uri_json_content(locus_uri)
            locus_name  = locus['id']
            alleles_uri = locus['alleles_fasta'] 
            file_name   = os.path.join(out_db_dir,locus_name+".fa")
            try:
                response = urllib.request.urlopen(locus_uri)
                with open(file_name, 'w') as locus_file:
                    locus_file.write(response.read().decode("utf-8") )
                    print(" - {} -> {}".format(alleles_uri, file_name))
                    config_file.write("{}\t{}\n".format(locus_name, os.path.basename(file_name)))
                    
            except urllib.error.HTTPError as e:
                res = json.loads(e.read())
                print("  Alleles not defined for this locus of this database:\n   Message: {}".format(res['message']))
                print("  Removing downloaded files of this database... ", end = '')
                res = rollback_db(out_db_dir)
                if res == 0:
                    print("Done!")
                else:
                    print("Manual removing of \"{}\" directory is recommended.")
            
        config_file.write("\n[profile]\n")
        # Fetch profile table
        file_name = os.path.join(out_db_dir, "profile.txt")
        print(" Fetching profiles: ")
        try:
            response = urllib.request.urlopen(profiles_uri)
            with open(file_name, 'w') as profile_file:
                profile_file.write(response.read().decode("utf-8") )
                print(" - {} -> {}".format(profiles_uri, file_name))
                config_file.write("{}\t{}\n".format(os.path.basename(out_db_dir), os.path.basename(file_name)))
                
        except urllib.error.HTTPError as e:
            res = json.loads(e.read())
            print("  Profiles not defined for this database:\n   Message: {}".format(res['message']))
            print("  Removing downloaded files of this database... ", end = '')
            res = rollback_db(out_db_dir)
            
            if res == 0:
                print("Done!")
            else: 
                print("Manual removing of \"{}\" directory is recommended.")


def download_cgmlst_dbs(db_list, out_dir, sting_bin_dir):
    # pp.pprint(db_list)
    for db in db_list:
        sp_name    = db['original_name'].rstrip()
        scheme_uri = db['resource_uri'].rstrip()
        out_db_dir  = os.path.join(out_dir, normalize_name(sp_name))
        # print(scheme_uri)
        scheme = get_uri_json_content(scheme_uri)
        if ('profiles_csv' in scheme.keys()):
            print("Database: \"{}\"".format(sp_name)) 
            fecth_cgmlst_db(scheme_uri, out_db_dir)
            
            print("Building STing index:")
            out_prefix       = os.path.abspath(os.path.join(out_db_dir, 'db', "index"))
            config_file_name = os.path.abspath(os.path.join(out_db_dir, "config.txt"))
            build_index(config_file_name, out_prefix, sting_bin_dir)


def replace_non_standard_chars(file_name):
    file_content = ""
    with open(file_name, 'r') as file_h:
        for line in file_h:
            if not line.startswith('>'):
                line = re.sub(r'[^ACTGNactgn\n]','N', line)
            file_content += line
    
    return file_content


def fetch_and_process_file(db_name, db_uri, edirect_path):
    file_content = ""
    with tempfile.TemporaryDirectory(prefix = 'sting_tmp') as tmp_dir:
        temp_filename = ""
        
        if db_name == "CARD":
            response           = urllib.request.urlopen(db_uri)
            card_nucl_filename = "nucleotide_fasta_protein_homolog_model.fasta"
            temp_filename      = os.path.join(tmp_dir, "card.tar.bz2")
            
            with open(temp_filename, 'wb') as temp_file:
                temp_file.write(response.read())
                
            with tarfile.open(temp_filename, "r:bz2") as tar_file:              
                tar_file.extractall(tmp_dir)
            
            target_filename = os.path.join(tmp_dir, card_nucl_filename)
            if (os.path.isfile(target_filename)):
                file_content = replace_non_standard_chars(target_filename)
        
        elif db_name == "NCBI-AMR":
            esearch_path = os.path.join(edirect_path, "esearch")
            efetch_path  = os.path.join(edirect_path, "efetch")
            query_str    = db_uri
            # query_str   = "NG_047055.1"
            edirect_cmd = "{} -db nucleotide -query \"{}\" | {} -db nucleotide -format fasta".format(esearch_path, query_str, efetch_path)
            target_filename = os.path.join(tmp_dir, "ndaro.fasta")
            # target_filename = "ndaro.fasta"
            with open(target_filename, 'wb') as target_file:
                print("  {}".format(edirect_cmd))
                out = subprocess.check_output(edirect_cmd, shell=True)
                target_file.write(out)
            
            if (os.path.isfile(target_filename)):
                file_content = replace_non_standard_chars(target_filename)
                
        else:
            response     = urllib.request.urlopen(db_uri)
            temp_filename = os.path.join(tmp_dir, "genes.fasta")
            with open(temp_filename, 'wb') as temp_file:
                temp_file.write(response.read())
                
            if (os.path.isfile(temp_filename)):
                file_content = replace_non_standard_chars(temp_filename)
            
    return file_content


def fetch_amr_db(db_dict, out_dir, edirect_path):
    db_name    = db_dict['original_name'].rstrip()
    db_uri     = db_dict['resource_uri'].rstrip()
    out_db_dir = os.path.join(out_dir, normalize_name(db_name))
    
    if not os.path.exists(out_db_dir):
        os.makedirs(out_db_dir)
   
    # Open config file for writing
    config_file_name = os.path.abspath(os.path.join(out_db_dir, "config.txt"))
    
    with open(config_file_name, 'w') as config_file:
        config_file.write("[loci]\n")
        print(" Fetching gene sequences: ")
        
        locus_name = "genes"
        file_name  = os.path.join(out_db_dir,locus_name+".fa")
        try:
            file_content = fetch_and_process_file(db_name, db_uri, edirect_path)
            with open(file_name, 'w') as locus_file:
                locus_file.write(file_content)
                print(" - {} -> {}".format(db_uri, file_name))
                config_file.write("{}\t{}\n".format(normalize_name(db_name), os.path.basename(file_name)))
                
        except urllib.error.HTTPError as e:
            res = json.loads(e.read())
            print("  Sequences are not available for this database:\n   Message: {}".format(res['message']))
            print("  Removing downloaded files of this database... ", end = '')
            res = rollback_db(out_db_dir)
            if res == 0:
                print("Done!")
            else:
                print("Manual removing of \"{}\" directory is recommended.")


def download_amr_dbs(db_list, out_dir, edirect_path, sting_bin_dir):
    # pp.pprint(db_list)
    for db in db_list:
        db_name    = db['original_name'].rstrip()
        db_uri     = db['resource_uri'].rstrip()
        out_db_dir = os.path.join(out_dir, normalize_name(db_name))
        
        print("Database: \"{}\"".format(db_name))
        fetch_amr_db(db, out_dir, edirect_path)
                   
        print("Building STing index:")
        out_prefix       = os.path.abspath(os.path.join(out_db_dir, 'db', "index"))
        config_file_name = os.path.abspath(os.path.join(out_db_dir, "config.txt"))
        build_index(config_file_name, out_prefix, sting_bin_dir, db_detection_mode = True)


def main():
    parser = setup_argument_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    data_file     = args.data_file
    out_dir       = args.out_dir
    edirect_path  = args.edirect_path
    sting_bin_dir = args.sting_path

    out_dir = os.path.join(os.getcwd(), args.out_dir)
        
    db_dict = load_db_file(data_file)
    # pp.pprint(db_dict)
    
    # Download all the MLST dbs from PubMLST
    # All the MLST dbs have the same URI, so getting the one from the first element
    mlst_resource_uri = db_dict['MLST'][1]['resource_uri']
    mlst_out_db_dir = os.path.join(out_dir, 'mlst')
    # download_mlst_dbs(mlst_resource_uri, mlst_out_db_dir, sting_bin_dir)
    
    # Download the AMR dbs
    amr_out_db_dir = os.path.join(out_dir, 'amr')
    download_amr_dbs(db_dict['AMR'], amr_out_db_dir, edirect_path, sting_bin_dir)
    
    # Download the cgMLST dbs from PubMLST
    cgmlst_out_db_dir = os.path.join(out_dir, 'cgmlst')
    # download_cgmlst_dbs(db_dict['cgMLST'], cgmlst_out_db_dir, sting_bin_dir)


if __name__ == '__main__':
    main()
