// ############################################################################
// STing
// Copyright (C) 2019 Georgia Institute of Technology
// STing is freely available to personal, academic and non-profit use only. You
// cannot redistribute STing to other users. No liability for software usage is
// assumed. For more information on how to obtain a commercial license please 
// conact Lavanya Rishishwar <lavanya.rishishwar@gatech.edu>
// ============================================================================ 
// (c) Georgia Institute of Technology 2019
// Author:  Hector F. Espitia-Navarro
//          hspitia@gatech.edu
//          School of Biological Sciences
//          Georgia Institute of Technology
//          
//          See AUTHORS file for more information
// 
// Contact: Lavanya Rishishwar
//          lavanya.rishishwar@gatech.edu
//          School of Biological Sciences
//          Georgia Institute of Technology
// ============================================================================ 
// Patent information
// Espitia, H., Chande, A. T., Jordan, I. K., & Rishishwar, L. (2017). 
//      A method of sequence typing with in silico aptamers from a next
//      generation sequencing platform. Google Patents. Patent application 
//      US15/726,005. Retrieved from 
//      https://patents.google.com/patent/US20190108308A1
// ############################################################################


#ifndef INDEXER_H
#define INDEXER_H

// ============================================================================
// Includes
// ============================================================================
// std
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
// ----------------------------------------------------------------------------
// c
// ----------------------------------------------------------------------------
#include <limits.h>
// ----------------------------------------------------------------------------
// SeqAn
// ----------------------------------------------------------------------------
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <seqan/stream.h>
// ----------------------------------------------------------------------------
// Application
// ----------------------------------------------------------------------------
#include "IndexerOptions.h"
#include "types.hpp"
#include "data_structures.hpp"
#include "dir_utils.hpp"
#include "custom_utils.hpp"
// ============================================================================
// Namespaces
// ============================================================================
using namespace seqan;
using namespace std;

class Indexer
{
        // memebers
        TConfigMap        config_map;
        TCharStringSet    ids;
        TUintVector       loci_alleles_index;
        TLociTable        loci_table;
        IndexerOptions    options;
        ArgumentParser    parser;
        TCharEsaIndex     profiles_esa_index;
        ProfilesTable     profiles_table;
        TSequenceEsaIndex alleles_esa_index;
        string            profile_name;
        string            program_name;
        TStringSet        seqs;
    
    public:
        // constructors / destructors
        Indexer();
        Indexer(const IndexerOptions & options);
        virtual ~Indexer();
        int run(int argc, char const ** argv);
    
    protected:
        // Class methods (different form constructors, accessors)
        void setupArgumentParser();
        ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv);
        void processMlstConfigFile();
        void processGeneDetectorConfigFile();
        void loadLociFiles();
        // void loadProfilesFile(string const &  filename);
        void loadProfilesFile();
        // void loadProfilesFile2(string const &  filename);
        void getProfileStrings(TCharStringSet & profiles);
        SaveFileResult createProfilesIndex(const TCharStringSet & profiles, 
                                           const string &         output_file_name);
        SaveFileResult createAllelesIndex(const string & output_file_name);
        SaveFileResult saveLociTable(const string & output_file_name);
        SaveFileResult saveLociAllelesIndex(const string & output_file_name);
        SaveFileResult saveProfilesFile(const string & output_file_name);
        SaveFileResult saveAlleleSequenceIds(const string & output_file_name);
        uint64_t getSeqsSize();
        int runMlstIndexer();
        int runGeneDetectorIndexer();
        
        template<typename T, typename U>
        SaveFileResult createIndex(T &            esa_index, 
                                   const U &      text, 
                                   const string & output_file_name);
        SaveFileResult copyFile(const string & source, const string & dest);
        
        // Accessors
        const TSequenceEsaIndex& getAllelesEsaIndex() const;
        void setAllelesEsaIndex(const TSequenceEsaIndex & alleles_esa_index);
        const TConfigMap & getConfigMap() const;
        void setConfigMap(const TConfigMap & config_map);
        const TCharStringSet & getIds() const;
        void setIds(const TCharStringSet & ids);
        TUintVector getLociAllelesIndex() const;
        void setLociAllelesIndex(TUintVector loci_alleles_index);
        const TLociTable & getLociTable() const;
        void setLociTable(const TLociTable & loci_table);
        const IndexerOptions & getOptions() const;
        void setOptions(const IndexerOptions & options);
        const ArgumentParser & getParser() const;
        void setParser(const ArgumentParser & parser);
        const string & getProfileName() const;
        void setProfileName(const string & profile_name);
        const TCharEsaIndex& getProfilesEsaIndex() const;
        void setProfilesEsaIndex(const TCharEsaIndex & profiles_esa_index);
        const ProfilesTable& getProfilesTable() const;
        void setProfilesTable(const ProfilesTable & profiles_table);
        const string & getProgramName() const;
        void setProgramName(const string & program_name);
        const TStringSet & getSeqs() const;
        void setSeqs(const TStringSet & seqs);
};

#endif /* INDEXER_H */
