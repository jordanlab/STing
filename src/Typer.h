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


#ifndef TYPER_H
#define TYPER_H

// ============================================================================
// Includes
// ============================================================================
// ----------------------------------------------------------------------------
// Application
// ----------------------------------------------------------------------------
#include "Timer.h"
#include "SequenceInfo.h"
#include "KmerMatcher.h"
// ============================================================================
// Namespaces
// ============================================================================
using namespace seqan;
using namespace std;

class Typer: public KmerMatcher
{
        // vector<int>    first_max_allele;                // indices of alleles with the top kmer hits
        // TStringVector  first_max_allele_numbers;        // allele numbers with top kmer hits (no ST)
        // TUintVector    first_max_kmer_counts;           // kmer counts for each locus of the inferred allelic profile
        
        // TDoubleVector  first_max_norm_kmer_counts;      // kmer counts for each locus of the inferred allelic profile
        // TStringVector  first_inferred_profile;          // inferred allelic profile with ST
        // string         first_inferred_profile_string;   // string representation to be found in profiles tree
        
        TUMapAlleleCountVector  top_alleles;            // Map with top alleles (SAlleleCount type) of each locus
        TUMapSeqInfoVector      top_alleles_info;
        TUintVector             top_alleles_indices;
        
        // vector<int>             best_max_allele;            // indices of the alleles with top kmer hits 
        TStringVector           best_max_allele_numbers;    // allele numbers with the top kmer hits (no ST) and best coverage and/or sd
        // TDoubleVector           best_max_norm_kmer_counts;
        TUMapSequenceInfo       best_allelic_profile;
        
        
        TCharEsaIndex       profiles_index;
        ProfilesTable       profiles_table;

        string              inferred_st;
        bool                fast_mode_on;
        bool                sensitive_mode_on;
        
        TSequenceEsaIndex   top_alleles_esa_index;
        vector<int>    reduced_index_allele_indices;
        
        // bool calledLocalProcessRead;

    public:
        enum StStatus
        {
            ST_PREDICTED   = 0,
            NO_KMER_HITS   = 1,
            NO_ST_IN_TABLE = 2
        };
        
        // constructors / destructors
        Typer();
        Typer(const KmerMatcherOptions & options);
        virtual ~Typer();
        int run(int argc, char const ** argv);
        
        StStatus status;
        
        // Class methods (different form constructors, accessors)
        void setupArgumentParser(const string & program_name,
                                 const string & program_description);
        ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv);
        void setupReducedIndexFromTopAlleles();
        void processReadsFiles(bool save_kmer_starts);
        // template <typename T>
        // void reprocessKmerizedReads(T &                       in_fastq_file,
        //                             uint64_t                  fastq_file_idx,
        //                             uint64_t                  k,
        //                             TSequenceEsaIndexFinder & finder_esa,
        //                             bool save_kmer_starts);
        KmerMatcher::ProcessReadResult 
        processRead(TSequence const &         read, 
                    uint64_t                  k, 
                    TSequenceEsaIndexFinder & finder_esa,
                    bool                      save_kmer_starts);
        void selectBestAlleles();
        void getTopNormAlleles(uint64_t n_top); // **
        // void computeKmerDepthAndAlleleCoverageNew();
        void computeKmerDepthAndAlleleCoverageForTopAlleles();
        void loadProfilesEsaIndex();
        void loadProfilesFile();
        void selectBestAllelesBasedOnCoverage();
        // StStatus getInferredSt();
        StStatus getInferredSt(string  &      inferred_st,
                               const string & profile_string);
        StStatus getInferredSt(TStringVector & first_inferred_profile,
                               const string & profile_string);
        SaveFileResult saveTopAllelesPerBaseKmerDepthData(const string & output_file_name);
        void updateBestAllelicProfileCovData();
        void printAllelicProfile(const TStringVector & inferred_profile);
        template<typename T>
        void printResultsLine(const string &    line_type,
                              const StStatus &  status,
                              const vector<T> & info_vector);
        template<typename T>
        string getResultsLine(const string &    line_type,
                              const StStatus &  status,
                              const vector<T> & info_vector);
        string getResultsHeader();
        string getResultsString(const StStatus & status);
        string getTidyResultsHeader();
        string getTidyResultsString(const StStatus & status);
        void printResults(const StStatus & status);
        void printTidyResults(const StStatus & status);
        
        // Accessors
        const TCharEsaIndex& getProfilesIndex() const;
        void setProfilesIndex(const TCharEsaIndex & profiles_index);
        const ProfilesTable& getProfilesTable() const;
        void setProfilesTable(const ProfilesTable & profiles_table);
    
        // string createProfileStringRep();
        string createProfileStringRep(TStringVector max_allele_numbers);
        string getStatusString(const Typer::StStatus & status);
        string getInfoMessage();
        void determineLowConfidence();
        void determineSensitiveMode();
        bool isCoverageCalculated();
        void printTopAllelesInfo();
        void printBestAllelicProfile();
};

#endif // TYPER_H

