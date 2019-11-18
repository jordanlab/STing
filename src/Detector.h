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


#ifndef DETECTOR_H
#define DETECTOR_H

// ============================================================================
// Includes
// ============================================================================
// ----------------------------------------------------------------------------
// Application
// ----------------------------------------------------------------------------
#include "Timer.h"
#include "KmerMatcher.h"
// ============================================================================
// Namespaces
// ============================================================================
using namespace seqan;
using namespace std;

class Detector: public KmerMatcher
{
    // TUintVector presence;
    vector<bool>   presence;
    vector<bool>   pres;
    vector<double> kmer_matches_perc;
    TUintVector    seqsWithCountsIndices;
    
        
    public:
        // constructors / destructors
        Detector();
        Detector(const KmerMatcherOptions & options);
        virtual ~Detector();
        int run(int argc, char const ** argv);
    private:
        // Class methods (different form constructors, accessors)
        void setupArgumentParser(string const & program_name,
                                 string const & program_description);
        ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv);
        // void processReadsFiles(bool save_kmer_starts);
        // vector<double> calculateMatchesPercentage();
        void calculateMatchesPercentage();
        // vector<bool> calculatePresence();
        void calculatePresence();
        void computeKmerDepthAndAlleleCoverage();
        void computeKmerDepthAndAlleleCoverageNew();
        void printResultsVertical();
        template<typename T>
        void printResultsLine(const string &    line_type, 
                              const vector<T> & info_vector);
        template<typename T>
        string getResultsLine(const string &    line_type, 
                              const vector<T> & info_vector);
        string getResultsHeader(uint64_t & header_length);
        string getResultsString();
        string getTidyResultsHeader(uint64_t & header_length);
        string getTidyResultsString(const uint64_t & header_length);
        void printResultsHorizontal();
        void printTidyResults();
        // void printPresence(TUintVector const & presences, uint64_t n_separators);
        // void printPresence(TBoolVector const & presences, uint64_t n_separators);
        // void printKmerPerc(TDoubleVector const & kmer_percs, uint64_t n_separators);
        // void printCounts();
        // void printAlleleCoverage(uint64_t const & idx, uint64_t n_separators);
        void getSeqsWithCounts();
        // SaveFileResult saveKmerDepthData(const string & output_file_name);
        // SaveFileResult saveKmerDepthData(const string & output_file_name, 
        //                                  TUintVector & gene_kmer_depths,
        //                                  uint64_t const & seq_idx,
        //                                  bool append);
};

#endif // DETECTOR_H
