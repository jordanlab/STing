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


#ifndef KMER_MATCHER_OPTIONS_H
#define KMER_MATCHER_OPTIONS_H

// ============================================================================
// Includes
// ============================================================================
// std
// ----------------------------------------------------------------------------
#include <string>
// #include <unordered_map>
#include <vector>
// ----------------------------------------------------------------------------
// SeqAn (https://www.seqan.de/)
// ----------------------------------------------------------------------------
#include <seqan/index.h>
// // ----------------------------------------------------------------------------
// // Binary Search Tree (https://gist.github.com/mgechev/5911348 with modifications)
// // ----------------------------------------------------------------------------
// #include <BST/BST.h>
// ----------------------------------------------------------------------------
// Application
// ----------------------------------------------------------------------------
#include "types.hpp"
// ============================================================================
// Namespaces
// ============================================================================
using namespace seqan;
using namespace std;

class KmerMatcherOptions
{
        string        index_filename;
        TStringVector in_fastq_1_files;
        TStringVector in_fastq_2_files;
        string        sample_name;
        uint64_t      kmer_length;
        uint64_t      verbosity_level;
        bool          verbose;
        bool          print_time;
        bool          print_kmer_counts;
        string        kmer_depth_out_file;
        string        ext_kmer_depth_out_file;
        string        output_file;
        bool          tidy_results;
        
        // typer specific
        bool          print_allele_cov;
        bool          print_kmer_depth;
        bool          fast_mode;
        bool          sensitive_mode;
        uint64_t      n_top_alleles;
        // uint64_t      depth_cutoff;
        
        // detector specific
        bool          vertical_results;
        bool          print_counts;
        bool          print_kmer_perc;
        double        threshold_perc;
        

    public:
        KmerMatcherOptions();
        virtual ~KmerMatcherOptions();
        string getIndexFilename() const;
        void setIndexFilename(const string & value);
        TStringVector getInFastq1Files() const;
        void setInFastq1Files(const TStringVector & value);
        TStringVector getInFastq2Files() const;
        void setInFastq2Files(const TStringVector & value);
        string getSampleName() const;
        void setSampleName(const string & value);
        uint64_t getKmerLength() const;
        void setKmerLength(const uint64_t & value);
        uint64_t getVerbosityLevel() const;
        uint64_t getNTopAlleles() const;
        void setNTopAlleles(const uint64_t & value);
        void setVerbosityLevel(const uint64_t & value);
        bool getVerticalResults() const;
        void setVerticalResults(const bool & value);
        bool getTidyResults() const;
        void setTidyResults(const bool & value);
        bool getVerbose() const;
        void setVerbose(const bool & value);
        bool getPrintTime() const;
        void setPrintTime(const bool & value);
        bool getPrintKmerCounts() const;
        void setPrintKmerCounts(const bool & value);
        bool getPrintCounts() const;
        void setPrintCounts(const bool & value);
        bool getPrintKmerPerc() const;
        void setPrintKmerPerc(const bool & value);
        bool getPrintAlleleCov() const;
        void setPrintAlleleCov(const bool & value);
        bool getPrintKmerDepth() const;
        void setPrintKmerDepth(const bool & value);
        bool getFastMode() const;
        void setFastMode(const bool & value);
        bool getSensitiveMode() const;
        void setSensitiveMode(const bool & value);
        double getThresholdPerc() const;
        void setThresholdPerc(double & value);
        string getKmerDepthOutFile() const;
        void setKmerDepthOutFile(const string & value);
        string getExtKmerDepthOutFile() const;
        void setExtKmerDepthOutFile(const string & value);
        string getOutputFile() const;
        void setOutputFile(const string & value);
        // double getDepthCutoff() const;
        // void setDepthCutoff(uint64_t & value);
};

#endif // KMER_MATCHER_OPTIONS_H
