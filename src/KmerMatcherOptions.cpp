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


#include "KmerMatcherOptions.h"

using namespace seqan;
using namespace std;

// -----------------------------------------------------------------------------------
// Constructors /destructors
// -----------------------------------------------------------------------------------

KmerMatcherOptions::KmerMatcherOptions()
{
    index_filename      = "";
    TStringVector v;
    in_fastq_1_files    = v;
    in_fastq_2_files    = v;
    sample_name         = "";
    kmer_length         = 30;
    verbosity_level     = 0;
    verbose             = false;
    print_kmer_counts   = false;
    kmer_depth_out_file = "";
    tidy_results        = false;
    
    print_allele_cov    = false;
    print_kmer_depth    = false;
    fast_mode           = false;
    
    vertical_results    = false;
    print_counts        = false;
    print_kmer_perc     = false;
    threshold_perc      = false;
}

KmerMatcherOptions::~KmerMatcherOptions()
{

}
// -----------------------------------------------------------------------------------
// Class methods (different from accessors)
// -----------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------
// Getters and setters
// -----------------------------------------------------------------------------------
string KmerMatcherOptions::getIndexFilename() const
{
    return index_filename;
}

void KmerMatcherOptions::setIndexFilename(const string & value)
{
    index_filename = value;
}

TStringVector KmerMatcherOptions::getInFastq1Files() const
{
    return in_fastq_1_files;
}

void KmerMatcherOptions::setInFastq1Files(const TStringVector & value)
{
    in_fastq_1_files = value;
}

TStringVector KmerMatcherOptions::getInFastq2Files() const
{
    return in_fastq_2_files;
}

void KmerMatcherOptions::setInFastq2Files(const TStringVector & value)
{
    in_fastq_2_files = value;
}

string KmerMatcherOptions::getSampleName() const
{
    return sample_name;
}

void KmerMatcherOptions::setSampleName(const string & value)
{
    sample_name = value;
}

uint64_t KmerMatcherOptions::getKmerLength() const
{
    return kmer_length;
}

void KmerMatcherOptions::setKmerLength(const uint64_t & value)
{
    kmer_length = value;
}

uint64_t KmerMatcherOptions::getNTopAlleles() const
{
    return n_top_alleles;
}

void KmerMatcherOptions::setNTopAlleles(const uint64_t & value)
{
    n_top_alleles = value;
}

uint64_t KmerMatcherOptions::getVerbosityLevel() const
{
    return verbosity_level;
}

void KmerMatcherOptions::setVerbosityLevel(const uint64_t & value)
{
    verbosity_level = value;
}

bool KmerMatcherOptions::getVerticalResults() const
{
    return vertical_results;
}

void KmerMatcherOptions::setVerticalResults(const bool & value)
{
    vertical_results = value;
}

bool KmerMatcherOptions::getTidyResults() const
{
    return tidy_results;
}

void KmerMatcherOptions::setTidyResults(const bool & value)
{
    tidy_results = value;
}

bool KmerMatcherOptions::getVerbose() const
{
    return verbose;
}

void KmerMatcherOptions::setVerbose(const bool & value)
{
    verbose = value;
}

bool KmerMatcherOptions::getPrintTime() const
{
    return print_time;
}

void KmerMatcherOptions::setPrintTime(const bool & value)
{
    print_time = value;
}        
        
bool KmerMatcherOptions::getPrintKmerCounts() const
{
    return print_kmer_counts;
}

void KmerMatcherOptions::setPrintKmerCounts(const bool & value)
{
    print_kmer_counts = value;
}

bool KmerMatcherOptions::getPrintKmerPerc() const 
{
    return print_kmer_perc;
}

void KmerMatcherOptions::setPrintKmerPerc(const bool & value)
{
    print_kmer_perc = value;
}

bool KmerMatcherOptions::getPrintCounts() const
{
    return print_counts;
}

void KmerMatcherOptions::setPrintCounts(const bool & value)
{
    print_counts = value;
}

bool KmerMatcherOptions::getPrintAlleleCov() const
{
    return print_allele_cov;
}

void KmerMatcherOptions::setPrintAlleleCov(const bool & value)
{
    print_allele_cov = value;
}

bool KmerMatcherOptions::getPrintKmerDepth() const
{
    return print_kmer_depth;
}

void KmerMatcherOptions::setPrintKmerDepth(const bool & value)
{
    print_kmer_depth = value;
}

bool KmerMatcherOptions::getFastMode() const
{
    return fast_mode;
}

void KmerMatcherOptions::setFastMode(const bool & value)
{
    fast_mode = value;
}

bool KmerMatcherOptions::getSensitiveMode() const
{
    return sensitive_mode;
}

void KmerMatcherOptions::setSensitiveMode(const bool & value)
{
    sensitive_mode = value;
}

double KmerMatcherOptions::getThresholdPerc() const
{
    return threshold_perc;
}

void KmerMatcherOptions::setThresholdPerc(double & value)
{
    threshold_perc = value;
}

string KmerMatcherOptions::getKmerDepthOutFile() const
{
    return kmer_depth_out_file;
}

void KmerMatcherOptions::setKmerDepthOutFile(const string & value)
{
    kmer_depth_out_file = value;
}


string KmerMatcherOptions::getExtKmerDepthOutFile() const
{
    return ext_kmer_depth_out_file;
}

void KmerMatcherOptions::setExtKmerDepthOutFile(const string & value)
{
    ext_kmer_depth_out_file = value;
}

string KmerMatcherOptions::getOutputFile() const
{
    return output_file;
}

void KmerMatcherOptions::setOutputFile(const string & value)
{
    output_file = value;
}
// double KmerMatcherOptions::getDepthCutoff() const
// {
//     return depth_cutoff;
// }

// void KmerMatcherOptions::setDepthCutoff(uint64_t & value)
// {
//     depth_cutoff = value;
// }