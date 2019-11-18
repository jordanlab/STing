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


#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

// #include <seqan/index.h>
// #include <seqan/seq_io.h>
// #include "types.hpp"
// #include <seqan/basic.h>

// ============================================================================
// Structs
// ============================================================================
struct SIndexerOptions
{
    CharString config_file_filename;
    CharString out_prefix_filename;
    
    SIndexerOptions() :
        config_file_filename(), 
        out_prefix_filename()
    {}
};
// ----------------------------------------------------------------------------
struct STyperOptions
{
    CharString    index_filename;
    
    TStringVector in_fastq_1_files;
    TStringVector in_fastq_2_files;
    
    uint64_t      kmer_length; 
    uint64_t      verbosity_level;
    
    STyperOptions() :
        index_filename(),
        in_fastq_1_files(),
        in_fastq_2_files(),
        kmer_length(30),
        verbosity_level(0)
    {}
};
// ----------------------------------------------------------------------------
struct SLocusRecord
{
    std::string id;
    std::string filename;
    uint64_t n_seqs;
    uint64_t last_seq_idx;
};
// ----------------------------------------------------------------------------
struct ProfilesTable
{
    std::vector<std::string> cols;
    std::map<std::string, std::vector<std::string> > table;
};
// ----------------------------------------------------------------------------
struct SAlleleCount
{
    int idx;
    double count;

    SAlleleCount():
        idx(-1),
        count(-1.0)
    {}
    
    bool operator<(const SAlleleCount & other) const
    {
        return count < other.count;
    }
    
};
// ----------------------------------------------------------------------------
struct SeqCoverage
{
    uint64_t cov_bases;
    uint64_t seq_size;
    double   coverage;
    
    SeqCoverage() :
        cov_bases(0),
        seq_size(0),
        coverage(0.0)
    {}
};
// ----------------------------------------------------------------------------

#endif /*DATA_STRUCTURES_HPP*/