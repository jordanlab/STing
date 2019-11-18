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


#ifndef KMER_MATCHER_H
#define KMER_MATCHER_H

// ============================================================================
// Includes
// ============================================================================
// std
// ----------------------------------------------------------------------------
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

// ----------------------------------------------------------------------------
// SeqAn (https://www.seqan.de/)
// ----------------------------------------------------------------------------
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/index.h>
// #include <seqan/misc/interval_tree.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <seqan/stream.h>
// ----------------------------------------------------------------------------
// Gzstream (http://www.cs.unc.edu/Research/compgeom/gzstream/ - https://gist.github.com/piti118/1508048)
// ----------------------------------------------------------------------------
#include <gzstream/gzstream.h>
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Application
// ----------------------------------------------------------------------------
#include "KmerMatcherOptions.h"
#include "types.hpp"
#include "data_structures.hpp"
#include "custom_utils.hpp"
#include "SequenceInfo.h"
// ============================================================================
// Namespaces
// ============================================================================
using namespace seqan;
using namespace std;

class KmerMatcher
{
    protected:
        TUintVector                    allele_counts;
        TDoubleVector                  norm_allele_counts;
        TCharStringSet                 allele_seq_ids;
        TSequenceEsaIndex              alleles_esa_index;
        TCharStringSet                 fastq_files;
        TLociTable                     loci_table;
        uint64_t                       total_hits;
        uint64_t                       total_kmers;
        uint64_t                       total_reads;
        
        // vector<int>                    first_max_allele;            // indices of alleles with the top kmer hits
        
        TUMapIntervalVector            seq_intervals;
        TUMapUintUintVector            kmer_starts;
        vector<TSeqCovVector>          coverages;
        vector<vector<TUintVector > >  perbase_kmer_depths;
        vector<TDoubleVector>          mean_kmer_depths;
        
        KmerMatcherOptions             options;
        ArgumentParser                 parser;
        TStringSet                     seqs;
        TUMapUintUintVector            kmerized_reads_indices;
        
        bool                           alleles_esa_index_exists;
        
    public:
        
        enum ProcessReadResult
        {
             PROCESSED_READ,       // Read was processed
             UNPROCESSED_READ      // Read was not processed
        };
        
        // constructors / destructors
        KmerMatcher();
        KmerMatcher(const KmerMatcherOptions & options);
        virtual ~KmerMatcher();
        
    protected:
        // Class methods (different form constructors, accessors)
        void setUpInFastqFiles();
        void loadSeqIdsFile(string const & filename);
        // void loadAllelesEsaIndex(string const & index_prefix);
        void loadAllelesEsaIndex();
        void initAlleleCounts();
        uint64_t countAllelesIndexSequences();
        void processReadsFiles(bool save_kmer_starts);
        virtual ProcessReadResult 
        processRead(TSequence const &         read, 
                    uint64_t                  k, 
                    TSequenceEsaIndexFinder & finder_esa,
                    bool save_kmer_starts);
        void normalizeCounts();
        void loadLociTable(string const & filename);
                
        TStringVector getAlleleCoverageForPrinting(uint64_t const & idx);
        // void printAlleleCoverage(uint64_t const & idx, uint64_t n_separators);
        TStringVector getMeanKmerDepthForPrinting(uint64_t const & idx);
        // void printMeanKmerDepth(uint64_t const & idx, uint64_t n_separators);
        SaveFileResult saveKmerDepthData(const string &   output_file_name, 
                                         TUintVector &    gene_kmer_depths,
                                         uint64_t const & seq_idx,
                                         bool             append);
        string fastqFilesToString();
        void printMessage(string const & message);
        void printTime(string const & message);
        template <typename T>
        void advanceNLines(T & file, uint64_t n_lines);
        template <typename T>
        void processReadsByLine(T &                       in_fastq_file, 
                                uint64_t                  fastq_file_idx, 
                                uint64_t                  k, 
                                TSequenceEsaIndexFinder & finder_esa,
                                bool save_kmer_starts);
        template <typename T>
        void reprocessKmerizedReads(T &                       in_fastq_file,
                                    uint64_t                  fastq_file_idx,
                                    uint64_t                  k,
                                    TSequenceEsaIndexFinder & finder_esa,
                                    bool save_kmer_starts);
};

#endif // KMER_MATCHER_H
