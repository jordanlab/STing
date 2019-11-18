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


#include "KmerMatcher.h"

using namespace seqan;
using namespace std;

// -----------------------------------------------------------------------------------
// Constructors /destructors
// -----------------------------------------------------------------------------------

KmerMatcher::KmerMatcher()
{
    total_hits               = 0;
    total_kmers              = 0;
    total_reads              = 0;
    alleles_esa_index_exists = false;
}

KmerMatcher::KmerMatcher(const KmerMatcherOptions & options)
{
    total_hits               = 0;
    total_kmers              = 0;
    total_reads              = 0;
    alleles_esa_index_exists = false;
    this->options            = options;
}

KmerMatcher::~KmerMatcher()
{

}
// -----------------------------------------------------------------------------------
// Class methods (different from accessors)
// -----------------------------------------------------------------------------------


// // -----------------------------------------------------------------------------------
// // Getters and setters
// // -----------------------------------------------------------------------------------

// TUintVector KmerMatcher::getAlleleCounts() const
// {
//     return allele_counts;
// }

// void KmerMatcher::setAlleleCounts(TUintVector allele_counts)
// {
//     this->allele_counts = allele_counts;
// }

// const TCharStringSet& KmerMatcher::getAlleleSeqIds() const
// {
//     return allele_seq_ids;
// }

// void KmerMatcher::setAlleleSeqIds(const TCharStringSet& allele_seq_ids)
// {
//     this->allele_seq_ids = allele_seq_ids;
// }

// const TSequenceEsaIndex& KmerMatcher::getAllelesEsaIndex() const
// {
//     return alleles_esa_index;
// }

// void KmerMatcher::setAllelesEsaIndex(const TSequenceEsaIndex& alleles_esa_index)
// {
//     this->alleles_esa_index = alleles_esa_index;
// }

// const TLociTable& KmerMatcher::getLociTable() const
// {
//     return loci_table;
// }

// void KmerMatcher::setLociTable(const TLociTable& loci_table)
// {
//     this->loci_table = loci_table;
// }

// const KmerMatcherOptions& KmerMatcher::getOptions() const
// {
//     return options;
// }

// void KmerMatcher::setOptions(const KmerMatcherOptions& options)
// {
//     this->options = options;
// }

// const ArgumentParser& KmerMatcher::getParser() const
// {
//     return parser;
// }

// void KmerMatcher::setParser(const ArgumentParser& parser)
// {
//     this->parser = parser;
// }

// const TCharStringSet& KmerMatcher::getFastqFiles() const
// {
//     return fastq_files;
// }

// void KmerMatcher::setFastqFiles(const TCharStringSet& fastq_files)
// {
//     this->fastq_files = fastq_files;
// }

// const TStringSet& KmerMatcher::getSeqs() const
// {
//     return seqs;
// }

// void KmerMatcher::setSeqs(const TStringSet& seqs)
// {
//     this->seqs = seqs;
// }

// ----------------------------------------------------------------------------
void KmerMatcher::setUpInFastqFiles()
{
    // Construct a string set with all of the input read files (forward and reverse)
    TCharStringSet tmp_fastq_files;
    for (uint64_t i = 0; i < options.getInFastq1Files().size(); ++i)
    {
        appendValue(tmp_fastq_files, options.getInFastq1Files()[i]);
    }

    for (uint64_t i = 0; i < options.getInFastq2Files().size(); ++i)
    {
        appendValue(tmp_fastq_files, options.getInFastq2Files()[i]);
    }

    fastq_files = tmp_fastq_files;
}
// ----------------------------------------------------------------------------
void KmerMatcher::loadSeqIdsFile(string const & filename)
{
    if (checkFile(filename) == FILE_NOT_FOUND)
    {
        cerr << "ERROR: Index file '" << filename << "' not found" << endl;
        exit (1);
    }

    fstream infile(filename, ios::in);
    string line;
    while (getline(infile, line))
    {
        appendValue(allele_seq_ids, line);
    }

    infile.close();
}
// ----------------------------------------------------------------------------
void KmerMatcher::initAlleleCounts()
{
    // TUintVector tmp_allele_counts(countAllelesIndexSequences(), 0);
    // allele_counts.clear();
    // allele_counts = tmp_allele_counts;
    
    // uint64_t allele_counts_size = allele_counts.size();
    allele_counts.clear();
    allele_counts.resize(countAllelesIndexSequences(), 0);
}

// ----------------------------------------------------------------------------
// load and index into a TSequenceEsaIndex variable
// void KmerMatcher::loadAllelesEsaIndex(string const & index_prefix)
void KmerMatcher::loadAllelesEsaIndex()
{
    string index_prefix = options.getIndexFilename();

    // Check for index file
    string index_filename = index_prefix + ESA_INDEX_EXT;
    // cout << index_filename << endl;
    if (checkFile(index_filename) == FILE_NOT_FOUND)
    {
        cerr <<"ERROR: Index file '" << index_filename << "' not found." << endl;
        exit(1);
    }

    // cerr << "Loading the input index..." << endl;

    // Load the allele sequence ids file
    string alleles_seqs_ids_filename = index_prefix + ALLELES_SEQ_IDS_EXT;
    loadSeqIdsFile(alleles_seqs_ids_filename);

    // Load the index file
    if (! open(alleles_esa_index, index_prefix.c_str()))
    {
        cerr << "ERROR: Could not open the input index '" << index_prefix << "'\n";
        exit(1);
    }

    alleles_esa_index_exists = true;

    // cerr << "  Index successfully loaded" << endl << endl;
}
// ----------------------------------------------------------------------------
uint64_t KmerMatcher::countAllelesIndexSequences()
{
    if (alleles_esa_index_exists)
        return countSequences(alleles_esa_index);

    return 0;
}
// ----------------------------------------------------------------------------
void KmerMatcher::processReadsFiles(bool save_kmer_starts)
{
    if(length(fastq_files) > 0)
    {
        uint64_t n_fastq_files = length(fastq_files);
        TSequenceEsaIndexFinder finder_esa(alleles_esa_index);

        for (uint64_t i = 0; i < n_fastq_files; ++i)
        {
            TStringVector filename_tokens;
            split(filename_tokens, getFileName(toCString(fastq_files[i])), ".");
            string extension = filename_tokens[filename_tokens.size() - 1];

            if (extension == "gz")  // process a gzipped fastq file
            {
                igzstream in_fastq_gz_file(toCString(fastq_files[i]));
                if ( ! in_fastq_gz_file.good()) {
                    cerr << "ERROR: Opening file '" << toCString(fastq_files[i]) << "' failed.\n";
                    exit(1);
                }
                // processReadsByLine(in_fastq_gz_file, options.getKmerLength(), finder_esa);
                // if (save_kmer_starts)
                // {
                //     reprocessKmerizedReads(in_fastq_gz_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                // }
                // else
                // {
                    processReadsByLine(in_fastq_gz_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                // }
                in_fastq_gz_file.close();
            }
            else if (extension == "fastq" || extension == "fq")     // process a uncompressed fastq file
            {
                fstream in_fastq_file(toCString(fastq_files[i]), ios::in);
                if ( ! in_fastq_file.good()) {
                    cerr << "ERROR: Opening file '" << toCString(fastq_files[i]) << "' failed.\n";
                    exit(1);
                }
                // processReadsByLine(in_fastq_file, options.getKmerLength(), finder_esa);
                // if (save_kmer_starts)
                // {
                //     reprocessKmerizedReads(in_fastq_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                // }
                // else
                // {
                    processReadsByLine(in_fastq_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                // }
                in_fastq_file.close();
            }
            else {
                cerr << "ERROR: The input file '" << toCString(fastq_files[i]) << "' has an invalid extension." << endl;
                exit(1);
            }
        }
    }
}
// ----------------------------------------------------------------------------
template <typename T>
void KmerMatcher::reprocessKmerizedReads(T &                       in_fastq_file,
                                         uint64_t                  fastq_file_idx,
                                         uint64_t                  k,
                                         TSequenceEsaIndexFinder & finder_esa,
                                         bool                      save_kmer_starts)
{
    uint64_t line_idx = 0;
    uint64_t n_lines  = 0;
    uint64_t start    = 0;
    // for(uint64_t line_idx : kmerized_reads_indices.at(fastq_file_idx))
    uint64_t n_kmerized_reads = kmerized_reads_indices.at(fastq_file_idx).size();
    // cerr << "DEBUG: KmerMatcher::reprocessKmerizedReadss() - file " << fastq_file_idx << " - n_kmerized_reads: " << n_kmerized_reads << endl;
    for(uint64_t i = 0; i < n_kmerized_reads; ++i)
    {
        string line;
        line_idx = kmerized_reads_indices.at(fastq_file_idx).at(i);
        
        if (i > 0) 
            start = kmerized_reads_indices.at(fastq_file_idx).at(i-1);
        
        n_lines = line_idx - start;
        
        advanceNLines(in_fastq_file, n_lines);
        getline(in_fastq_file, line);
        TCharStringSet ids;
        StringSet<Dna5String> seqs;

        TSequence seq = line;
        // process the current read
        
        if (k <= length(seq))
        {
            processRead(seq, k, finder_esa, save_kmer_starts);
            // cerr << "DEBUG: KmerMatcher::reprocessKmerizedReqads() - 300 - read: " << i << " " << seq << endl;
            // find reverse-complement of the current read
            reverseComplement(seq);     // SeqAn function
            processRead(seq, k, finder_esa, save_kmer_starts);
            // cerr << "DEBUG: KmerMatcher::reprocessKmerizedReqads() - 304 - read: " << i << " " << seq << endl;
        }
        else
        {
            cerr << "ERROR: The k-mer length is greater than the length of reads."
                    " Please check the length of the reads and choose a shorter k-mer size." << endl;
            exit(1);
        }
    } // for kmerized_reads_indices
}
// ----------------------------------------------------------------------------
template <typename T>
void KmerMatcher::processReadsByLine(T &                       in_fastq_file,
                                     uint64_t                  fastq_file_idx,
                                     uint64_t                  k,
                                     TSequenceEsaIndexFinder & finder_esa,
                                     bool save_kmer_starts)
{
    uint64_t counter = 0;
    string line;
    while (getline(in_fastq_file, line))
    {
        if (++counter % 4 == 2) { // get only the sequence line

            TCharStringSet ids;
            StringSet<Dna5String> seqs;

            TSequence seq = line;
            // process the current read
            
            if (k <= length(seq))
            {
                KmerMatcher::ProcessReadResult readRes1;
                KmerMatcher::ProcessReadResult readRes2;
                
                readRes1 = processRead(seq, k, finder_esa, save_kmer_starts);
                // find reverse-complement of the current read
                reverseComplement(seq);     // SeqAn function
                readRes2 = processRead(seq, k, finder_esa, save_kmer_starts);
                
                // if (!save_kmer_starts &&
                //     (readRes1 == KmerMatcher::PROCESSED_READ ||
                //      readRes2 == KmerMatcher::PROCESSED_READ))
                if (readRes1 == KmerMatcher::PROCESSED_READ ||
                    readRes2 == KmerMatcher::PROCESSED_READ)
                {
                    ++total_reads;
                    
                    // Add the line of the k-merized read
                    // if (kmerized_reads_indices.find(fastq_file_idx) == kmerized_reads_indices.end()) 
                    if (kmerized_reads_indices.find(fastq_file_idx) == kmerized_reads_indices.end()) 
                        kmerized_reads_indices[fastq_file_idx] = {};

                    kmerized_reads_indices.at(fastq_file_idx).push_back(counter);
                }
            }
            else
            {
                cerr << "ERROR: The k-mer length is greater than the length of reads."
                        " Please check the length of the reads and choose a shorter k-mer size." << endl;
                exit(1);
            }
        }
    }
}
// ----------------------------------------------------------------------------
template <typename T>
void KmerMatcher::advanceNLines(T & file, uint64_t n_lines)
{   
    uint64_t stop = n_lines - 1;
    for(uint64_t i = 0; i < stop; ++i)
    {
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
}

// ----------------------------------------------------------------------------
KmerMatcher::ProcessReadResult 
KmerMatcher::processRead(TSequence const &           read,
                         uint64_t                    k,
                         Finder<TSequenceEsaIndex> & finder_esa,
                         bool                        save_kmer_starts)
{
    ProcessReadResult result = KmerMatcher::UNPROCESSED_READ;
    // look for middle kmer
    uint64_t read_length = length(read);
    uint64_t start_pos     = (read_length - k) / 2;
    TSequence middle_kamer = infix(read, start_pos, (start_pos + k + 1));
    clear(finder_esa);
    
    // cout << "Read: " << read << endl;
    // cout << "   Middle kmer: " << middle_kamer << endl;
    // cout << "   Found: " << find(finder_esa, middle_kamer) << endl;
    // clear(finder_esa);

    if (find(finder_esa, middle_kamer))   // middle k-mer found
    {
        // ++total_reads;
        result = KmerMatcher::PROCESSED_READ;
        clear(finder_esa);
        uint64_t stop_base = read_length - k;
        for (uint64_t j = 0; j <= stop_base; ++j) // k-merize all the read
        {
            TSequence kmer       = infix(read, j, j+k);
            uint64_t pos_counter = 0;
            while (find(finder_esa, kmer)) // search each k-mer in the index
            {
                TSequenceEsaIndexPosition pos = position(finder_esa);
                ++allele_counts[pos.i1];
                ++total_hits;
                ++pos_counter;
                
                if (save_kmer_starts) 
                {
                    // Add the start of the hit in the corresponding allele
                    if (kmer_starts.find(pos.i1) == kmer_starts.end()) 
                        kmer_starts[pos.i1] = {};

                    kmer_starts[pos.i1].push_back(pos.i2);
                }
            }
            
            if (pos_counter > 0) ++total_kmers;
            
            clear(finder_esa);
        }
    }
    return result;
}
// // ----------------------------------------------------------------------------
// template <typename T>
// void KmerMatcher::reprocessKmerizedReads(T &                       in_fastq_file,
//                                          uint64_t                  fastq_file_idx,
//                                          uint64_t                  k,
//                                          TSequenceEsaIndexFinder & finder_esa,
//                                          bool                      save_kmer_starts)
// {
//     uint64_t line_idx = 0;
//     uint64_t n_lines  = 0;
//     uint64_t start    = 0;
//     // for(uint64_t line_idx : kmerized_reads_indices.at(fastq_file_idx))
//     uint64_t n_kmerized_reads = kmerized_reads_indices.at(fastq_file_idx).size();
//     // cerr << "DEBUG: reprocessKmerizedReads()" << endl;
//     for(uint64_t i = 0; i < n_kmerized_reads; ++i)
//     {
//         string line;
//         line_idx = kmerized_reads_indices.at(fastq_file_idx).at(i);
        
//         if (i > 0) 
//             start = kmerized_reads_indices.at(fastq_file_idx).at(i-1);
        
//         n_lines = line_idx - start;
        
//         advanceNLines(in_fastq_file, n_lines);
//         getline(in_fastq_file, line);
//         TCharStringSet ids;
//         StringSet<Dna5String> seqs;

//         TSequence seq = line;
//         // process the current read
        
//         if (k <= length(seq))
//         {
//             processRead(seq, k, finder_esa, save_kmer_starts);
//             // find reverse-complement of the current read
//             reverseComplement(seq);     // SeqAn function
//             processRead(seq, k, finder_esa, save_kmer_starts);
//         }
//         else
//         {
//             cerr << "ERROR: The k-mer length is greater than the length of reads."
//                     " Please check the length of the reads and choose a shorter k-mer size." << endl;
//             exit(1);
//         }
//     } // for kmerized_reads_indices
// }
// ----------------------------------------------------------------------------
void KmerMatcher::normalizeCounts()
{
    uint64_t allele_counts_size = allele_counts.size();
    for (uint64_t i = 0; i < allele_counts_size; ++i) {
        uint64_t allele_lenght = length(indexText(alleles_esa_index)[i]);
        double norm_count = (double)allele_counts.at(i)/(double)(allele_lenght);
        norm_allele_counts.push_back(norm_count);
        
        // string seq_id = toCString(allele_seq_ids[i]);
    }
}
// ----------------------------------------------------------------------------
void KmerMatcher::loadLociTable(string const & filename)
{
    // string filename = options.getIndexFilename() + string(LOCI_TABLE_EXT)

    fstream loci_file(filename, ios::in);

    if (! loci_file.good())
    {
        cerr << "ERROR: Opening file '" << filename << "' failed.\n";
        exit(1);
    }

    string line;
    while (getline(loci_file, line))
    {
        vector<string> fields;
        if (!line.empty() && line.at(0) != COMMENT_CHAR) // The current line is not an empty/comment line
        {
          // cout << line << endl;
          split(fields, line, "\t");
          // SLocusRecord loci = {fields[0], fields[1], fields[2], (uint64_t)stoi(fields[3]), (uint64_t)stoi(fields[4])};
          SLocusRecord loci = {fields[0], fields[1], (uint64_t)stoi(fields[2]), (uint64_t)stoi(fields[3])};
          loci_table.push_back(loci);
        }
    }
}

// ----------------------------------------------------------------------------
TStringVector KmerMatcher::getAlleleCoverageForPrinting(uint64_t const & idx)
{
    TStringVector out_vector;
    for (auto a_cov : coverages.at(idx)) 
    {
        stringstream cov;
        
        cov  << a_cov.cov_bases 
             << "/" 
             << a_cov.seq_size 
             << "(" 
             << getNumberStringWithPrecision(a_cov.coverage, 1) 
             << "%)";
             
        out_vector.push_back(cov.str());
    }
    return out_vector;
}
// // ----------------------------------------------------------------------------
// void KmerMatcher::printAlleleCoverage(uint64_t const & idx, uint64_t n_separators = 4)
// {
//     cout << endl
//          << "allele_coverage";
    
//     // room for the initial fields in the previous line (profile)
//     for(unsigned i = 0; i < n_separators; ++i) {
//         cout << FIELD_SEPARATOR;
//     }

//     for (auto a_cov : coverages.at(idx)) 
//     {
//         cout << a_cov.cov_bases 
//              << "/" 
//              << a_cov.seq_size 
//              << "(" 
//              << getNumberStringWithPrecision(a_cov.coverage, 1) 
//              << "%)"
//              << FIELD_SEPARATOR;
//     }
// }
// ----------------------------------------------------------------------------
TStringVector KmerMatcher::getMeanKmerDepthForPrinting(uint64_t const & idx)
{
    TStringVector out_vector;
    
    for (auto a_depth : mean_kmer_depths.at(idx)) 
        out_vector.push_back(getNumberStringWithPrecision(a_depth, 1));
    
    return out_vector;
}
// ----------------------------------------------------------------------------
// void KmerMatcher::printMeanKmerDepth(uint64_t const & idx, uint64_t n_separators = 4)
// {
//     cout << endl
//          << "mean_kmer_depth";

//     // room for the initial fields in the previous line (profile)
//     for(unsigned i = 0; i < n_separators; ++i) {
//         cout << FIELD_SEPARATOR;
//     }
    
//     for (auto a_depth : mean_kmer_depths.at(idx)) 
//     {
//         cout << getNumberStringWithPrecision(a_depth, 1)
//              << FIELD_SEPARATOR;
//     }
// }
// ----------------------------------------------------------------------------

/**
 * @brief      Saves the k-mer depth data of a sequence (allele/gene), to a file
 *
 * @param[in]  output_file_name    The output file name
 * @param      allele_kmer_depths  The k-mer depths of the sequence
 * @param      seq_idx             The sequence index of the sequence
 * @param[in]  append              Flag to indicate if data must be appended or not.
 *
 * @return     A SaveFileResult that indicates success (SAVE_SUCCESS) or fail (SAVE_ERROR).
 */
SaveFileResult KmerMatcher::saveKmerDepthData(const string &   output_file_name, 
                                              TUintVector &    seq_kmer_depths,
                                              uint64_t const & seq_idx,
                                              bool             append = true)
{
    fstream out_file;
    
    if (!append) 
    {
        out_file.open(output_file_name, ios::trunc |  ios::out);
    } else
    {
        out_file.open(output_file_name, ios::out | ios::app);
    }
    
    if (!out_file)
    {
        cerr << "ERROR: Per-base k-mer depth file could not be saved." << endl;
        return SAVE_ERROR;
    }
    
    if (!append) 
    {
        // FIX: remove n_best_prof column when all traces of 1st and 2nd 
        // best alleles are removed
        // Write header
        // out_file << "n_best_prof" << FIELD_SEPARATOR
        out_file << "allele"              << FIELD_SEPARATOR
                 << "position"            << FIELD_SEPARATOR
                 << "kmer_depth"          << FIELD_SEPARATOR
                 << "weighted_kmer_depth" 
		         << endl;
    }
    
    string seq_id = toCString(allele_seq_ids[seq_idx]);
    uint64_t seq_kmer_depths_size = seq_kmer_depths.size();
    
    for(uint64_t i = 0; i < seq_kmer_depths_size; ++i)
    {
        // FIX: remove hardcoded 1 when all traces of 1st and 2nd 
        // best alleles are removed
        // out_file << 1 << FIELD_SEPARATOR
        
        float weighted_kmer_depth = 0.0;
        uint64_t depth = seq_kmer_depths.at(i);
        if (depth > 0)
        {
            weighted_kmer_depth = float (total_kmers / depth);
        }
        
        out_file << seq_id                << FIELD_SEPARATOR
                 << i+1                   << FIELD_SEPARATOR
                 // << seq_kmer_depths.at(i)
                 << depth << FIELD_SEPARATOR
                 << weighted_kmer_depth
                 << endl;
    } // End per-base depths loop
    
    out_file.close();
    return SAVE_SUCCESS;
}
// ----------------------------------------------------------------------------
string KmerMatcher::fastqFilesToString()
{
    // Fastq files to string
    string fq_files = "";
    string sep = ",";
    for (uint i = 0; i < length(fastq_files); ++i)
    {
        if (i == length(fastq_files) - 1) sep = "";
        fq_files += getFileName(toCString(fastq_files[i])) + sep;
    }
    return fq_files;
}
// ----------------------------------------------------------------------------
void KmerMatcher::printMessage(string const & message)
{
    if (options.getVerbose())
        cerr << message;
}
// ----------------------------------------------------------------------------
void KmerMatcher::printTime(string const & message)
{
    if (options.getPrintTime())
        cerr << message;
}
// ----------------------------------------------------------------------------
// template void KmerMatcher::advanceNLines(igzstream&, uint64_t);
// template void KmerMatcher::advanceNLines(fstream&, uint64_t);
// template void KmerMatcher::processReadsByLine(igzstream&,
//                                  uint64_t,
//                                  uint64_t,
//                                  TSequenceEsaIndexFinder &,
//                                  bool);
// template void KmerMatcher::processReadsByLine(fstream&,
//                                  uint64_t,
//                                  uint64_t,
//                                  TSequenceEsaIndexFinder &,
//                                  bool);
template void KmerMatcher::reprocessKmerizedReads(igzstream&,
                                     uint64_t,
                                     uint64_t,
                                     TSequenceEsaIndexFinder &,
                                     bool);
template void KmerMatcher::reprocessKmerizedReads(fstream&,
                                     uint64_t,
                                     uint64_t,
                                     TSequenceEsaIndexFinder &,
                                     bool);
