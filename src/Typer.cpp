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


#include "Typer.h"

using namespace seqan;
using namespace std;

// -----------------------------------------------------------------------------------
// Constructors /destructors
// -----------------------------------------------------------------------------------

Typer::Typer()
{
    alleles_esa_index_exists   = false;
}

Typer::Typer(const KmerMatcherOptions & options)
{
    alleles_esa_index_exists = false;
    this->options            = options;
}

Typer::~Typer()
{

}
// -----------------------------------------------------------------------------------
// Class methods (different from accessors)
// -----------------------------------------------------------------------------------
void Typer::setupArgumentParser(const string & program_name,
                                const string & program_description)
{
    // string program_name             = TYPER_APP_NAME;

    string in_idx_param_name             = "INDEX_PREFIX_FILENAME";
    string in_fastq_1_param_name         = "FASTQ1";
    string in_fastq_2_param_name         = "FASTQ2";
    string sample_name_name              = "SAMPLE_NAME";
    string kmer_length_param_name        = "KMER_LENGTH";
    string n_top_alleles_param_name      = "N_TOP_ALLELES";
    string kmer_depth_out_param_name     = "KMER_DEPTH_FILENAME";
    string output_file_param_name        = "OUTPUT_FILENAME";
    string ext_kmer_depth_out_param_name = "EXT_KMER_DEPTH_FILENAME";
    string kmer_depth_cutoff_param_name  = "KMER_DEPTH_CUTOFF";
    string verbose_param_name            = "VERBOSITY_LEVEL";

    // Setup ArgumentParser.
    setAppName(parser, program_name);
    setShortDescription(parser, "");
    setCategory(parser, "Bacterial Typing");
    setVersion(parser, TYPER_APP_VERSION);
    setDate(parser, TYPER_APP_UPDATE);

    // Usage message
    addUsageLine(parser, "-x <\\fI" +  in_idx_param_name + "\\fP> "
                 "-1 <\\fI" +  in_fastq_1_param_name + "\\fP> [\\fIOPTIONS\\fP]");

    // Description
    addDescription(parser, program_description);

    // Options
    addSection(parser, "Required input parameters");
    addOption(parser, ArgParseOption(
                  "x", "index-prefix", "Index prefix filename.",
                  ArgParseArgument::INPUT_FILE, in_idx_param_name));
    setRequired(parser, "index-prefix");

    addOption(parser, ArgParseOption(
                  "1", "fastq-1-files", "Files with #1 mates, paired with "
                  "files in \\fI" + in_fastq_2_param_name + "\\fP. Valid file "
                  "extensions are \\fI.fq\\fP, \\fI.fastq\\fP (uncompressed fastq), "
                  "and \\fI.gz\\fP (gzip'ed fastq).",
                  ArgParseArgument::INPUT_FILE, in_fastq_1_param_name, true));
    setRequired(parser, "fastq-1-files");
    
    // Input options
    addSection(parser, "Input options");
    addOption(parser, ArgParseOption(
                  "2", "fastq-2-files", "Files with #2 mates, paired with "
                  "files in \\fI" + in_fastq_1_param_name + "\\fP. Valid file "
                  "extensions are \\fI.fq\\fP, \\fI.fastq\\fP (uncompressed fastq), "
                  "and \\fI.gz\\fP (gzip'ed fastq).",
                  ArgParseArgument::INPUT_FILE, in_fastq_2_param_name, true));

    addOption(parser, ArgParseOption(
                  "s", "sample-name", "Name of the sample to be analized.",
                  ArgParseArgument::STRING, sample_name_name));
    // setDefaultValue(parser, "sample-name", 30);

    addOption(parser, ArgParseOption(
                  "k", "kmer-length", "Length of the k-mers to process the input reads.",
                  ArgParseArgument::INTEGER, kmer_length_param_name));
    setDefaultValue(parser, "kmer-length", DEFAULT_KMER_LENGHT);
    
    addOption(parser, ArgParseOption(
        "n", "n-top-alleles", "Number of top alleles by k-mer frequencies, from which best alleles will be called.",
        ArgParseArgument::INTEGER, n_top_alleles_param_name));
    setDefaultValue(parser, "n-top-alleles", DEFAULT_TOP_N_ALLELES);

    // Output options
    addSection(parser, "Output options");
    addOption(parser, ArgParseOption(
                  "c", "kmer-counts", "Select and print the normalized k-mer counts at "
                  "each locus (k-mer hits divided by allele length)."));

    addOption(parser, ArgParseOption(
                  "a", "allele-cov", "Select to calculate and print the percent of the length "
                  "of each allele that is covered by k-mer hits. Note that this calculation will "
                  "require more execution time."));

    // addOption(parser, ArgParseOption(
    //     "e", "depth-cutoff", "Minimum number of k-mers at the first and last base of an allele to cosider them as covered.",
    //     ArgParseArgument::INTEGER, kmer_depth_cutoff_param_name));
    // setDefaultValue(parser, "depth-cutoff", 3);

    addOption(parser, ArgParseOption(
                  "d", "kmer-depth", "Select to calculate and print the average k-mer depth of each allele. "
                  "Note that this calculation will require more execution time."));
    
    addOption(parser, ArgParseOption(
        "t", "k-depth-file", "Output filename to save the detailed per-base k-mer depth data.",
        ArgParseArgument::OUTPUT_FILE, kmer_depth_out_param_name));
    
   // addOption(parser, ArgParseOption(
   //      "e", "ext-k-depth-file", "Output filename to save the extended per-base k-mer depth data.",
   //      ArgParseArgument::OUTPUT_FILE, ext_kmer_depth_out_param_name));
    
    addOption(parser, ArgParseOption(
        "o", "output-file", "Output filename to save the typing results (instead of stdout).",
        ArgParseArgument::OUTPUT_FILE, output_file_param_name));

    addOption(parser, ArgParseOption(
                  "y", "print-tidy", "Select to print results in a tidy format."));
    
    addOption(parser, ArgParseOption(
                  "v", "verbose", "Select to print informative messages (to stderr)."));
    
    addOption(parser, ArgParseOption(
                  "m", "time", "Select to print time in seconds (to stderr) of each step of the typing process."));
    hideOption(parser, "time");
    
    // Presets
    addSection(parser, "Presets");
    addOption(parser, ArgParseOption(
                  "", "fast", "(Default). Select to set the \\fIfast\\fP mode to call alleles "
                  "based solely on k-mer frequencies. The best allele of each locus is that with "
                  "the highest k-mer hit frequency."));
    addOption(parser, ArgParseOption(
                  "", "sensitive", "Select to set the \\fIsensitive\\fP mode to call alleles "
                  "based on k-mer frequencies and coverage information. The best allele of each " 
                  "locus is that with the highest k-mer hit frequency and the highest allele coverage. "
                  "Allele ties are solved selecting the allele with the minimum k-mer depth standard "
                  "deviation."));
    
    // addOption(parser, ArgParseOption(
    //     "f", "fast-mode", "Select to set the fast mode which excludes allele coverage and k-mer depth (mean and per-base) calculation. Note that this option is mutually exclusive with -a, -d and -t options."));

    // Additional information
    addTextSection(parser, "Note on " + in_fastq_1_param_name + " and "
                   "" + in_fastq_2_param_name);
    addText(parser, "\n\\fI" + in_fastq_1_param_name + "\\fP and "
            "\\fI" + in_fastq_2_param_name + "\\fP can be comma-separated lists "
            "(no whitespace) and can be specified many times, "
            "e.g., -1 file1.fq,file2.fq -1 file3.fq");

    addTextSection(parser, "Output conventions");
    
    addTextSubSection(parser, "Status values:");
    addListItem(parser, "\\fBst_predicted\\fP", "ST inferred from k-mer "
                "hits and profiles table.");

    addListItem(parser, "\\fBno_kmer_hits\\fP", "There are no k-mer "
                "hits for one or more of the loci to\n infer a profile and its "
                "associated ST. Probable causes:");
    addListItem(parser, "", "1) k-mer length not adequate (too long),");
    addListItem(parser, "", "2) low quality data/too many N's in the data.");

    addListItem(parser, "\\fBno_st_in_table\\fP", "There is no ST associated "
                "to the inferred allelic profile. Probable causes:");
    addListItem(parser, "", "1) k-mer length not adequate (too short),");
    addListItem(parser, "", "2) this could be a new allelic profile.");

    addTextSubSection(parser, "Other values/symbols:");

    addListItem(parser, "\\fB" + string(NA_STRING) + "\\fP", "No ST associated "
                "or no k-mer hits for any allele of a locus.");

    addListItem(parser, "\\fB" + string(LOW_CONFIDENCE_SYMBOL) + "\\fP", "Low "
                "confidence. Allele predicted with a length coverage "
                "below 100%. Probable causes:");
    addListItem(parser, "", "1) No enough k-mer hits to cover the whole length "
                "of the allele");
    addListItem(parser, "", "2) No enough k-mer depth in some part along the "
                "range [10, length-10] bp of the allele to consider it as covered.");
}
// ----------------------------------------------------------------------------
ArgumentParser::ParseResult
Typer::parseCommandLine(int argc, char const ** argv)
{
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    string index_filename;
    getOptionValue(index_filename, parser, "index-prefix");
    options.setIndexFilename(index_filename);

    TStringVector fastq1 = getOptionValues(parser, "fastq-1-files");

    TStringVector fastq2;
    if (isSet(parser, "fastq-2-files"))
        fastq2 = getOptionValues(parser, "fastq-2-files");

    options.setInFastq1Files(mergeOptionValues(fastq1));
    options.setInFastq2Files(mergeOptionValues(fastq2));

    string sample_name = "";
    if (isSet(parser, "sample-name"))
    {
        getOptionValue(sample_name, parser, "sample-name");
        options.setSampleName(sample_name);
    }

    uint64_t kmer_length;
    getOptionValue(kmer_length, parser, "kmer-length");
    options.setKmerLength(kmer_length);
    
    uint64_t n_top_alleles;
    getOptionValue(n_top_alleles, parser, "n-top-alleles");
    options.setNTopAlleles(n_top_alleles);

    options.setPrintKmerCounts(isSet(parser, "kmer-counts"));

    options.setPrintAlleleCov(isSet(parser, "allele-cov"));

    // uint64_t depth_cutoff;
    // getOptionValue(depth_cutoff, parser, "depth-cutoff");
    // options.setDepthCutoff(depth_cutoff);

    options.setPrintKmerDepth(isSet(parser, "kmer-depth"));

    string kmer_depth_out_filename = "";
    if (isSet(parser, "k-depth-file"))
        getOptionValue(kmer_depth_out_filename, parser, "k-depth-file");

    options.setKmerDepthOutFile(kmer_depth_out_filename);
    
    // string ext_kmer_depth_out_filename = "";
    // if (isSet(parser, "ext-k-depth-file")) 
    //     getOptionValue(ext_kmer_depth_out_filename, parser, "ext-k-depth-file");
    
    // options.setExtKmerDepthOutFile(ext_kmer_depth_out_filename);
    
    string output_filename = "";
    if (isSet(parser, "output-file")) 
        getOptionValue(output_filename, parser, "output-file");
    
    options.setOutputFile(output_filename);
    
    // if (isSet(parser, "fast-mode") && (isSet(parser, "allele-cov") ||
    //                                    isSet(parser, "kmer-depth") ||
    //                                    isSet(parser, "k-depth-file")))
    // {
    //     cerr << "ERROR: The option -f (--fast-mode) can not be used along with -a (--allele-cov), -d (--kmer-depth) -t (--k-depth-file)." << endl;
    //     return ArgumentParser::PARSE_ERROR;
    // }

    options.setFastMode(isSet(parser, "fast"));
    
    options.setSensitiveMode(isSet(parser, "sensitive"));

    options.setTidyResults(isSet(parser, "print-tidy"));

    options.setVerbose(isSet(parser, "verbose"));
    
    options.setPrintTime(isSet(parser, "time"));

    return ArgumentParser::PARSE_OK;
}
// ----------------------------------------------------------------------------
void Typer::setupReducedIndexFromTopAlleles()
{
    // cerr << "DEBUG: setupReducedIndexFromTopAlleles() - loci_table.size() " << loci_table.size() << endl;
    TStringSet best_allele_sequences;
    string allele_seq_id;
    // uint64_t counter = 0;
    // vector<int> allele_indices;
    for (auto&& locus : loci_table)
    {
        TAlleleCountVector alleles = top_alleles.at(locus.id);
        // cerr << "DEBUG: setupReducedIndexFromTopAlleles() - " << ++counter << " - locus: " << locus.id << " - alleles: " << alleles.size() << endl;
        if (alleles.size() > 0)
        {
            for (auto&& allele : alleles) 
            {
                TSequence seq = indexText(alleles_esa_index)[allele.idx];
                // cerr << "DEBUG: setupReducedIndexFromTopAlleles() - " << counter << " - allele: " << allele_seq_ids[allele.idx] << " - size: " << length(seq) << endl;
                // cerr << "DEBUG: setupReducedIndexFromTopAlleles() - " << counter << " - allele: " << allele_seq_ids[allele.idx] << endl;
                appendValue(best_allele_sequences, seq);
                // cerr << "DEBUG: setupReducedIndexFromTopAlleles() - Allele seq added to the StringSet"<< endl;
                reduced_index_allele_indices.push_back(allele.idx);
                // allele_indices.push_back(allele.idx);
                // cerr << "DEBUG: setupReducedIndexFromTopAlleles() - Allele idx added to vector" << endl;
            }
        }
    }
    // reduced_index_allele_indices = allele_indices;
    // if (length(best_allele_sequences) > 0)
    // {
    // cout << "DEBUG 1.2" << endl;
    indexText(top_alleles_esa_index) = best_allele_sequences;
    // indexRequire(esa_index, FibreSA());
    indexRequire(top_alleles_esa_index, EsaSA());
    // }
    
    // cout << "DEBUG 1.3" << endl;
    // int index_size = length(indexText(top_alleles_esa_index));
    
    // for(int i = 0; i < index_size; ++i) {
    //     allele_seq = indexText(top_alleles_esa_index)[i];
    //     allele_seq_id = toCString(allele_seq_ids[first_max_allele.at(i)]);
    //     cout << "allele: " << allele_seq_id << "\tallele_seq: " << allele_seq << endl;
    // }
    
}
// ----------------------------------------------------------------------------
void Typer::processReadsFiles(bool save_kmer_starts)
{
    // cerr << "DEBUG: Typer implementation of processReadsFiles()" << endl;
    uint64_t n_fastq_files = length(fastq_files);
    if(n_fastq_files > 0)
    {
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
                if (save_kmer_starts)
                {
                    // Are there useful reads in the surrent i-th reads file?
                    if (kmerized_reads_indices.find(i) != kmerized_reads_indices.end())
                    {
                        TSequenceEsaIndexFinder finder(top_alleles_esa_index);
                        // reprocessKmerizedReads(in_fastq_gz_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                        reprocessKmerizedReads(in_fastq_gz_file, i, options.getKmerLength(), finder, save_kmer_starts);
                    }
                }
                else
                {
                    processReadsByLine(in_fastq_gz_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                }
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
                if (save_kmer_starts)
                {
                    // Are there useful reads in the surrent i-th reads file?
                    if (kmerized_reads_indices.find(i) != kmerized_reads_indices.end())
                    {
                        // cerr << "qDEBUG: fastq file: " << i << endl;
                        TSequenceEsaIndexFinder finder(top_alleles_esa_index);
                        // cerr << "DEBUG: processReadsFiles - Before reprocessKmerizedReads" << endl;
                        // reprocessKmerizedReads(in_fastq_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                        // cerr << "DEBUG: processReadsFiles - After reprocessKmerizedReads" << endl;
                        reprocessKmerizedReads(in_fastq_file, i, options.getKmerLength(), finder, save_kmer_starts);
                    }
                }
                else
                {
                    processReadsByLine(in_fastq_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
                }
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
KmerMatcher::ProcessReadResult 
Typer::processRead(TSequence const &           read,
                   uint64_t                    k,
                   Finder<TSequenceEsaIndex> & finder_esa,
                   bool                        save_kmer_starts)
{
    // calledLocalProcessRead = true;
    ProcessReadResult result = KmerMatcher::UNPROCESSED_READ;
    // look for middle k-mer
    uint64_t read_length = length(read);
    uint64_t start_pos     = (read_length - k) / 2;
    TSequence middle_kamer = infix(read, start_pos, (start_pos + k + 1));
    // cerr << "DEBUG: processRead - Typer - read: "<< read << endl;
    clear(finder_esa);
    if (find(finder_esa, middle_kamer))   // middle k-mer found
    {
        // ++total_reads;
        result = KmerMatcher::PROCESSED_READ;
        clear(finder_esa);
        uint64_t stop_base = read_length - k;
        uint64_t allele_idx;
        uint64_t hit_base;

        for (uint64_t j = 0; j <= stop_base; ++j) // k-merize all the read
        {
            TSequence kmer       = infix(read, j, j + k);
            uint64_t pos_counter = 0;
            while (find(finder_esa, kmer)) // search each k-mer in the index
            {
                TSequenceEsaIndexPosition pos = position(finder_esa);
                
                if (save_kmer_starts)
                {
                    allele_idx = reduced_index_allele_indices.at(pos.i1);
                    hit_base   = pos.i2;

                    // Add the start of the hit in the corresponding allele
                    if (kmer_starts.find(allele_idx) == kmer_starts.end())
                        kmer_starts[allele_idx] = {};

                    kmer_starts[allele_idx].push_back(hit_base);
                    
                    // // Temp
                    // ++allele_counts[allele_idx];
                    // ++total_hits;
                    // ++pos_counter;
                    // // Temp
                    
                }
                else 
                {
                    ++allele_counts[pos.i1];
                    ++total_hits;
                    ++pos_counter;
                }
            }

            if (pos_counter > 0) ++total_kmers;

            clear(finder_esa);
        }
    }
    return result;
}
// // ----------------------------------------------------------------------------
// KmerMatcher::ProcessReadResult 
// Typer::processRead(TSequence const &           read,
//                    uint64_t                    k,
//                    Finder<TSequenceEsaIndex> & finder_esa,
//                    bool                        save_kmer_starts)
// {
//     ProcessReadResult result = KmerMatcher::UNPROCESSED_READ;
//     // look for middle kmer
//     uint64_t start_pos     = (length(read) - k) / 2;
//     TSequence middle_kamer = infix(read, start_pos, (start_pos + k + 1));
//     clear(finder_esa);

//     if (find(finder_esa, middle_kamer))   // middle k-mer found
//     {
//         // ++total_reads;
//         result = KmerMatcher::PROCESSED_READ;
//         clear(finder_esa);

//         for (uint64_t j = 0; j <= length(read) - k; ++j) // k-merize all the read
//         {
//             TSequence kmer       = infix(read, j, j + k);
//             uint64_t pos_counter = 0;
//             while (find(finder_esa, kmer)) // search each k-mer in the index
//             {
//                 TSequenceEsaIndexPosition pos = position(finder_esa);
//                 ++allele_counts[pos.i1];
//                 ++total_hits;
//                 ++pos_counter;

//                 if (save_kmer_starts)
//                 {
//                     uint64_t allele_idx = pos.i1;
//                     uint64_t hit_base   = pos.i2;

//                     // Store the start of the k-mer hit if the match
//                     // was in one of the best alleles
//                     if (find(first_max_allele.begin(), first_max_allele.end(), allele_idx) != first_max_allele.end())
//                     {
//                         // Add the start of the hit in the corresponding allele
//                         if (kmer_starts.find(allele_idx) == kmer_starts.end())
//                             kmer_starts[allele_idx] = {};

//                         kmer_starts[allele_idx].push_back(hit_base);
//                     }
//                 }
//                 else 
//                 {
//                     ++allele_counts[pos.i1];
//                     ++total_hits;
//                     ++pos_counter;
//                 }
//             }

//             if (pos_counter > 0) ++total_kmers;

//             clear(finder_esa);
//         }
//     }
//     return result;
// }
// ----------------------------------------------------------------------------
void Typer::selectBestAlleles()
{
    for(auto&& locus : loci_table) 
    {
        SAlleleCount first_allele;
        SequenceInfo best_allele;
        best_allele.setAlleleNumber(NA_STRING);
        
        string locus_id    = locus.id;
        uint64_t n_alleles = 0;
        
        if (top_alleles.find(locus.id) != top_alleles.end())
        {
            n_alleles = top_alleles.at(locus_id).size();
        }
        

        if (n_alleles > 0) // there are k-mer hits for the first top allele
        {
            first_allele  = top_alleles.at(locus_id).at(0);
            
            if (first_allele.count > 0) // there are k-mer hits for the first top allele
            {
                uint64_t allele_idx = first_allele.idx;
                string allele_desc = toCString(allele_seq_ids[allele_idx]);
                TStringVector fields;
                string allele_number;
                split(fields, allele_desc, ALLELE_SEQ_ID_SEPARATOR);
                if (fields.size() >= 2)
                {
                    allele_number = fields.at(1);
                }
                else
                {
                    cerr << "ERROR: Bad format of the allele ids in a locus sequence file. "
                            "Each sequence id must be separated by '" << ALLELE_SEQ_ID_SEPARATOR << "' (e.g. >abcZ" << ALLELE_SEQ_ID_SEPARATOR << "34)" << endl;
                    exit(1);
                }
                
                best_allele.setIdx(allele_idx);
                best_allele.setDesc(allele_desc);
                best_allele.setAlleleNumber(allele_number);
                best_allele.setKmerCount(allele_counts.at(allele_idx));
                best_allele.setNormKmerCount(first_allele.count);
                best_allele.setLength(length(indexText(alleles_esa_index)[allele_idx]));
                
                
                best_max_allele_numbers.push_back(allele_number);
            }
        }
        else 
        {
            // best_max_allele.push_back(-1);
            // cout << "DEBUG:   aquI " << endl;
            best_max_allele_numbers.push_back(NA_STRING);
            // best_max_norm_kmer_counts.push_back(0);
        }
        
        // Add current best allele to the best allelic profile map
        best_allelic_profile[locus_id] = best_allele;
    }
}
// ----------------------------------------------------------------------------
void Typer::getTopNormAlleles(uint64_t n_top)
{
    uint64_t start = 0;
    for(auto&& locus : loci_table) 
    {
        uint64_t end            = locus.last_seq_idx;
        // double first_max_count  = -1;
        // int first_max_idx       = -1;
        
        // create a empty vector for storing top alleles of the current locus
        if (top_alleles.find(locus.id) == top_alleles.end()) 
            top_alleles[locus.id] = {};
        
        // loop through the allele counts that correspond to the current locus
        for (uint64_t j = start; j <= end; ++j)
        {
            // create a new allele count instance (struct) for the current allele
            SAlleleCount current_alele_count;
            
            if (norm_allele_counts[j] > 0) 
            {
                current_alele_count.idx   = (int) j;
                current_alele_count.count = norm_allele_counts[j];
                // add the current allele count 
                top_alleles.at(locus.id).push_back(current_alele_count);
                top_alleles_indices.push_back(j);
            }
            
        }
        
        // If there were some alleles with counts
        uint64_t n_alleles_with_counts = top_alleles.at(locus.id).size();
        if (n_alleles_with_counts > 0)
        {
            // Sort the vector of allele counts for the current locus
            std::sort(top_alleles.at(locus.id).begin(), top_alleles.at(locus.id).end(), 
                greaterSAlleleCount);
            // cout << "--- " << locus.id << " --- " << endl;
            // for (auto &&x: top_alleles[locus.id])
            //     cout << x.idx << "\t" << x.count << endl;
            
            // cout << "---" << endl;
            // for (auto &&x: top_alleles[locus.id])
            //     cout << x.idx << "\t" << x.count << endl;
            // cout << "DEBUG:     before resize" << endl;
            // cout << "n_top:" << n_top << endl;
            uint64_t new_size = (n_alleles_with_counts < n_top) ? n_alleles_with_counts : n_top;
            top_alleles.at(locus.id).resize(new_size);
        }
        
        start = end + 1;
    }
}
// // ----------------------------------------------------------------------------
// void Typer::computeKmerDepthAndAlleleCoverageNew()
// {
//     vector<vector<int>> max_alleles; // = {first_max_allele, second_max_allele};
//     max_alleles.push_back(first_max_allele);
//     // max_alleles.push_back(second_max_allele);

//     // for(auto profile : max_alleles) {
//     for (uint64_t i = 0; i < max_alleles.size(); ++i)
//     {
//         vector<TUintVector> allele_kmer_depths;
//         TDoubleVector allele_mean_kmer_depths;
//         TSeqCovVector allele_coverages;
//         // uint64_t n_kmer_depth_saved = 0;
//         // for (auto&& allele_idx : max_alleles[i])
//         for (unsigned j = 0; j < max_alleles.at(i).size(); ++j)
//         {
//             int allele_idx = max_alleles.at(i).at(j);

//             SeqCoverage allele_cov;
//             TUintVector per_base_depths;
//             double mean_depth = 0.0;
//             if (allele_idx > -1) // is not an NA allele
//             {
//                 // Sort the kmer hit starts of the current allele
//                 std::sort(kmer_starts[allele_idx].begin(), kmer_starts[allele_idx].end());

//                 // // Construct the intervals of the k-mer hits and store them in a vector
//                 // TIntervalVector intervals;
//                 // for (unsigned j = 0; j < kmer_starts[allele_idx].size(); ++j)
//                 // {
//                 //     uint64_t start = kmer_starts[allele_idx].at(j);
//                 //     intervals.push_back(TInterval(start, start + options.getKmerLength(), allele_idx));
//                 // }
//                 // // Construct the interval tree to calculate depths
//                 // IntervalTree<uint64_t, uint64_t> tree(intervals);
//                 // Initialize variables
//                 uint64_t allele_length     = length(indexText(alleles_esa_index)[allele_idx]);
//                 uint64_t depth_sum         = 0;
//                 uint64_t total_bases       = 0;
//                 uint64_t allele_starts_size = kmer_starts[allele_idx].size();
//                 // cout << "DEBUG " << allele_idx << " " << allele_starts_size << endl;
//                 bool check_for_depth_ratio = true;
//                 per_base_depths.resize(allele_length, 0);
                
//                 // Loop through k-mer starts and increment depth
//                 for (unsigned startIdx = 0; startIdx < allele_starts_size; ++startIdx) 
//                 {
//                     uint64_t start = kmer_starts[allele_idx].at(startIdx);
//                     // cout << "DEBUG " << allele_idx << " " << allele_starts_size << " " << startIdx << " " << start <<  kmer_starts.size() << endl;
                    
//                     for(unsigned base = start; base < (start + options.getKmerLength()); ++base)
//                         ++per_base_depths.at(base);
//                 }

//                 // Loop through the bases of the current allele to calculate per-base depth and allele coverage
//                 for (uint64_t base = 0; base < allele_length; ++base)
//                 {
//                     if (per_base_depths.at(base) == 0) check_for_depth_ratio = false;

//                     // Look for depth gaps
//                     if (check_for_depth_ratio &&
//                             (base >= BASES_MARGIN - 1) && (base <= allele_length - BASES_MARGIN - 1))
//                     {
//                         double depth_ratio = (double)per_base_depths.at(base) / (double)per_base_depths.at(base - 1);

//                         if (depth_ratio > DEPTH_DIFF_MAX_PEAK || depth_ratio < DEPTH_DIFF_MIN_PEAK)
//                         {
//                             per_base_depths.at(base) = 0;
//                             check_for_depth_ratio    = false;
//                         }
//                     }

//                     // per_base_depths.push_back(depth);
//                     depth_sum += per_base_depths.at(base);

//                     if (per_base_depths.at(base) > 0)
//                         ++total_bases;
//                 }
//                 // Calculate mean k-mer depth for the current allele
//                 mean_depth = (double)depth_sum / (double)allele_length;

//                 // Calculate coverage for the current allele
//                 double coverage      = 100 * ((double)total_bases / (double)allele_length);
//                 allele_cov.cov_bases = total_bases;
//                 allele_cov.seq_size  = allele_length;
//                 allele_cov.coverage  = coverage;

//                 // Save the current allele per-base k-mer depth data to a file
//                 if (options.getKmerDepthOutFile() != "")
//                 {
//                     string out_filename = options.getKmerDepthOutFile();

//                     // bool append_depth_data = (n_kmer_depth_saved > 0);
//                     bool append_depth_data = (j > 0);

//                     if (saveKmerDepthData(options.getKmerDepthOutFile(),
//                                           per_base_depths,
//                                           allele_idx,
//                                           append_depth_data) == SAVE_ERROR)
//                     {
//                         exit(1);
//                     }
//                     else
//                     {
//                         // ++n_kmer_depth_saved;

//                         if (j == max_alleles.at(i).size() - 1) {
//                             string message = "Per-base k-mer depth data saved to '" + out_filename + "'.";
//                             printMessage(message);
//                         }
//                     }
//                 }
//             } // End if NA allele
//             allele_kmer_depths.push_back(per_base_depths);
//             allele_mean_kmer_depths.push_back(mean_depth);
//             allele_coverages.push_back(allele_cov);
//         } // End alleles loop
//         perbase_kmer_depths.push_back(allele_kmer_depths);
//         mean_kmer_depths.push_back(allele_mean_kmer_depths);
//         coverages.push_back(allele_coverages);
//     } // End profiles loop
// }
// ----------------------------------------------------------------------------
void Typer::computeKmerDepthAndAlleleCoverageForTopAlleles()
{
    // cout << "DEBUG: computeKmerDepthAndAlleleCoverageForTopAlleles - top_alleles.size(): " << top_alleles.size() << endl;
    for (auto const& top_allele : top_alleles)
    {
        string locus_id            = top_allele.first;
        TAlleleCountVector alleles = top_allele.second;
        // cout << "top_alleles.size(): " << top_alleles.size() << endl;// DEBUG
        if (top_alleles_info.find(locus_id) == top_alleles_info.end())
            top_alleles_info[locus_id] = {};
        
        // cout << "DEBUG: " << locus_id << endl;
        uint64_t n_top = alleles.size();
        for (uint64_t i = 0; i < n_top; ++i)
        {
            // cout << "  DEBUG: alleles.size(): " << alleles.size() << endl;
            int allele_idx = alleles.at(i).idx;
            
            SequenceInfo current_allele;
            // current_allele.idx        = allele_idx;
            current_allele.setIdx(allele_idx);
            // current_allele.norm_count = alleles.at(i).count;
            current_allele.setNormKmerCount(alleles.at(i).count);
            // cout << "DEBUG: " << alleles.at(i).idx << "\t" << alleles.at(i).count << endl;
            // cout << "DEBUG: i=" << i << "  " << toCString(allele_seq_ids[allele_idx]) << "  current_allele.norm_count: " << current_allele.norm_count << endl;
            if (current_allele.getNormKmerCount() > 0) // there are kmer hits for the corresponding locus
            {   
                uint64_t allele_length      = length(indexText(alleles_esa_index)[allele_idx]);
                TUintVector per_base_depths = vector<uint64_t>(allele_length, 0);
                uint64_t depth_sum          = 0;
                uint64_t total_bases        = 0;
                bool check_for_depth_ratio  = true;
                
                // cerr << "DEBUG - kmer_starts.at(allele_idx).size(): " << kmer_starts.at(allele_idx).size() << endl;
                // Sort the kmer hit starts of the current allele
                std::sort(kmer_starts[allele_idx].begin(), kmer_starts[allele_idx].end());
                
                // Loop through kmer starts and increment depth
                for (unsigned j = 0; j < kmer_starts[allele_idx].size(); ++j) 
                {
                    uint64_t start = kmer_starts[allele_idx].at(j);
                    
                    for(unsigned base = start; base < (start + options.getKmerLength()); ++base)
                        ++per_base_depths.at(base);
                        // ++current_allele.per_base_depths.at(base);
                }
                
                // Look for depth gaps
                for(unsigned base = 0; base < allele_length; ++base) 
                {
                    uint64_t depth = per_base_depths.at(base);
                    
                    if (depth == 0) 
                        check_for_depth_ratio = false;
                    
                    if (check_for_depth_ratio && 
                        (base >= BASES_MARGIN - 1) && (base <= allele_length - BASES_MARGIN - 1))
                    {
                        double depth_ratio = (double)depth/(double)per_base_depths.at(base-1);
                        
                        if (depth_ratio > DEPTH_DIFF_MAX_PEAK || depth_ratio < DEPTH_DIFF_MIN_PEAK)
                        {
                            per_base_depths.at(base) = 0;
                            check_for_depth_ratio    = false;
                        }
                    }
                    
                    // per_base_depths.push_back(depth);
                    depth_sum += per_base_depths.at(base);
                    
                    if (per_base_depths.at(base) > 0)
                        ++total_bases;
                }
                current_allele.setDesc(toCString(allele_seq_ids[allele_idx]));
                current_allele.setPerBaseDepths(per_base_depths);
                // Calculate mean k-mer depth for the current allele
                current_allele.setAvgKmerDepth((double)depth_sum / (double)allele_length);

                // Calculate coverage for the current allele
                current_allele.setNBasesCovered(total_bases);
                current_allele.setLength(allele_length);
                // current_allele.coverage  = 100 * ((double)total_bases / (double)allele_length);
                
                // cout << "DEBUG: i=" << i << "  " << toCString(allele_seq_ids[allele_idx]) << "   current_allele.coverage: " << current_allele.coverage << endl;
                // cout << "DEBUG: i=" << i << "  " << toCString(allele_seq_ids[allele_idx]) << "   current_allele.mean_kmer_depth: " << current_allele.mean_kmer_depth << endl << endl;
                // top_alleles_info[locus_id].push_back(current_allele);
                top_alleles_info[locus_id].push_back(current_allele);
            } // end if current_allele.norm_count
            
            // cerr << "DEBUG - current_allele: \n" << current_allele << endl;
        } // end for loop n_top alleles
    } // end for loop top alleles (for each loci)
}
// ----------------------------------------------------------------------------
string Typer::createProfileStringRep(TStringVector max_allele_numbers)
{
    // construct the profile (max alleles) string representation
    string inferred_profile_string = "*";
    string sep                     = "-";
    for (uint64_t j = 0; j < (max_allele_numbers.size()); ++j)
    {
        inferred_profile_string += max_allele_numbers[j] + sep;
    }

    return inferred_profile_string;
}
// ----------------------------------------------------------------------------
void Typer::loadProfilesEsaIndex()
{
    string profiles_index_filename = options.getIndexFilename() + PROFILE_INDEX_SUFFIX;
    string esa_index_name          = profiles_index_filename + ESA_INDEX_EXT;

    if (checkFile(esa_index_name) == FILE_NOT_FOUND)
    {
        cerr << "ERROR: Index file '" << esa_index_name << "' not found." << endl;
        exit(1);
    }

    if (!open(profiles_index, profiles_index_filename.c_str()))
    {
        cerr << "ERROR: Could not open the allelic profiles index '" << profiles_index_filename << "'\n";
        exit(1);
    }
}
// ----------------------------------------------------------------------------
void Typer::loadProfilesFile()
{
    string filename = options.getIndexFilename() + ALLELIC_PROFILE_EXT;

    unordered_map<string, vector<string>> table;
    fstream profiles_file(filename, ios::in);

    if (profiles_file.is_open())
    {
        string line;
        vector<string> col_names;
        uint64_t current_line = 0;
        while (getline(profiles_file, line))
        {
            vector<string> fields;

            if (!line.empty() && line.at(0) != COMMENT_CHAR)    // The current line is not an empty/comment line
            {
                split(fields, line, FIELD_SEPARATOR);

                if (++current_line == 1)    // get header and create the map
                {
                    profiles_table.cols = fields;
                    profiles_table.cols.at(0) = "ST";
                    for (string col : profiles_table.cols) {
                        table[col] = vector<string>();
                        profiles_table.table[col] = table[col];
                    }
                }
                else {
                    for (uint64_t i = 0; i < profiles_table.cols.size(); ++i) {
                        string col = profiles_table.cols[i];
                        profiles_table.table[col].push_back(fields[i]);
                    }
                }
            }
        }
        profiles_file.close();
    }
    else {
        cerr << "ERROR: Could not open the profiles file '" << filename << "'\n";
        exit(1);
    }
}
// ----------------------------------------------------------------------------
void Typer::selectBestAllelesBasedOnCoverage()
{
    // best_max_allele.clear();
    // best_max_norm_kmer_counts.clear();
    best_max_allele_numbers.clear();
    best_allelic_profile.clear();
    
    for(auto&& locus : loci_table) 
    {
        SequenceInfo best_allele;
        string locus_id        = locus.id;
        TSeqInfoVector alleles = top_alleles_info.at(locus_id);
        // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << endl;
        // cerr << "DEBUG: " << endl;
        // Update allele numbers
        uint64_t n_alleles = alleles.size();
        // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << "  - n. alleles: " << n_alleles << endl;
        // for(auto&& allele : alleles)
        // {
        //     TStringVector fields;
        //     split(fields, toCString(allele_seq_ids[allele.getIdx()]), ALLELE_SEQ_ID_SEPARATOR);
        //     if (fields.size() >= 2)
        //     {
        //         allele.setAlleleNumber(fields.at(1));
        //     }
        //     else
        //     {
        //         cerr << "ERROR: Bad format of the allele ids in a locus sequence file. "
        //                 "Each sequence id must be separated by '" << ALLELE_SEQ_ID_SEPARATOR << "' (e.g. >abcZ" << ALLELE_SEQ_ID_SEPARATOR << "34)" << endl;
        //         exit(1);
        //     }
        // }
        if (n_alleles > 0 && alleles.at(0).getNormKmerCount() > 0) // there are kmer hits for the first top allele
        {
            // cout << "DEBUG: locus_id: " << locus_id << endl;
            // cerr << "DEBUG: Before sorting by Score" << endl;
            // printVector(alleles, "\n", cerr);
            
            // // Sort top alleles by coverage
            // std::sort(alleles.begin(), alleles.end(), SSequenceInfoGreater(SequenceInfo::COVERAGE));
            // cerr << "DEBUG: Before sorting by Coverage" << endl;
            // cerr << "DEBUG: Before sorting by Count, Cov and sd" << endl;
            // printVector(alleles, "\n", cerr);
            // uint64_t count = 0;
            // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << " - Before updating allele info" << endl;
            for(auto&& allele : alleles) {
                
                // Update allele numbers
                TStringVector fields;
                // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - allele: " << allele.getIdx() << " norm count: " << allele.getNormKmerCount() << " Before split " << ++count << endl;
                string allele_desc = toCString(allele_seq_ids[allele.getIdx()]);
                split(fields, allele_desc, ALLELE_SEQ_ID_SEPARATOR);
                // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - allele: " << allele_desc << " After split " << count << endl;
                if (fields.size() >= 2)
                {
                    allele.setAlleleNumber(fields.at(1));
                }
                else
                {
                    cerr << "ERROR: Bad format of the allele ids in a locus sequence file. "
                            "Each sequence id must be separated by '" << ALLELE_SEQ_ID_SEPARATOR << "' (e.g. >abcZ" << ALLELE_SEQ_ID_SEPARATOR << "34)" << endl;
                    exit(1);
                }
                
                // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << " - allele: " << allele.getAlleleNumber() << " - Before calculating SD" << endl;
                allele.calculateKmerDepthSd();
                // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << " - allele: " << allele.getAlleleNumber() << " - After calculating SD" << endl;
                // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << " - allele: " << allele.getAlleleNumber() << " - Before calculating score" << endl;
                allele.calculateScore();
                // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << " - allele: " << allele.getAlleleNumber() << " - After calculating score" << endl;
            }
            // cerr << "DEBUG: selectBestAllelesBasedOnCoverage - locus: " << locus_id << " - After updating allele info" << endl;
            
            // Sort top alleles by score
            // cerr << "DEBUG: After sorting by Count, Cov and sd" << endl;
            // std::sort(alleles.begin(), alleles.end(), SSequenceInfoGreater(SequenceInfo::SCORE));
            std::sort(alleles.begin(), alleles.end(), SSequenceInfoSort());
            // printVector(alleles, "\n", cerr);
            
            // TStringVector fields;
            // vector<uint64_t> tied_alleles_indices = {0}; // tied alleles by coverage
            // // SSeqInfo first_allele = alleles.at(0);
            // SequenceInfo first_allele = alleles.at(0);
            // // check for tied alleles by coverage
            // if (n_alleles > 1)
            // {
            //     for(unsigned i = 1; i < n_alleles; ++i)
            //     {
            //         if (first_allele.getCoverage() == alleles.at(i).getCoverage()) 
            //         {
            //             tied_alleles_indices.push_back(i);
            //         }
            //     }
            // }
            
            // cout << "DEBUG: indices of tied alleles: " << locus_id << endl;
            // printVector(tied_alleles_indices, ", ", cerr);
            // // calculate k-mer depth sd if there are tied alleles by coverage
            // if (tied_alleles_indices.size() > 1)
            // {
            //     for(auto const& i : tied_alleles_indices)
            //     {
            //         alleles.at(i).calculateKmerDepthSd();
            //     }
            //     // Sort top alleles by depth sd
            //     std::sort(alleles.begin(), alleles.end(), SSequenceInfoLess(SequenceInfo::DEPTH_SD));
            // }
            
            // Update arrays of values from the top allele
            // best_max_allele.push_back(alleles.at(0).getIdx());
            best_max_allele_numbers.push_back(alleles.at(0).getAlleleNumber());
            // best_max_norm_kmer_counts.push_back(alleles.at(0).getNormKmerCount());
            
            // Take the top allele and copy its information to the best allelic profile map  
            best_allele = alleles.at(0);
        }
        else 
        {
            // best_max_allele.push_back(-1);
            // best_max_norm_kmer_counts.push_back(0);
            best_max_allele_numbers.push_back(NA_STRING);
        }
        
        // Update alleles (with the correct order) for the current locus
        top_alleles_info.at(locus_id)  = alleles;
        // Save the new best allele at the best allelic profile map
        best_allelic_profile[locus_id] = best_allele;
    }
}
// ----------------------------------------------------------------------------
Typer::StStatus Typer::getInferredSt(string &        inferred_st,
                                     const string &  profile_string)
{
    StStatus status;
    // Search the profile (max alleles) in the allelic profile index
    uint64_t profile_idx_found;
    Finder<TCharEsaIndex> finder(profiles_index);

    if (find(finder, profile_string))
    {
        TSequenceEsaIndexPosition pos = position(finder);
        profile_idx_found             = pos.i1;

        // Get the predicted ST
        string st_col_name = profiles_table.cols.at(0);
        inferred_st        = profiles_table.table[st_col_name][profile_idx_found];
        status             = Typer::ST_PREDICTED;
    }
    else {
        inferred_st = NA_STRING;
        if (profile_string.find(NA_STRING) != std::string::npos)
            status = Typer::NO_KMER_HITS;
        else
            status = Typer::NO_ST_IN_TABLE;
    }

    return status;
}
// ----------------------------------------------------------------------------
SaveFileResult Typer::saveTopAllelesPerBaseKmerDepthData(const string & output_file_name)
{
    fstream out_file;
    out_file.open(output_file_name, ios::trunc |  ios::out);
    
    if (!out_file)
    {
        cerr << "ERROR: Per-base k-mer depth file could not be saved." << endl;
        return SAVE_ERROR;
    }
    
    // Write header
    out_file << "n"          << FIELD_SEPARATOR
             << "locus"      << FIELD_SEPARATOR
             << "allele"     << FIELD_SEPARATOR
             << "position"   << FIELD_SEPARATOR
             << "kmer_depth" << FIELD_SEPARATOR
             << "norm_count" << FIELD_SEPARATOR
             << "coverage"   << FIELD_SEPARATOR
             << "score"      << FIELD_SEPARATOR
             << "mean_kmer_depth"
             // << "weighted_kmer_depth" 
             << endl;
    
    // cout << "DEBUG: top_alleles_info.size(): " << top_alleles_info.size() << endl;
    
    for(auto&& locus : loci_table) 
    // for (auto const& top_allele : top_alleles_info)
    {
        // string locus_id        = top_allele.first;
        // TSeqInfoVector alleles = top_allele.second;
        string locus_id        = locus.id;
        TSeqInfoVector alleles = top_alleles_info.at(locus_id);
        
        // cout << "DEBUG - Writing k-mer depth data: " << locus_id << endl;
        // cout << "  DEBUG: alleles.size(): " << alleles.size() << endl;
        for (uint64_t i = 0; i < alleles.size(); ++i)
        {
            // int allele_idx   = alleles.at(i).getIdx();
            // string allele_id = toCString(allele_seq_ids[allele_idx]);
            string allele_id = alleles.at(i).getDesc();
            
            if (alleles.at(i).getNormKmerCount() > 0) // there are kmer hits for the corresponding locus
            {
                vector<uint64_t> per_base_depths = alleles.at(i).getPerBaseDepths();
                uint64_t per_base_depths_size    = per_base_depths.size();
                double norm_kmer_count           = alleles.at(i).getNormKmerCount();
                double coverage                  = alleles.at(i).getCoverage();
                double score                     = alleles.at(i).getScore();
                double avg_depth                 = alleles.at(i).getAvgKmerDepth();
                for(uint64_t j = 0; j < per_base_depths_size; ++j)
                {
                    out_file << i+1
                             << FIELD_SEPARATOR
                             << locus_id
                             << FIELD_SEPARATOR
                             << allele_id
                             << FIELD_SEPARATOR
                             << j+1
                             << FIELD_SEPARATOR
                             << per_base_depths.at(j)
                             << FIELD_SEPARATOR
                             << norm_kmer_count
                             << FIELD_SEPARATOR
                             << coverage
                             << FIELD_SEPARATOR
                             << score
                             << FIELD_SEPARATOR
                             << avg_depth
                             << endl;
                }
            }
        }
    }
    
    out_file.close();
    return SAVE_SUCCESS;
}
// ----------------------------------------------------------------------------
void Typer::updateBestAllelicProfileCovData()
{
    for(auto&& locus : loci_table) 
    {
        // cerr << "DEBUG: updateBestAllelicProfileCovData antes "  << endl;
        TSeqInfoVector alleles = top_alleles_info.at(locus.id);
        if (alleles.size() > 0) {
            uint64_t n_alleles = 1;
            for (uint64_t i = 0; i < n_alleles; ++i)
            {
                // best_allelic_profile.at(locus.id).setLength(alleles.at(i).getLength());
                best_allelic_profile.at(locus.id).setNBasesCovered(alleles.at(i).getNBasesCovered());
                best_allelic_profile.at(locus.id).setAvgKmerDepth(alleles.at(i).getAvgKmerDepth());
                best_allelic_profile.at(locus.id).calculateKmerDepthSd();
                // TUintVector per_base_depths = alleles.at(i).getPerBaseDepths();
                // best_allelic_profile.at(locus.id).setPerBaseDepths(per_base_depths);
            }
        }
    }
}
// ----------------------------------------------------------------------------
string Typer::getInfoMessage()
{
    string margin = "  ";
    string info;
    info  = "\n\x1B[1mSTing TYPER OUTPUT DETAILS\x1B[0m\n";
    info += "    \x1B[1mStatus:\n";
    info += "      st_predicted:\x1B[0m\n"
            "        ST inferred from k-mer matches and profiles table.\n\n";
    info += "      \x1B[1mNO_KMER_HITS:\x1B[0m\n"
            "        There are no k-mer matches for one or more of the loci to infer\n"
            "        a profile and its associated ST. Probable causes: \n"
            "         1) k-mer size not adequate (too long),\n"
            "         2) low quality data/too many N's in the data.\n\n";
    info += "      \x1B[1mno_st_in_table:\x1B[0m\n"
            "        There is no ST associated to the inferred alleleic profile.\n"
            "        Probable causes: \n"
            "         1) k-mer size not adequate (too short),\n"
            "         2) this could be a new allelic profile.\n\n";
    info += "    \x1B[1m" + string(NA_STRING) + ":\x1B[0m\n"
            "      No ST associated or no k-mer matches for any allele of a locus.\n";
    info += "    \x1B[1m*:\x1B[0m\n"
            "      Allele predicted with less than 100% of length coverage (low confidence).\n"
            "      Probable causes:\n"
            "        1) No enough k-mer hits to cover the whole lenght of the allele.\n"
            "        2) No enough k-mer depth in some part along the allele ([10bp, lenght-10bp]\n"
            "           bases range) to consider it as covered.\n";
    // info += "===============================================================================\n";

    return info;
}
// ----------------------------------------------------------------------------
template<typename T>
string Typer::getResultsLine(const string &    line_type,
                             const StStatus &  status,
                             const vector<T> & info_vector)
{
    stringstream out_string;
    if (options.getSampleName() != "")
        out_string << options.getSampleName() << FIELD_SEPARATOR;
    
    out_string << line_type
               << FIELD_SEPARATOR
               << getStatusString(status)
               << FIELD_SEPARATOR
               << inferred_st
               << FIELD_SEPARATOR;
    printVector(info_vector, FIELD_SEPARATOR, out_string);
    out_string << FIELD_SEPARATOR
               << total_kmers
               << FIELD_SEPARATOR
               << total_reads
               << FIELD_SEPARATOR
               << fastqFilesToString()
               << endl;
    
    return out_string.str();
}
// ----------------------------------------------------------------------------
template<typename T>
void Typer::printResultsLine(const string &    line_type,
                             const StStatus &  status,
                             const vector<T> & info_vector)
{
    if (options.getSampleName() != "")
        cout << options.getSampleName() << FIELD_SEPARATOR;
    
    cout << line_type
         << FIELD_SEPARATOR
         << getStatusString(status)
         << FIELD_SEPARATOR
         << inferred_st
         << FIELD_SEPARATOR;
    printVector(info_vector, FIELD_SEPARATOR, std::cout);
    cout << FIELD_SEPARATOR
         << total_kmers
         << FIELD_SEPARATOR
         << total_reads
         << FIELD_SEPARATOR
         << fastqFilesToString()
         << endl;
}
// ----------------------------------------------------------------------------
string Typer::getResultsHeader()
{
    // string out_string = "";
    stringstream out_string;
    // Get the loci ids list from the loci_table
    TStringVector loci_ids;
    for (auto && locus : loci_table) {
        loci_ids.push_back(locus.id);
    }

    // Prepare ST and Allelic profile header
    TStringVector header;
    if (options.getSampleName() != "")
        header.push_back("Sample");
    header.push_back("Line_type");
    header.push_back("Status");
    header.push_back("ST");
    header.insert(header.end(), loci_ids.begin(), loci_ids.end());
    header.push_back("Total_k-mers");
    header.push_back("Total_reads");
    header.push_back("Input_files");

    // print header
    printVector(header, FIELD_SEPARATOR, out_string);
    out_string << endl;
    
    return out_string.str();
}
// ----------------------------------------------------------------------------
string Typer::getTidyResultsHeader()
{
    stringstream out_string;
    // Prepare header
    TStringVector header;
    if (options.getSampleName() != "")
        header.push_back("Sample");
    header.push_back("Status");
    header.push_back("ST");
    header.push_back("Locus");
    header.push_back("Predicted_allele");
    if (options.getPrintKmerCounts())
        header.push_back("Norm_k-mer_counts");
    if (options.getPrintAlleleCov())
    {
        header.push_back("Allele_cov(%)");
        header.push_back("Covered_bases");
        header.push_back("Allele_length");
    }
    if (options.getPrintKmerDepth())
        header.push_back("Avg_k-mer_depth");
    header.push_back("Total_k-mers");
    header.push_back("Total_reads");
    header.push_back("Input_files");

    // print header
    printVector(header, FIELD_SEPARATOR, out_string);
    out_string << endl;
    
    return out_string.str();
}
// ----------------------------------------------------------------------------
string Typer::getResultsString(const StStatus & status)
{
    string out_string = "";
    
    // Construct and populate vectors for printing 
    vector<string> best_inferred_profile;
    vector<string> best_norm_kmer_counts;
    vector<string> best_coverages;
    vector<string> best_avg_kmer_depths;
     
    for (auto && locus : loci_table)
    {
        SequenceInfo allele = best_allelic_profile.at(locus.id);
        best_inferred_profile.push_back(allele.getAlleleNumber());
        best_norm_kmer_counts.push_back(getNumberStringWithPrecision(allele.getNormKmerCount(),1));
        best_coverages.push_back(allele.getCoverageStr());
        best_avg_kmer_depths.push_back(getNumberStringWithPrecision(allele.getAvgKmerDepth(),1));
    }
    
    // Print results to a string
    out_string += getResultsLine("allelic_profile", status, best_inferred_profile);
    if (options.getPrintKmerCounts())
        out_string += getResultsLine("norm_kmer_counts", status, best_norm_kmer_counts);
    if (options.getPrintAlleleCov())
        out_string += getResultsLine("allele_coverage", status, best_coverages);
    if (options.getPrintKmerDepth())
        out_string += getResultsLine("avg_kmer_depth", status, best_avg_kmer_depths);
    
    return out_string;
}
// ----------------------------------------------------------------------------
void Typer::printResults(const StStatus & status)
{
    // bool save_to_file = true;
    string output_filename = options.getOutputFile();
    
    if (output_filename != "") {
        fstream out_file;
        out_file.open(output_filename, ios::trunc |  ios::out);
        
        if (!out_file)
        {
            cerr << "ERROR: Output file '" << output_filename << "'' could not be saved." << endl;
            cerr << getResultsHeader();
            cout << getResultsString(status);
        }
        else {
            out_file << getResultsHeader();
            out_file << getResultsString(status);
        }
    }
    else 
    {
        cerr << getResultsHeader();
        cout << getResultsString(status);
    }
    
}

// ----------------------------------------------------------------------------
string Typer::getTidyResultsString(const StStatus & status)
{
    stringstream out_string;
    
    for (auto && locus : loci_table)
    {
        SequenceInfo allele = best_allelic_profile.at(locus.id);
        
        if (options.getSampleName() != "")
            out_string << options.getSampleName() << FIELD_SEPARATOR;
        out_string << getStatusString(status)
                   << FIELD_SEPARATOR
                   << inferred_st
                   // << predicted_st
                   << FIELD_SEPARATOR
                   << locus.id
                   << FIELD_SEPARATOR
                   << allele.getAlleleNumber()
                   << FIELD_SEPARATOR;
        if (options.getPrintKmerCounts())
        {
            out_string << getNumberStringWithPrecision(allele.getNormKmerCount(), 1)
                       << FIELD_SEPARATOR;
        }
        if (options.getPrintAlleleCov()) {
            out_string << getNumberStringWithPrecision(allele.getCoverage(), 1)
                       << FIELD_SEPARATOR
                       << allele.getNBasesCovered()
                       << FIELD_SEPARATOR
                       << allele.getLength()
                       << FIELD_SEPARATOR;
        }
        if (options.getPrintKmerDepth())
        {
            out_string << getNumberStringWithPrecision(allele.getAvgKmerDepth(), 1)
                       << FIELD_SEPARATOR;
        }
        out_string << total_kmers
                   << FIELD_SEPARATOR
                   << total_reads
                   << FIELD_SEPARATOR
                   << fastqFilesToString()
                   << endl;
    }
    
    return out_string.str();
}
// ----------------------------------------------------------------------------
void Typer::printTidyResults(const StStatus & status)
{
    string output_filename = options.getOutputFile();
    
    if (output_filename != "") {
        fstream out_file;
        out_file.open(output_filename, ios::trunc |  ios::out);
        
        if (!out_file)
        {
            cerr << "ERROR: Output file '" << output_filename << "'' could not be saved." << endl;
            cerr << getTidyResultsHeader();
            cout << getTidyResultsString(status);
        }
        else {
            out_file << getTidyResultsHeader();
            out_file << getTidyResultsString(status);
        }
    }
    else 
    {
        cerr << getTidyResultsHeader();
        cout << getTidyResultsString(status);
    }
}
// ----------------------------------------------------------------------------
void Typer::determineLowConfidence()
{
    for(auto&& item : best_allelic_profile)
    {
        double cov  = item.second.getCoverage();
        if(cov > 0.0 && cov < 100.0 )
        {
            item.second.setAlleleNumber(item.second.getAlleleNumber() + string(LOW_CONFIDENCE_SYMBOL));
        }
    }
}
// ----------------------------------------------------------------------------
string Typer::getStatusString(const Typer::StStatus & status)
{
    if (status == Typer::ST_PREDICTED)
        return "st_predicted";
    else if(status == Typer::NO_KMER_HITS)
        return "no_kmer_hits";
        
    return "no_st_in_table";
}
// ----------------------------------------------------------------------------
void Typer::determineSensitiveMode()
{
    // sensitive_mode_on = (options.getPrintAlleleCov() ||
    //                      options.getPrintKmerDepth() ||
    //                      options.getKmerDepthOutFile() != "");
    sensitive_mode_on = options.getSensitiveMode();
    fast_mode_on      = !sensitive_mode_on;
}
// ----------------------------------------------------------------------------
bool Typer::isCoverageCalculated()
{
    return (options.getPrintAlleleCov() || 
            options.getPrintKmerDepth() || 
            options.getKmerDepthOutFile() != "");
            // options.getKmerDepthOutFile() != "" || 
            // options.getExtKmerDepthOutFile() != "");
}
// ----------------------------------------------------------------------------
void Typer::printTopAllelesInfo()
{
    cerr << "top_alleles_info:" << endl;
    for (auto&& locus : loci_table)
    {
        TSeqInfoVector alleles = top_alleles_info.at(locus.id);
        cerr << "Locus: " << locus.id << endl;;
        printVector(alleles, "\n", std::cerr);
    }
}
// ----------------------------------------------------------------------------
void Typer::printBestAllelicProfile()
{
    cerr << "best_allelic_profile:" << endl;
    for (auto&& locus : loci_table)
    {
        cerr << "Locus: " << locus.id << endl;
        cerr << best_allelic_profile.at(locus.id) << endl;
    }
}
// ----------------------------------------------------------------------------
int Typer::run(int argc, char const ** argv)
{
    Timer timer;
    Timer general_timer;
    general_timer.start();
    stringstream time_message;
    stringstream command_line;
    string program_description = "STing \\fBtyper\\fP is an ultrafast "
                                 "assembly- and alignment-free program for sequence typing directly from "
                                 "whole-genome raw sequence reads. STing \\fBtyper\\fP is based on k-mer "
                                 "frequencies, and works with locus-based typing schemes like those defined "
                                 "in the traditional MLST method and its derivatives (e.g, rMLST, cgMLST or "
                                 "wgMLST). STing \\fBtyper\\fP requires an index (DB) created with the "
                                 "STing \\fBindexer\\fP program.";

    string program_name = TYPER_APP_NAME;

    setupArgumentParser(program_name, program_description);

    // Parse command line arguments and options
    ArgumentParser::ParseResult res = parseCommandLine(argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (options.getInFastq1Files().size() > 0 ||
            options.getInFastq2Files().size() > 0)
    {
        if ((options.getInFastq1Files().size() > 0 &&
                checkFiles(options.getInFastq1Files()) == FILE_NOT_FOUND) ||
                (options.getInFastq2Files().size() > 0 &&
                 checkFiles(options.getInFastq2Files()) == FILE_NOT_FOUND))
            return 1;

        determineSensitiveMode();
        setUpInFastqFiles();
        
        for(int i = 0; i < argc; ++i)
            command_line << argv[i] << ' ';

        command_line << endl;
        string message = "Command called: " + command_line.str();
        printMessage(message);
        
        timer.start();
        // Load alleles ESA index
        printMessage("Loading the input index...  ");
        loadAllelesEsaIndex();
        printMessage("Done!\n");
        timer.end();
        // time_message.str("");
        // time_message << "loadAllelesEsaIndex\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        
        time_message << endl << "============================================================" << endl;
        time_message << "Time table:" << endl;
        time_message << endl << "Process\tTime(s)" << endl;
        time_message << "LoadAllelesEsaIndex\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        timer.start(); 
        // Load the loci table
        loadLociTable(options.getIndexFilename() + string(LOCI_TABLE_EXT));
        timer.end();
        // time_message.str("");
        time_message << "loadLociTable\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        
        timer.start();
        // Process reads files: search k-mers in the loaded allele sequences index
        printMessage("Processing reads...  ");
        // Initialize the allele counts
        initAlleleCounts();
        // Initialize offset of top allele indices for the reduced index
        timer.end();
        // time_message.str("");
        time_message << "initAlleleCounts\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        
        timer.start();
        // Process all of the input read files
        // calledLocalProcessRead = false;
        // cerr << "DEBUG: calledLocalProcessRead: " << boolalpha << calledLocalProcessRead << endl;
        bool save_kmer_starts = false;
        processReadsFiles(save_kmer_starts);
        timer.end();
        // time_message.str("");
        time_message << "processReadsFiles (1st pass)\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        printMessage("  Done!\n");
        // cerr << "DEBUG: calledLocalProcessRead: " << boolalpha << calledLocalProcessRead << endl;
        
        // cerr << "DEBUG\t" << "total_reads\t" << total_reads << endl;
        // cerr << "DEBUG\t" << "total_kmers\t" << total_kmers << endl;
        // cerr << "DEBUG\t" << "total_hits\t" << total_hits << endl;
        // for(auto&& fqFile : kmerized_reads_indices) {
        //     cerr << "DEBUG\t" << "fq file: " << fqFile.first << "  reads: " << fqFile.second.size() << endl;
        //     // printVector(fqFile.second, ",", cout);
        // }
        
        timer.start();
        // Normalize kmer hit counts by the length of each allele
        normalizeCounts();
        timer.end();
        // time_message.str("");
        time_message << "normalizeCounts\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        
        stringstream msg;
        msg << "Getting the top " << options.getNTopAlleles() << " alleles of each locus...";
        printMessage(msg.str());
        timer.start();
        // Get the top N Alleles for each locus
        getTopNormAlleles(options.getNTopAlleles()); // **
        timer.end();
        // time_message.str("");
        time_message << "getTopNormAlleles\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        printMessage("  Done!\n"); 
        // cerr << "DEBUG: top_alleles_indices: ";
        // // cout << "DEBUG: " << top_alleles_info.size() << endl;
        // printVector(top_alleles_indices, string(", "), std::cerr);
        // cerr << endl;
        
        if (sensitive_mode_on || isCoverageCalculated())
        {
            // cerr << "DEBUG\t" << "total_reads\t" << total_reads << endl;
            // cerr << "DEBUG\t" << "total_kmers\t" << total_kmers << endl;
            // cerr << "DEBUG\t" << "total_hits\t" << total_hits << endl;
            // cerr << "DEBUG: Allele counts based on kmer_starts vector:" << endl;
            // for(auto&& aidx : kmer_starts) {
            //     cerr << "DEBUG\t"  << aidx.first << "\t" << allele_seq_ids[aidx.first] << "\t" << aidx.second.size() << endl;
            // }
            // re-initialize count variables
            // total_kmers = 0;
            // // total_reads = 0;
            // total_hits  = 0;
            // timer.start();
            // initAlleleCounts();
            // timer.end();
            // cerr << "initAlleleCounts (2nd)\t" << timer.getElapsedTime().count() << endl;
            // timer.reset();
            
            printMessage("Calculating k-mer depth and allele coverage... ");
            // cerr << "DEBUG: Before calling setupReducedIndexFromTopAlleles()" << endl;
            timer.start();
            // Build a ESA index from the sequences of the top N allele 
            setupReducedIndexFromTopAlleles();
            // cerr << "DEBUG: After calling setupReducedIndexFromTopAlleles()" << endl;
            timer.end();
            // time_message.str("");
            time_message << "setupReducedIndexFromTopAlleles\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
            
            // // DEBUG
            // for(auto&& locus : loci_table) {
            //     cerr << "DEBUG: top_alleles.at(locus.id).size(): " << top_alleles.at(locus.id).size() << endl;
            // }
            // cerr << "DEBUG: total_reads: " << total_reads << endl;
            
            // cerr << "DEBUG: run - Before process reads at the 2nd pass" << endl;
            timer.start();
            save_kmer_starts = true;
            processReadsFiles(save_kmer_starts);
            timer.end();
            // time_message.str("");
            time_message << "processReadsFiles (2nd pass)\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
            
            timer.start();
            // Compute k-mer depth and allele coverage
            computeKmerDepthAndAlleleCoverageForTopAlleles();
            timer.end();
            // time_message.str("");
            time_message << "computeKmerDepthAndAlleleCoverageForTopAlleles\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
            printMessage("  Done!\n");
            
            // DEBUG
            // for(auto&& locus : loci_table) 
            // {
            //     string locus_id        = locus.id;
            //     TSeqInfoVector alleles = top_alleles_info.at(locus_id);
            //     cout << "\nlocus_id: " << locus_id << endl;
            //     for(auto&& allele : alleles) {
            //         cerr << allele << endl;
            //     }
            // }
        }
        
        if (fast_mode_on)
        {   
            timer.start();
            printMessage("Calling alleles based on k-mer frequencies... ");
            // Get the alleles with the maximum counts
            selectBestAlleles();
            
            // cerr << "DEBUG: antes" << endl;
            // printBestAllelicProfile();
            // printTopAllelesInfo();
            
            if (isCoverageCalculated()){
                updateBestAllelicProfileCovData();
            }
            // cerr << "DEBUG: despues" << endl;
            // printBestAllelicProfile();
            // printTopAllelesInfo();
            
            timer.end();
            // time_message.str("");
            time_message << "selectBestAlleles\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
            printMessage("  Done!\n");
        }
        else if (sensitive_mode_on)
        {   
            timer.start();
            printMessage("Calling alleles based on k-mer frequencies and coverage info... ");
            // Select best alleles
            selectBestAllelesBasedOnCoverage();
            timer.end();
            // time_message.str("");
            time_message << "selectBestAllelesBasedOnCoverage\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
            printMessage("  Done!\n");
        }
        
        // if (options.getExtKmerDepthOutFile() != "") {
        if (options.getKmerDepthOutFile() != "") {
            printMessage(string("Saving per base k-mer depth data to '" + options.getKmerDepthOutFile() +  "'... "));
            timer.start();
            // saveTopAllelesPerBaseKmerDepthData(options.getExtKmerDepthOutFile());
            saveTopAllelesPerBaseKmerDepthData(options.getKmerDepthOutFile());
            timer.end();
            // time_message.str("");
            time_message << "saveTopAllelesPerBaseKmerDepthData\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
            printMessage("  Done!\n");
        }
        
        timer.start();
        // load the ESA index of the profiles table
        loadProfilesEsaIndex();
        timer.end();
        // time_message.str("");
        time_message << "loadProfilesEsaIndex\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        
        timer.start();
        // Load the file with the allelic profiles table
        loadProfilesFile();
        timer.end();
        // time_message.str("");
        time_message << "loadProfilesFile\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        
        printMessage("Predicting ST... ");
        string profile_string = "";
        profile_string        = createProfileStringRep(best_max_allele_numbers);
        
        timer.start();
        status = getInferredSt(inferred_st, profile_string);
        // cerr << "DEBUG: Typer::run - status: " << status << endl;
        timer.end();
        // time_message.str("");
        time_message << "getInferredSt\t" << timer.getElapsedTime().count() << endl;
        // printTime(time_message.str());
        timer.reset();
        printMessage("  Done!\n");
        
        if (sensitive_mode_on)
        {   
            timer.start();
            // Get low confidence marks on alleles depending on coverage
            determineLowConfidence();
            timer.end();
            // time_message.str("");
            time_message << "determineLowConfidence\t" << timer.getElapsedTime().count() << endl;
            // printTime(time_message.str());
            timer.reset();
        }
        
        general_timer.end();
        // time_message.str("");
        time_message << "Total time\t" << general_timer.getElapsedTime().count() << endl;
        time_message << "============================================================" << endl << endl;
        general_timer.reset();
        printTime(time_message.str());
        
        printMessage("Inferred allelic profile and its associated ST: \n\n");
        
        if (options.getTidyResults())
            printTidyResults(status);
        else
            printResults(status);

    }
    return 0;
}
// ----------------------------------------------------------------------------

// ============================================================================
// Main function
// ============================================================================
int main(int argc, char const ** argv)
{
    Typer typer_app;
    return typer_app.run(argc, argv);
    
    // // Test goToLine method
    // igzstream in_fastq_gz_file("test.txt.gz");
    // fstream in_fastq_file("test.txt");
    // string line;
    // vector<uint64_t> lines = {4, 6, 8};
    
    // cout << "plain text file:" << endl;
    // uint64_t start = 0;
    // for(unsigned i = 0; i < lines.size(); ++i) {
    //     if (i > 0) start = lines.at(i-1);
    //     uint64_t n_lines = lines.at(i) - start;
    //     advanceNLines(in_fastq_file, n_lines);
    //     getline(in_fastq_file, line);
    //     cout << "line: " << line << endl;
    // }
    
    // cout << "gzip'ed file:" << endl;
    // start = 0;
    // for(unsigned i = 0; i < lines.size(); ++i) {
    //     if (i > 0) start = lines.at(i-1);
    //     uint64_t n_lines = lines.at(i) - start;
    //     advanceNLines(in_fastq_gz_file, n_lines);
    //     getline(in_fastq_gz_file, line);
    //     cout << "line: " << line << endl;
    // }
    
}

// -----------------------------------------------------------------------------------
// Getters and setters
// -----------------------------------------------------------------------------------

const TCharEsaIndex& Typer::getProfilesIndex() const
{
    return profiles_index;
}

void Typer::setProfilesIndex(const TCharEsaIndex& profiles_index)
{
    this->profiles_index = profiles_index;
}

const ProfilesTable& Typer::getProfilesTable() const
{
    return profiles_table;
}

void Typer::setProfilesTable(const ProfilesTable & profiles_table)
{
    this->profiles_table = profiles_table;
}
