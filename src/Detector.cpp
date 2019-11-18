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


#include "Detector.h"

using namespace seqan;
using namespace std;

// -----------------------------------------------------------------------------------
// Constructors /destructors
// -----------------------------------------------------------------------------------

Detector::Detector()
{
    alleles_esa_index_exists = false;
}

Detector::Detector(const KmerMatcherOptions & options)
{
    alleles_esa_index_exists = false;
    this->options            = options;
}

Detector::~Detector()
{

}
// -----------------------------------------------------------------------------------
// Class methods (different from accessors)
// -----------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------
// Getters and setters
// -----------------------------------------------------------------------------------

/**
 * @brief      Setup configuration for argument parser 
 *
 * @param      program_name         The program name
 * @param      program_description  The program description
 */
void Detector::setupArgumentParser(string const & program_name,
                                   string const & program_description)
{
    string in_idx_param_name         = "INDEX_PREFIX_FILENAME";
    string in_fastq_1_param_name     = "FASTQ1";
    string in_fastq_2_param_name     = "FASTQ2";
    string sample_name_name          = "SAMPLE_NAME";
    string kmer_length_param_name    = "KMER_LENGTH";
    string output_file_param_name    = "OUTPUT_FILENAME";
    string threshold_param_name      = "THRESHOLD";
    string kmer_depth_out_param_name = "KMER_DEPTH_FILENAME";
    string verbose_param_name        = "VERBOSITY_LEVEL";
    
    // Setup ArgumentParser.
    setAppName(parser, program_name);
    setShortDescription(parser, "");
    setCategory(parser, "Gene Detection");
    setVersion(parser, DETECTOR_APP_VERSION);
    setDate(parser, DETECTOR_APP_UPDATE);
    
    // Usage message
    addUsageLine(parser, "-x <\\fI" +  in_idx_param_name + "\\fP>"
                         " -1 <\\fI" +  in_fastq_1_param_name + "\\fP> [\\fIOPTIONS\\fP]");
    
    // Description
    addDescription(parser, program_description);
    // addDescription(parser,
    //     "This program executes an MLST analysis from a set of "
    //     "short reads and an MLST database.");
    
    // Options
    addSection(parser, "Required input parameters");
    addOption(parser, ArgParseOption(
        "x", "index-prefix", "Database prefix filename.",
        ArgParseArgument::INPUT_FILE, in_idx_param_name));
    setRequired(parser, "index-prefix");
    
    addOption(parser, ArgParseOption(
        "1", "fastq-1-files", "Files with #1 mates, paired with "
        "files in \\fI" + in_fastq_1_param_name + "\\fP.",
        ArgParseArgument::INPUT_FILE, in_fastq_1_param_name, true));
    setRequired(parser, "fastq-1-files");
    
    // Input options
    addSection(parser, "Input options");
    addOption(parser, ArgParseOption(
        "2", "fastq-2-files", "Files with #2 mates, paired with "
        "files in \\fI" + in_fastq_2_param_name + "\\fP.",
        ArgParseArgument::INPUT_FILE, in_fastq_2_param_name, true));
    
    addOption(parser, ArgParseOption(
                  "s", "sample-name", "Name of the sample to be analized.",
                  ArgParseArgument::STRING, sample_name_name));
    
    addOption(parser, ArgParseOption(
        "k", "kmer-length", "Length of the k-mers to process the input reads.",
        ArgParseArgument::INTEGER, kmer_length_param_name));
    setDefaultValue(parser, "kmer-length", 30);
    
    addOption(parser, seqan::ArgParseOption(
        "r", "threshold", "Minimum length coverage (%) required to consider a gene as present in a sample.",
        ArgParseArgument::DOUBLE, threshold_param_name));
    setDefaultValue(parser, "threshold", 75.0);
    setMinValue(parser, "threshold", "1.0");
    setMaxValue(parser, "threshold", "100.0");
    // addOption(parser, seqan::ArgParseOption(
    //     "l", "long-results", "Select to print results in a long format (vertically). Default\twide format (horizontally)."));
    // Output options
    addSection(parser, "Output options");
    addOption(parser, seqan::ArgParseOption(
        "c", "kmer-counts", "Select to print the number of k-mer matches at each gene."));
        
    addOption(parser, seqan::ArgParseOption(
        "p", "kmer-perc", "Select to print the percentage of k-mer matches from the total of processed k-mers."));
    
    // addOption(parser, seqan::ArgParseOption(
    //     "r", "threshold", "Minimum percentage of k-mer matches required to consider a gene as present in a sample.",
    //     ArgParseArgument::DOUBLE, threshold_param_name));
    // setDefaultValue(parser, "threshold", 30.0);
    // setMinValue(parser, "threshold", "1.0");
    // setMaxValue(parser, "threshold", "100.0");
    
    addOption(parser, ArgParseOption(
        "g", "gene-cov", "Select to print the percent of the gene length that is covered by the corresponding k-mer matches."));
    
    addOption(parser, ArgParseOption(
        "d", "kmer-depth", "Select to print the mean k-mer depth of each gene."));
    
    addOption(parser, ArgParseOption(
        "o", "output-file", "Output filename to save the typing results (instead of stdout).",
        ArgParseArgument::OUTPUT_FILE, output_file_param_name));

    addOption(parser, ArgParseOption(
        "y", "print-tidy", "Select to print results in a tidy format."));
    
    addOption(parser, ArgParseOption(
        "t", "k-depth-file", "Output filename to save the detailed per-base k-mer depth data.",
        ArgParseArgument::OUTPUT_FILE, kmer_depth_out_param_name));
    
    addOption(parser, seqan::ArgParseOption(
        "v", "verbose", "Select to print informative messages (to stderr). Default\tnon verbose."));
    
    // Additional information
    addText(parser, "");
    addText(parser, "\n\\fI" + in_fastq_1_param_name + "\\fP and "
        "\\fI" + in_fastq_2_param_name + "\\fP can be comma-separated lists "
        "(no whitespace) and can be specified many times. "
        "E.g. -1 file1.fq,file2.fq -1 file3.fq");
}
// ----------------------------------------------------------------------------

/**
 * @brief      Parses the command line arguments and options
 *
 * @param[in]  argc  The number of arguments from std in
 * @param      argv  The arguments from std in
 *
 * @return     An ArgumentParser::ParseResult result of the parsing
 */
ArgumentParser::ParseResult 
Detector::parseCommandLine(int argc, char const ** argv)
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
    
    double threshold_perc;
    getOptionValue(threshold_perc, parser, "threshold");
    options.setThresholdPerc(threshold_perc);
    
    options.setPrintKmerCounts(isSet(parser, "kmer-counts"));
    
    options.setPrintKmerPerc(isSet(parser, "kmer-perc"));
    
    // options.setVerticalResults(isSet(parser, "long-results"));
    
    options.setPrintAlleleCov(isSet(parser, "gene-cov"));
    
    options.setPrintKmerDepth(isSet(parser, "kmer-depth"));
    
    string kmer_depth_out_filename = "";
    if (isSet(parser, "k-depth-file")) 
        getOptionValue(kmer_depth_out_filename, parser, "k-depth-file");
    
    options.setKmerDepthOutFile(kmer_depth_out_filename);
    
    options.setVerbose(isSet(parser, "verbose"));
    
    string output_filename = "";
    if (isSet(parser, "output-file")) 
        getOptionValue(output_filename, parser, "output-file");
    
    options.setOutputFile(output_filename);
    
    
    options.setTidyResults(isSet(parser, "print-tidy"));
    
    // uint64_t verbosity_level;
    // getOptionValue(verbosity_level, parser, "verbosity-level");
    // options.setVerbosityLevel(verbosity_level);
    
    return ArgumentParser::PARSE_OK;
}
// // ----------------------------------------------------------------------------
// void Detector::processReadsFiles(bool save_kmer_starts)
// {
//     if(length(fastq_files) > 0)
//     {
//         uint64_t n_fastq_files = length(fastq_files);
//         TSequenceEsaIndexFinder finder_esa(alleles_esa_index);

//         for (uint64_t i = 0; i < n_fastq_files; ++i)
//         {
//             TStringVector filename_tokens;
//             split(filename_tokens, getFileName(toCString(fastq_files[i])), ".");
//             string extension = filename_tokens[filename_tokens.size() - 1];

//             if (extension == "gz")  // process a gzipped fastq file
//             {
//                 igzstream in_fastq_gz_file(toCString(fastq_files[i]));
//                 if ( ! in_fastq_gz_file.good()) {
//                     cerr << "ERROR: Opening file '" << toCString(fastq_files[i]) << "' failed.\n";
//                     exit(1);
//                 }
//                 // processReadsByLine(in_fastq_gz_file, options.getKmerLength(), finder_esa);
//                 // if (save_kmer_starts)
//                 // {
//                 //     reprocessKmerizedReads(in_fastq_gz_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
//                 // }
//                 // else
//                 // {
//                     processReadsByLine(in_fastq_gz_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
//                 // }
//                 in_fastq_gz_file.close();
//             }
//             else if (extension == "fastq" || extension == "fq")     // process a uncompressed fastq file
//             {
//                 fstream in_fastq_file(toCString(fastq_files[i]), ios::in);
//                 if ( ! in_fastq_file.good()) {
//                     cerr << "ERROR: Opening file '" << toCString(fastq_files[i]) << "' failed.\n";
//                     exit(1);
//                 }
//                 // processReadsByLine(in_fastq_file, options.getKmerLength(), finder_esa);
//                 // if (save_kmer_starts)
//                 // {
//                 //     reprocessKmerizedReads(in_fastq_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
//                 // }
//                 // else
//                 // {
//                     processReadsByLine(in_fastq_file, i, options.getKmerLength(), finder_esa, save_kmer_starts);
//                 // }
//                 in_fastq_file.close();
//             }
//             else {
//                 cerr << "ERROR: The input file '" << toCString(fastq_files[i]) << "' has an invalid extension." << endl;
//                 exit(1);
//             }
//         }
//     }
// }
// ----------------------------------------------------------------------------

/**
 * @brief      Gets the indices of sequences with k-mer counts.
 */
void Detector::getSeqsWithCounts() 
{
    seqsWithCountsIndices = {};
    
    for(unsigned i = 0; i < allele_counts.size(); ++i)
    {
        if (allele_counts.at(i) > 0)
            seqsWithCountsIndices.push_back(i);
    }
}
// ----------------------------------------------------------------------------


/**
 * @brief      Calculates the percentage of k-mer matches from the total processed, for each gene.
 */
// vector<double> Detector::calculateMatchesPercentage() {
void Detector::calculateMatchesPercentage() {
    // vector<double> percentages;
    // vector<string> percentages;
    for(auto&& count : allele_counts) {
        // float percentage = pround(float(count)/float(total_kmers), 2);
        double percentage = 0.0;
        
        if (total_kmers > 0) 
        {
            // percentage = ((double)count/(double)total_kmers) * 100.0;
            percentage = ((double)count/(double)total_hits) * 100.0;
        }
        
        // percentages.push_back(percentage);
        kmer_matches_perc.push_back(percentage);
    }
    // return percentages;
}
// ----------------------------------------------------------------------------

/**
 * @brief      Calculates the presence/absence of each gene with k-mer counts based on its coverage.
 */
// vector<bool> Detector::calculatePresence(vector<double> const percentages) {
void Detector::calculatePresence() {
    // vector<bool> presence;
    // for (auto&& perc : kmer_matches_perc) 
    // { 
    //     presence.push_back(perc >= options.getThresholdPerc());
    // }
    
    // vector<bool> pres;
    for(auto&& seq_cov : coverages.at(0)) {
        presence.push_back(seq_cov.coverage >= options.getThresholdPerc());
        // pres.push_back(seq_cov.coverage >= options.getThresholdPerc());
    }
    
    // cout << "coverages.size()\t" << coverages.size() << endl;
    // cout << "coverages.at(0).size()\t" << coverages.at(0).size() << endl;
    // cout << "pres.size()\t" << pres.size() << endl;
    
    // return presence;
}
// ----------------------------------------------------------------------------

/**
 * @brief      Calculates the k-mer depth and allele coverage for each gene with counts.
 */
void Detector::computeKmerDepthAndAlleleCoverage()
{
    // vector<TUintVector> seq_kmer_depths;
    TUintVector seq_kmer_depths;
    TDoubleVector seq_mean_kmer_depths;
    TSeqCovVector seq_coverages;
    
    for(unsigned i = 0; i < seqsWithCountsIndices.size(); ++i)
    {
        uint64_t seq_idx = seqsWithCountsIndices.at(i);
        SeqCoverage seq_cov;
        // sort the kmer starts of the current seq
        std::sort(kmer_starts[seq_idx].begin(), kmer_starts[seq_idx].end());
        
        // Construct and store the intervals in a vector
        TIntervalVector intervals;
        for (unsigned j = 0; j < kmer_starts[seq_idx].size(); ++j) 
        {
            uint64_t start = kmer_starts[seq_idx].at(j);
            intervals.push_back(TInterval(start, start+options.getKmerLength(), seq_idx));
        }
        
        // Construct an interval tree to compute depth and coverage
        IntervalTree<uint64_t, uint64_t> tree(intervals);
        uint64_t seq_length     = length(indexText(alleles_esa_index)[seq_idx]);
        uint64_t depth_sum         = 0;
        uint64_t total_bases       = 0;
        bool check_for_depth_ratio = true;
        TUintVector depths;
        
        // Loop trhough each base of the current gene
        for (uint64_t base = 0; base < seq_length; ++base)
        {
            // Compute depth for the current base
            String<uint64_t> results;
            findIntervals(results, tree, base, base+1);
            uint64_t depth = length(results);
            
            if (depth == 0) check_for_depth_ratio = false;
            
            // Look for depth gaps
            if (check_for_depth_ratio && 
                (base >= BASES_MARGIN - 1) && (base <= seq_length - BASES_MARGIN - 1)) {
                
                double depth_ratio = (double)depth/(double)depths.at(base-1);
                
                // Look for depth disruptions
                // is there a gap?
                if (depth_ratio > DEPTH_DIFF_MAX_PEAK || depth_ratio < DEPTH_DIFF_MIN_PEAK)
                {
                    depth                 = 0;      // brake the coverage
                    check_for_depth_ratio = false;
                }
            }
            
            depths.push_back(depth);
            depth_sum += depth;
            
            if (depth > 0)
                ++total_bases;   // to calculate coverage
        }
        
        // seq_kmer_depths.push_back(depths);
        double mean_depth = (double)depth_sum/(double)seq_length;
        seq_mean_kmer_depths.push_back(mean_depth);
        
        // Compute the coverage for the current gene
        double coverage = 100 * ((double)total_bases / (double)seq_length);
        
        seq_cov.cov_bases = total_bases;
        seq_cov.seq_size  = seq_length;
        seq_cov.coverage  = coverage;
            
        seq_coverages.push_back(seq_cov);
        
        // Save the current gene per-base k-mer depth data to a file 
        if (options.getKmerDepthOutFile() != "") 
        {
            string out_filename = options.getKmerDepthOutFile();
            
            bool append_depth_data = true;
            if (i == 0) // first gene in 
                append_depth_data = false;
            
            if (saveKmerDepthData(options.getKmerDepthOutFile(), 
                                  // seq_kmer_depths,
                                  depths,
                                  seq_idx, 
                                  append_depth_data) == SAVE_ERROR) 
            {
                exit(1);
            } 
            else 
            { 
                if (i == seqsWithCountsIndices.size() - 1) {
                    string message = "Per-base k-mer depth data saved to '" + out_filename + "'.";
                    printMessage(message);
                }
            }
        }
        
    } // End genes loop
        
    // perbase_kmer_depths.push_back(seq_kmer_depths);
    mean_kmer_depths.push_back(seq_mean_kmer_depths);
    coverages.push_back(seq_coverages);
}
// ----------------------------------------------------------------------------
void Detector::computeKmerDepthAndAlleleCoverageNew()
{
    // vector<TUintVector> seq_kmer_depths;
    TUintVector seq_kmer_depths;
    TDoubleVector seq_mean_kmer_depths;
    TSeqCovVector seq_coverages;
    
    for(unsigned i = 0; i < seqsWithCountsIndices.size(); ++i)
    {
        uint64_t seq_idx = seqsWithCountsIndices.at(i);
        SeqCoverage seq_cov;
        // sort the kmer starts of the current seq
        std::sort(kmer_starts[seq_idx].begin(), kmer_starts[seq_idx].end());
        
        uint64_t seq_length        = length(indexText(alleles_esa_index)[seq_idx]);
        uint64_t depth_sum         = 0;
        uint64_t total_bases       = 0;
        bool check_for_depth_ratio = true;
        TUintVector depths(seq_length, 0);
        
        // Loop through k-mer starts and increment depth
        for (unsigned j = 0; j < kmer_starts[seq_idx].size(); ++j) 
        {
            uint64_t start = kmer_starts[seq_idx].at(j);
            
            for(unsigned k = start; k < (start + options.getKmerLength()); ++k)
                ++depths.at(k);
        }
        
        // Loop through each base of the current gene
        for (uint64_t base = 0; base < seq_length; ++base)
        {
            if (depths.at(base) == 0) check_for_depth_ratio = false;
            
            // Look for depth gaps
            if (check_for_depth_ratio && 
                (base >= BASES_MARGIN - 1) && (base <= seq_length - BASES_MARGIN - 1)) {
                
                double depth_ratio = (double)depths.at(base)/(double)depths.at(base-1);
                
                // Look for depth disruptions
                // is there a gap?
                if (depth_ratio > DEPTH_DIFF_MAX_PEAK || depth_ratio < DEPTH_DIFF_MIN_PEAK)
                {
                    depths.at(base)       = 0;      // brake the coverage
                    check_for_depth_ratio = false;
                }
            }
            
            // depths.push_back(depth);
            depth_sum += depths.at(base);
            
            if (depths.at(base) > 0)
                ++total_bases;   // to calculate coverage
        }
        
        // seq_kmer_depths.push_back(depths);
        double mean_depth = (double)depth_sum / (double)seq_length;
        seq_mean_kmer_depths.push_back(mean_depth);
        
        // Compute the coverage for the current gene
        double coverage = 100 * ((double)total_bases / (double)seq_length);
        
        seq_cov.cov_bases = total_bases;
        seq_cov.seq_size  = seq_length;
        seq_cov.coverage  = coverage;
            
        seq_coverages.push_back(seq_cov);
        
        // Save the current gene per-base k-mer depth data to a file 
        if (options.getKmerDepthOutFile() != "") 
        {
            string out_filename = options.getKmerDepthOutFile();
            
            // bool append_depth_data = true;
            // if (i == 0) // first gene in 
            //     append_depth_data = false;
            bool append_depth_data = (i > 0);
            
            if (saveKmerDepthData(options.getKmerDepthOutFile(), 
                                  // seq_kmer_depths,
                                  depths,
                                  seq_idx, 
                                  append_depth_data) == SAVE_ERROR) 
            {
                exit(1);
            } 
            else 
            { 
                if (i == seqsWithCountsIndices.size() - 1) {
                    string message = "Per-base k-mer depth data saved to '" + out_filename + "'.";
                    printMessage(message);
                }
            }
        }
        
    } // End genes loop
        
    // perbase_kmer_depths.push_back(seq_kmer_depths);
    mean_kmer_depths.push_back(seq_mean_kmer_depths);
    coverages.push_back(seq_coverages);
}
// ----------------------------------------------------------------------------
template<typename T>
void Detector::printResultsLine(const string &    line_type,
                                const vector<T> & info_vector)
{
    if (options.getSampleName() != "")
        cout << options.getSampleName() << FIELD_SEPARATOR;
    
    cout << line_type << FIELD_SEPARATOR;
    printVector(info_vector, FIELD_SEPARATOR, std::cout);
    cout << FIELD_SEPARATOR 
         << total_hits
         << FIELD_SEPARATOR
         << total_kmers
         << FIELD_SEPARATOR
         << total_reads
         << FIELD_SEPARATOR
         << fastqFilesToString()
         << endl;
}

// ----------------------------------------------------------------------------
template<typename T>
string Detector::getResultsLine(const string &    line_type,
                                const vector<T> & info_vector)
{
    stringstream out_string;
    if (options.getSampleName() != "")
        out_string << options.getSampleName() << FIELD_SEPARATOR;
    
    out_string << line_type << FIELD_SEPARATOR;
    
    if (info_vector.size() > 0)
        printVector(info_vector, FIELD_SEPARATOR, out_string);
    else 
        out_string << NA_STRING;
    
    out_string << FIELD_SEPARATOR 
               << total_hits
               << FIELD_SEPARATOR
               << total_kmers
               << FIELD_SEPARATOR
               << total_reads
               << FIELD_SEPARATOR
               << fastqFilesToString()
               << endl;
               
    return out_string.str();
}
// ----------------------------------------------------------------------------

/**
 * @brief      Prints results in an horizontal way.
 */
void Detector::printResultsHorizontal()
{
    uint64_t header_length;
    string output_filename = options.getOutputFile();
    
    if (output_filename != "") {
        fstream out_file;
        out_file.open(output_filename, ios::trunc |  ios::out);
        
        if (!out_file)
        {
            cerr << "ERROR: Output file '" << output_filename << "'' could not be saved." << endl;
            cerr << getResultsHeader(header_length);
            cout << getResultsString();
        }
        else {
            out_file << getResultsHeader(header_length);
            out_file << getResultsString();
        }
    }
    else 
    {
        cerr << getResultsHeader(header_length);
        cout << getResultsString();
    }
    
}
// ----------------------------------------------------------------------------
string Detector::getResultsHeader(uint64_t & header_length)
{
    TStringVector header;
    stringstream out_string;
    // Prepare header
    if (options.getSampleName() != "")
        header.push_back("Sample");
    header.push_back("Line_type");
    
    if (seqsWithCountsIndices.size() > 0)
    {
        // Get the seq ids from genes with a k-mer count greater than 0
        for(auto&& i : seqsWithCountsIndices)
        {
            header.push_back(toCString(allele_seq_ids[i]));
        }
    }
    else {
        header.push_back("Gene");
    }
    
    // finishing the header line
    header.push_back("Total_hits");
    header.push_back("Total_kmers");
    header.push_back("Total_reads");
    header.push_back("Input_files");
    
    header_length = header.size();
    
    // print header
    printVector(header, FIELD_SEPARATOR, out_string);
    out_string << endl;
    
    return out_string.str();
}
// ----------------------------------------------------------------------------
string Detector::getResultsString()
{
    stringstream out_string;
    TUintVector counts;
    TStringVector percs;
    
    // Get the numbers from genes with a k-mer count greater than 0
    for(auto&& i : seqsWithCountsIndices)
    {
        counts.push_back(allele_counts.at(i));
        percs.push_back(getNumberStringWithPrecision(kmer_matches_perc.at(i),1));
    }
        
    out_string << getResultsLine("presence", presence);
    if (options.getPrintKmerCounts())
        out_string << getResultsLine("kmer_counts", counts);
    if(options.getPrintKmerPerc())
        out_string << getResultsLine("kmer_perc", percs);
    if (options.getPrintAlleleCov())
        out_string << getResultsLine("allele_coverage", getAlleleCoverageForPrinting(0));
    if (options.getPrintKmerDepth())
        out_string << getResultsLine("mean_kmer_depth", getMeanKmerDepthForPrinting(0));
    
    return out_string.str();
}

// ----------------------------------------------------------------------------
string Detector::getTidyResultsHeader(uint64_t & header_length)
{
    TStringVector header;
    stringstream out_string;
    // Prepare header
    if (options.getSampleName() != "")
        header.push_back("Sample");
    header.push_back("Gene");
    header.push_back("Presence");
    if (options.getPrintKmerCounts())
        header.push_back("k-mer_counts");
    if(options.getPrintKmerPerc())
        header.push_back("k-mer_perc");
    if (options.getPrintAlleleCov())
    {
        header.push_back("Gene_cov(%)");
        header.push_back("Covered_bases");
        header.push_back("Gene_length");
    }
    if (options.getPrintKmerDepth())
        header.push_back("Mean_k-mer_depth");
    header.push_back("Total_hits");
    header.push_back("Total_k-mers");
    header.push_back("Total_reads");
    header.push_back("Input_files");
    
    header_length = header.size();
    
    // print header
    printVector(header, FIELD_SEPARATOR, out_string);
    out_string << endl;
    
    return out_string.str();
}

// ----------------------------------------------------------------------------
string Detector::getTidyResultsString(const uint64_t & header_length)
{
    stringstream out_string;
    // string out_string = "";
    if (seqsWithCountsIndices.size() > 0)
    {
        for(unsigned i = 0; i < seqsWithCountsIndices.size(); ++i)
        {
            unsigned seq_idx       = seqsWithCountsIndices.at(i);
            SeqCoverage seq_cov    = coverages.at(0).at(i);
            double mean_kmer_depth = mean_kmer_depths.at(0).at(i);
            
            if (options.getSampleName() != "")
                out_string << options.getSampleName() << FIELD_SEPARATOR;
            out_string << toCString(allele_seq_ids[seq_idx])
                       << FIELD_SEPARATOR
                       << presence.at(i)
                       << FIELD_SEPARATOR;
            if (options.getPrintKmerCounts())
            {
                out_string << allele_counts.at(seq_idx)
                     << FIELD_SEPARATOR;
            }
            if(options.getPrintKmerPerc())
            {
                out_string << getNumberStringWithPrecision(kmer_matches_perc.at(seq_idx), 1)
                     << FIELD_SEPARATOR;
            }
            if (options.getPrintAlleleCov())
            {
                out_string << getNumberStringWithPrecision(seq_cov.coverage, 1)
                           << FIELD_SEPARATOR
                           << seq_cov.cov_bases
                           << FIELD_SEPARATOR
                           << seq_cov.seq_size
                           << FIELD_SEPARATOR;
            }
            if (options.getPrintKmerDepth())
            {
                out_string << getNumberStringWithPrecision(mean_kmer_depth, 1)
                     << FIELD_SEPARATOR;
            }
            // out_string << total_hits
            //            << FIELD_SEPARATOR
            //            << total_kmers
            //            << FIELD_SEPARATOR
            //            << total_reads
            //            << FIELD_SEPARATOR
            //            << fastqFilesToString()
            //            << endl;
        }
    }
    else 
    {
        uint64_t n_items = header_length - 4;
        if (options.getSampleName() != "")
        {
            out_string << options.getSampleName() << FIELD_SEPARATOR;
            --n_items;
        }
        string delim = "";
        for(unsigned i = 0; i < n_items; ++i)
        {
            out_string << NA_STRING << FIELD_SEPARATOR;
        }
        
        // out_string << fastqFilesToString()
        //            << endl;
    }
    
    out_string << total_hits
               << FIELD_SEPARATOR
               << total_kmers
               << FIELD_SEPARATOR
               << total_reads
               << FIELD_SEPARATOR
               << fastqFilesToString()
               << endl;
    
    return out_string.str();
}

// ----------------------------------------------------------------------------
/**
 * @brief      Prints results in an horizontal way.
 */
void Detector::printTidyResults()
{
    uint64_t header_length;
    string output_filename = options.getOutputFile();
    
    if (output_filename != "") {
        fstream out_file;
        out_file.open(output_filename, ios::trunc |  ios::out);
        
        if (!out_file)
        {
            cerr << "ERROR: Output file '" << output_filename << "'' could not be saved." << endl;
            cerr << getTidyResultsHeader(header_length);
            cout << getTidyResultsString(header_length);
        }
        else {
            out_file << getTidyResultsHeader(header_length);
            out_file << getTidyResultsString(header_length);
        }
    }
    else 
    {
        cerr << getTidyResultsHeader(header_length);
        cout << getTidyResultsString(header_length);
    }
}

// ----------------------------------------------------------------------------
/**
 * @brief      Runs the detector application. Contains the main logic of the detector and calls all of the rest of protected/private methods.
 *
 * @param[in]  argc  The number of argments from sdtin.
 * @param      argv  The arguments from stdin
 *
 * @return     An int that indicates success (0) or error (1) in the app execution.
 */
int Detector::run(int argc, char const ** argv)
{
    Timer timer;
    string program_description = "STing \\fBdetector\\fP is an ultrafast "
                                 "assembly- and alignment-free program for detecting genes directly from "
                                 "NGS raw sequence reads. STing \\fBdetector\\fP is based on k-mer "
                                 "frequencies. STing \\fBdetector\\fP requires an index (DB) created with the "
                                 "STing \\fBindexer\\fP program (using the GDETECT mode).";;
    
    string program_name = DETECTOR_APP_NAME;
                                  
    setupArgumentParser(program_name, program_description);

    // cerr << "Running class-based gene detector app!" << endl << endl;
    
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
        
        timer.start();
        setUpInFastqFiles();
        timer.end();
        // cerr << "setUpInFastqFiles\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        timer.start();
        // Load alleles ESA index
        printMessage("Loading the input index...");
        loadAllelesEsaIndex();
        printMessage("  Done!\n");
        // printMessage("  Index successfully loaded\n");
        timer.end();
        // cerr << "loadAllelesEsaIndex\t" << timer.getElapsedTime().count() << endl;
        timer.reset();

        timer.start();
        // Process reads files\tsearch k-mers in the loaded allele sequences index
        // cerr << "Processing reads..." << endl;
        printMessage("Processing reads...");
        // Initialize the allele counts
        initAlleleCounts();
        timer.end();
        // cerr << "initAlleleCounts\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        timer.start();
        bool save_kmer_starts = true;
        // Process all of the input read files
        processReadsFiles(save_kmer_starts);
        printMessage("  Done!\n");
        timer.end();
        // cerr << "processReadsFiles\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        // cout << "total_reads: " << total_reads << endl;
        // // Check that allele_counts have been updated after process all of the reads
        // TUintVector counts = getAlleleCounts();
        // for(uint64_t i = 0; i < counts.size(); ++i) {
        //     if (counts[i] > 0)
        //         cout << i << "\t"<< counts[i] << endl;
        // }

        timer.start();
        // Load the loci table
        loadLociTable(options.getIndexFilename() + string(LOCI_TABLE_EXT));
        timer.end();
        // cerr << "loadLociTable\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        
        timer.start();
        getSeqsWithCounts();
        timer.end();
        // cerr << "getSeqsWithCounts\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        timer.start();
        // Normalize counts by gene length
        normalizeCounts();
        timer.end();
        // cerr << "normalizeCounts\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        
        timer.start();
        // computeKmerDepthAndAlleleCoverage();
        computeKmerDepthAndAlleleCoverageNew();
        timer.end();
        // cerr << "computeKmerDepthAndAlleleCoverage\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        timer.start();
        calculateMatchesPercentage();
        timer.end();
        // cerr << "calculateMatchesPercentage\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        timer.start();
        calculatePresence();
        timer.end();
        // cerr << "calculatePresence\t" << timer.getElapsedTime().count() << endl;
        timer.reset();
        
        printMessage("\nResults:\n\n");
        
        // if (options.getVerticalResults())
        //     printResultsVertical();
        // else
        
        // timer.start();
        
        if(options.getTidyResults())
            printTidyResults();
        else
            printResultsHorizontal();
        
        // timer.end();
        // cerr << "Printing results\t" << timer.getElapsedTime().count() << endl;
        // timer.reset();
        
    }
    return 0;
}
// ----------------------------------------------------------------------------
// ============================================================================
// Main function
// ============================================================================

/**
 * @brief      Main function
 *
 * @param[in]  argc  The number of arguments from stdin
 * @param      argv  The arguments from stdin
 *
 * @return     An int that indicates success (0) or error (1) in the app execution.
 */
int main(int argc, char const ** argv)
{
    Detector gene_detector_app;
    return gene_detector_app.run(argc, argv);
}
