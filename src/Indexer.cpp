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


#include "Indexer.h"

// ============================================================================
// Namespaces
// ============================================================================
using namespace seqan;
using namespace std;

// -----------------------------------------------------------------------------------
// Constructors /destructors
// -----------------------------------------------------------------------------------

Indexer::Indexer()
{
    program_name = INDEXER_APP_NAME;
}

/**
 * @brief      Constructs the object.
 *
 * @param[in]  options  The program options
 */
Indexer::Indexer(const IndexerOptions & options)
{
    this->options = options;
}

// -----------------------------------------------------------------------------------
// Class methods (different from accessors)
// -----------------------------------------------------------------------------------


/**
 * @brief      Configures the ArgumentParser object of the Indexer class
 */
void Indexer::setupArgumentParser()
{
    string program_name           = INDEXER_APP_NAME;
    string config_file_param_name = "CONFIG_FILE";
    string out_prefix_param_name  = "PREFIX";
    string mode_param_name        = "MODE";

    // Setup ArgumentParser.
    setAppName(parser, program_name);
    setShortDescription(parser, "");
    setCategory(parser, "Bacterial Typing and Gene Detection");
    setVersion(parser, INDEXER_APP_VERSION);
    setDate(parser, INDEXER_APP_UPDATE);

    // Usage message
    addUsageLine(parser, "[\\fIOPTIONS\\fP] -c <\\fI" + config_file_param_name + "\\fP>");

    // Description
    addDescription(parser,
        "STing \\fBindexer\\fP creates indexes (DBs) required by the STing \\fBtyper\\fP "
        "and \\fBdetector\\fP programs for loci-based typing analysis and "
        "detecting genes, respectively, from NGS raw sequence reads.");

    // Options
    addOption(parser, ArgParseOption(
        "c", "config-file", "A tab delimited file whith names and paths to "
        "the typing scheme files (see the FILE FORMAT DETAILS section below).",
        ArgParseArgument::INPUT_FILE, config_file_param_name));
    setRequired(parser, "config-file");

    addOption(parser, ArgParseOption(
        "p", "db-prefix", "Filename prefix for the DB files to be created. "
        "You can specify a folder structure here to store your DB at a "
        "particular location, e.g., path/to/my/db/prefix. "
        "Default: name of the config file \\fI" + config_file_param_name + "\\fP",
        ArgParseArgument::OUTPUT_FILE , out_prefix_param_name));
    
    addOption(parser, ArgParseOption(
        "m", "mode", "Indexing mode. Valid options: MLST, GDETECT. "
        "Select MLST to create a database for MLST analysis or "
        "GDETECT to create a database for gene detection.",
        ArgParseArgument::STRING , mode_param_name));
    setDefaultValue(parser, "mode", "MLST");
    
    // Additional information
    addTextSection(parser, "File format details");

    addText(parser, "");
    addTextSubSection(parser, config_file_param_name);
    addText(parser, "A tab separated file with the name and location of files "
                    "for creating a DB.");
    addText(parser, "Config files for MLST DBs (MLST mode) must have two sections: "
                    "\\fB[loci]\\fP that describes names and paths to alleles " 
                    "sequence files for each locus, and \\fB[profile]\\fP that " 
                    "describes the name and path to the profile file.");
    addText(parser, "Config files for gene detection DBs (GDETECT mode), only "
                    "require the [loci] section.");
    addText(parser, "An example of a config file for a MLST DB is as follows:");
    addListItem(parser, "", "");
    addListItem(parser, "[loci]", "");
    addListItem(parser, "locus1   relative/path/to/locusFile1", "");
    addListItem(parser, "locus2   relative/path/to/locusFile2", "");
    addListItem(parser, "locusN   relative/path/to/locusFileN", "");
    addListItem(parser, "[profile]", "");
    addListItem(parser, "profile  relative/path/to/profileFile", "");
    addText(parser, "");
    addText(parser, "Paths are relative to the config file itself. Blank lines "
                    "and comments (lines starting with '#') will be ignored.");

    addText(parser, "");
    addTextSubSection(parser, "Allele sequence file");
    addText(parser, "A standard multi-FASTA file (.fa or .fasta) in which each "
                    "sequence description must be the locus name and the allele "
                    "number separated by '_':");
                    
    addListItem(parser, "", "");
    addListItem(parser, ">abcZ_1", "");
    addListItem(parser, "TTTGATACTGTTGCCGA...", "");
    addListItem(parser, ">abcZ_2", "");
    addListItem(parser, "TTTGATACTGTTGCCGA...", "");

    addTextSubSection(parser, "Profile file");
    addText(parser, "A tab separated file that contains the ST and the "
        "corresponding allelic profile:");
    addListItem(parser, "", "");
    addListItem(parser, "ST  abcZ  adk  aroE  fumC  gdh  pdhC  pgm  clonal_complex", "");
    addListItem(parser, "1   1     3    1     1     1    1     3    ST-1 complex/subgroup I/II", "");
    addListItem(parser, "2   1     3    4     7     1    1     3    ST-1 complex/subgroup I/II", "");
    addListItem(parser, "3   1     3    1     1     1    23    13   ST-1 complex/subgroup I/II", "");
    addListItem(parser, "4   1     3    3     1     4    2     3    ST-4 complex/subgroup IV", "");
    addText(parser, "");
}
// ----------------------------------------------------------------------------

/**
 * @brief      Parses the command line arguments of the Indexer tool
 *
 * @param[in]  argc  The argc
 * @param      argv  The argv
 *
 * @return     PARSE_OK if everything is okay, otherwise PARSE_ERROR
 */
ArgumentParser::ParseResult
Indexer::parseCommandLine(int argc, char const ** argv)
{
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    string config_filename;
    getOptionValue(config_filename, parser, "config-file");
    options.setConfigFileFilename(config_filename);

    string prefix_filename;
    getOptionValue(prefix_filename, parser, "db-prefix");
    options.setOutPrefixFilename(prefix_filename);

    string mode;
    getOptionValue(mode, parser, "mode");
    
    if (mode == MLST_MODE_OPTION)
        options.setMode(MLST_MODE);
    else if (mode == GENE_DETECTOR_MODE_OPTION)
        options.setMode(GENE_DETECTOR_MODE);
    else
    {
        printHelp(parser);
        cerr << "\nERROR: Invalid option for -m (--mode) argument. "
                "Valid options are '"<< MLST_MODE_OPTION<< "' and '" << GENE_DETECTOR_MODE_OPTION <<"'"<< endl;
        return ArgumentParser::PARSE_ERROR;
    }
    
    if (!isSet(parser, "db-prefix"))
        options.setOutPrefixFilename(getFileName(options.getConfigFileFilename()));

    return ArgumentParser::PARSE_OK;
}

// -----------------------------------------------------------------------------------

/**
 * @brief      Process the input config file
 */
void Indexer::processMlstConfigFile()
{
    // Process the config file
    fstream infile(options.getConfigFileFilename(), ios::in);
    
    string line;
    string sep        = FIELD_SEPARATOR;
    string config_key = "";

    while(getline(infile, line))
    {
        if (!line.empty() && line.at(0) != COMMENT_CHAR) // The current line is not an empty/comment line
        {
            if(line == LOCI_CONFIG_HEADER)          // loci section header
                config_key = LOCI_CONFIG_KEY;

            else if(line == PROFILE_CONFIG_HEADER)  // profile section header
                config_key = PROFILE_CONFIG_KEY;
            
            else                                    // loci/profile item line
            {
                TStringVector fields;
                split(fields, line, FIELD_SEPARATOR);

                if(fields.size() > 1)
                    config_map[config_key][fields[0]] = fields[1];

                else
                {
                    cerr << "ERROR: Config file is not properly formatted!. "
                            "Please check for tabs as separators instead of spaces." << endl;
                    exit(1);
                }
            }
        }
    }
    // check if there were the loci and profile sections in the config file
    
    if (config_map.find(LOCI_CONFIG_KEY) == config_map.end())
    {
        cerr << "ERROR: The config file does not contain the [loci] section." 
             << endl;
        exit(1);
    }
    if (config_map.find(PROFILE_CONFIG_KEY) == config_map.end())
    {
        cerr << "ERROR: The config file does not contain the [profile] section." 
             << endl;
        exit(1);
    }
        
        
}
// ----------------------------------------------------------------------------

/**
 * @brief      Process the input config file
 */
void Indexer::processGeneDetectorConfigFile()
{
    // Process the config file
    fstream infile(options.getConfigFileFilename(), ios::in);

    string line;
    string sep        = FIELD_SEPARATOR;
    string config_key = "";

    // cout << "Loading sequences from loci files:"<< endl << endl;
    // cout << "#" << sep << "Seqs." << sep << "File" << endl;

    while(getline(infile, line))
    {
        if (!line.empty() && line.at(0) != COMMENT_CHAR) // The current line is not an empty/comment line
        {
            if(line == LOCI_CONFIG_HEADER)          // loci section header
                config_key = LOCI_CONFIG_KEY;
            
            else if(line == PROFILE_CONFIG_HEADER)  // profile section header
            {
                cerr << "ERROR: The config file contains the [profile] section. "
                        "It seems like you are using an MLST config file." 
                     << endl;
                exit(1);
            }

            else                                    // loci/profile item line
            {
                TStringVector fields;
                split(fields, line, FIELD_SEPARATOR);

                if(fields.size() > 1)
                    config_map[config_key][fields[0]] = fields[1];

                else
                {
                    printVector(fields, "\t", cout);
                    cerr << "ERROR: Config file is not properly formatted!. "
                            "Please check for tabs as separators instead of spaces." << endl;
                    exit(1);
                }
            }
        }
    }
    // check if there was the loci in the config file
    
    if (config_map.find(LOCI_CONFIG_KEY) == config_map.end())
    {
        cerr << "ERROR: The config file does not contain the [loci] section." 
             << endl;
        exit(1);
    }
}
// ----------------------------------------------------------------------------
void Indexer::loadLociFiles()
{
    
    // Separar la construcción de la tabla de loci 
    // de la lectura de secuencias de cada loci
    // Así se asegura la construcción del índice de secuencias 
    // en el mismo orden de las columnas del perfil.
    // 
    //  1. Se construye la tabla de loci 
    //      --> void buildLociTable();
    //  2. Se carga el archivo de perfiles 
    //      --> loadProfilesFile(string const &  filename)
    //  3. Se cargan las secuencias de cada loci en el orden de las columnas de los perfiles
    //      --> loadLociFiles()

    uint64_t seqs_counter = 0;
    uint64_t total_seqs   = 0;
    uint64_t file_counter = 0;
    string sep            = FIELD_SEPARATOR;
    
    cout << "Loading sequences from sequences files:"<< endl << endl;
    cout << "N"       << sep 
         << "Loci"    << sep
         << "#Seqs." << sep
         << "File"    << endl;
    
    for (auto & locus_item : config_map[LOCI_CONFIG_KEY])
    {
        string locus_full_path;
        string config_file_dir;
        string config_file_base;
        relativeDirBaseSplit(config_file_dir, config_file_base, options.getConfigFileFilename());
        
        if (config_file_dir == "")
            config_file_dir = ".";
        
        locus_full_path = pathAppend(config_file_dir, locus_item.second);
        
        // SLocusRecord locus = {locus_item.first, locus_item.second, 0, 0};
        SLocusRecord locus = {locus_item.first, locus_full_path, 0, 0};
        
        // Check the current locus file
        if(checkFile(locus.filename) == FILE_NOT_FOUND) 
            exit(1);
        
        // Open ref_file (fasta file).
        SeqFileIn seq_file_in(toCString(locus.filename));
        
        // Read all sequences from the current file
        seqs_counter = length(seqs);
        readRecords(ids, seqs, seq_file_in);
        
        total_seqs         = length(seqs);
        locus.n_seqs       = total_seqs - seqs_counter;
        locus.last_seq_idx = length(seqs) - 1;
        
        loci_table.push_back(locus);
        ++file_counter;
        
        for (uint64_t i = 0; i < locus.n_seqs; ++i)
        {
            loci_alleles_index.push_back(file_counter);
        }
        
        cout << file_counter   << sep
             << locus.id       << sep 
             << locus.n_seqs   << sep 
             << locus.filename << endl;
    }
    
    toUpper(seqs);
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Loads the profiles file defined in the input config file.
 *
 * @param      filename  The filename of the profiles file
 */
// void Indexer::loadProfilesFile(const string & filename)
void Indexer::loadProfilesFile()
{
    
    TConfigMap tmp_config_map = getConfigMap();
    for (auto profile_item : tmp_config_map[PROFILE_CONFIG_KEY])
    {
        string profile_full_path;
        string config_file_path = getPathName(options.getConfigFileFilename());
        
        if (config_file_path.compare("") == 0)
            config_file_path = CURRENT_DIRECTORY;
            
        profile_full_path = pathAppend(config_file_path, profile_item.second);
        
        setProfileName(profile_item.first);
        
        unordered_map<string, vector<string>> table;
        fstream profiles_file(profile_full_path, ios::in);
        
        if (profiles_file.is_open()) 
        {
            string line;
            vector<string> col_names;
            uint64_t current_line = 0;
            while(getline(profiles_file, line))
            {
                vector<string> fields;
                
                if (!line.empty() && line.at(0) != COMMENT_CHAR)    // The current line is not an empty/comment line
                {
                    split(fields, line, FIELD_SEPARATOR);
                    
                    if (++current_line == 1)    // get header and create the map
                    {   
                        profiles_table.cols = fields;
                        profiles_table.cols.at(0) = "ST";
                        for(string col : profiles_table.cols) {
                            table[col] = vector<string>();
                            profiles_table.table[col] = table[col];
                        }
                    }
                    else {
                        
                        for(uint64_t i = 0; i < profiles_table.cols.size(); ++i) {
                            string col = profiles_table.cols[i];
                            profiles_table.table[col].push_back(fields[i]);
                        }
                    }
                }
            }
            
            profiles_file.close();
        } 
        else {
            cerr << "ERROR: Could not open the profiles file '" << profile_full_path << "'\n";
            exit(1);
        }
    }
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Gets the string representation of all the profiles defined in the
 *             profiles file, to be indexed using an ESA index
 *
 * @param[out] profiles  The StringSet of the profiles
 */
void Indexer::getProfileStrings(TCharStringSet & profiles)
{
    TStringVector loci_ids;
    // get the loci ids from the lociTable to preserve the same order in which
    // the loci sequences were loaded
    for(auto&& locus : loci_table) {
        loci_ids.push_back(locus.id);
    }
    
    // printVector(loci_ids, ", ", cout); cout << endl;
    
    // Check if the number of columns in the profiles file is at least the number of
    // loci plus one which are specified in the config fle
    if (profiles_table.cols.size() > loci_ids.size()) 
    {
        // loop through the loci table merging the allele numbers of each 
        // profile in a single string
        string sep = "-";
        for(uint64_t i = 0; i < profiles_table.table["ST"].size(); ++i) 
        {
            string current_prof = "*";
            for(auto&& clocus : loci_table)
            {
                string locus = clocus.id;
                if (profiles_table.table.find(locus) != profiles_table.table.end())
                {
                    current_prof = current_prof + profiles_table.table[locus][i] + sep;
                }
                else
                {
                    cerr << endl 
                         << "ERROR: The locus '" << locus << "' in the config "
                            "file, does not exist in the profiles file.\n" 
                            "Please check that each locus in the config file has a "
                            "corresponding column in the profiles file with the same name." << endl;
                    exit(1);
                }
                    
            }
            appendValue(profiles, current_prof);
        }
    }
    else 
    {
        cerr << endl 
             << "ERROR: At least " << loci_ids.size() + 1 
             << " columns (a ST column + # loci in config file) are required in the profiles file but only "
             << profiles_table.cols.size() << " were found."
             << endl;
        exit(1);
    }
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Creates an ESA index of the string representation of the allelic
 *             profiles.
 *
 * @param[in]  profiles          The StringSet of the string representation of
 *                               the allelic profiles
 * @param[in]  output_file_name  The output file name of the index
 *
 * @return     SAVE_SUCCES if the index was saved correctly, otherwise SAVE_ERROR.
 */
SaveFileResult Indexer::createProfilesIndex(const TCharStringSet & profiles, 
                                            const string &         output_file_name)
{
    return createIndex(profiles_esa_index, profiles, output_file_name);
}
// ----------------------------------------------------------------------------

/**
 * @brief      Creates an ESA index of the alleles of each locus defined in the config file.
 *
 * @param      output_file_name  The output file name of the index
 *
 * @return     SAVE_SUCCES if the index was saved correctly, otherwise SAVE_ERROR.
 */
SaveFileResult Indexer::createAllelesIndex(const string & output_file_name)
{
    return createIndex(alleles_esa_index, seqs, output_file_name);
}
// ----------------------------------------------------------------------------

/**
 * @brief      Saves a loci table which contains each locus name, its filename
 *             the number of alleles (sequences) and the index of its last
 *             allele.
 *
 * @param      output_file_name  The output file name
 *
 * @return     SAVE_SUCCES if the file was saved correctly, otherwise SAVE_ERROR.
 */
SaveFileResult Indexer::saveLociTable(const string & output_file_name)
{
    fstream out_file(output_file_name, ios::out);
    string sep = FIELD_SEPARATOR;
    
    if (!out_file)
    {
        cerr << "ERROR: The loci sequences data file could not be saved." << endl;
        return SAVE_ERROR;
    }
    
    for (uint64_t i = 0; i < loci_table.size(); ++i)
    {
        out_file << loci_table[i].id << sep;
        out_file << loci_table[i].filename << sep;
        out_file << loci_table[i].n_seqs << sep;
        out_file << loci_table[i].last_seq_idx << endl;
    }
    out_file.close();
    
    return SAVE_SUCCESS;
}
// ----------------------------------------------------------------------------

/**
 * @brief      Saves a file in wich each allele index is related the
 *             corresponding locus index (loci-alleles index). For example, in
 *             the lines: 
 *             
 *                 1<tab>2 
 *                 2<tab>3
 *             
 *             The allele sequence 1 corresponds to the locus 2, and the 
 *             allele sequence 2 corresponds to the locus 3
 *
 * @param      output_file_name  The output file name
 *
 * @return     SAVE_SUCCES if the file was saved correctly, otherwise
 *             SAVE_ERROR.
 */
SaveFileResult Indexer::saveLociAllelesIndex(const string & output_file_name)
{
    fstream out_file(output_file_name, ios::out);
    
    if (!out_file)
    {
        cerr << "ERROR: Loci-Alleles index file could not be saved." << endl;
        return SAVE_ERROR;
    }
    
    for (uint64_t i = 0; i < loci_alleles_index.size(); ++i)
    {
        out_file << loci_alleles_index[i] << endl;
    }
    out_file.close();
    
    return SAVE_SUCCESS;
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Saves the allelic profiles file.
 *
 * @param      output_file_name  The output file name
 *
 * @return     SAVE_SUCCES if the file was saved correctly, otherwise SAVE_ERROR.
 */
SaveFileResult Indexer::saveProfilesFile(const string & output_file_name)
{
    string profile_full_path;
    string config_file_path = getPathName(options.getConfigFileFilename());
    
    if (config_file_path.compare("") == 0)
        config_file_path = CURRENT_DIRECTORY;
        
    profile_full_path = pathAppend(config_file_path, 
                                   config_map[PROFILE_CONFIG_KEY][profile_name]);
    
    return copyFile(profile_full_path, output_file_name);
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Saves a file with all the sequence identifiers of each allele of
 *             all of the loci loaded.
 *
 * @param      output_file_name  The output file name
 *
 * @return     SAVE_SUCCES if the file was saved correctly, otherwise
 *             SAVE_ERROR.
 */
SaveFileResult Indexer::saveAlleleSequenceIds(const string & output_file_name)
{
    fstream out_file(output_file_name, ios::out);
    
    if (!out_file)
    {
        cerr << "ERROR: Sequence identifiers file could not be saved" << endl;
        return SAVE_ERROR;
    }
    
    for (uint64_t i = 0; i < length(ids); ++i)
    {
        out_file << ids[i] << endl;
    }
    out_file.close();
    
    return SAVE_SUCCESS;
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Gets the size (total number) of the alleles sequences of the
 *             index.
 *
 * @return     The allele size.
 */
uint64_t Indexer::getSeqsSize()
{
    return length(seqs);
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Creates and saves an ESA index.
 *
 * @param[out] esa_index         The ESA index
 * @param[in]  text              The set of sequences/texts to be indexed
 * @param[in]  output_file_name  The output file name
 *
 * @tparam     T                 Specific type of the ESA index
 * @tparam     U                 Specific type of the sequences/texts to be
 *                               indexed
 *
 * @return     SAVE_SUCCES if the index was saved correctly, otherwise
 *             SAVE_ERROR.
 */
template<typename T, typename U>
SaveFileResult Indexer::createIndex(T &            esa_index, 
                                    const U &      text, 
                                    const string & output_file_name)
{
    // Construct an ESA index from a StringSet
    indexText(esa_index) = text;
    // indexRequire(esa_index, FibreSA());
    indexRequire(esa_index, EsaSA());
    // indexRequire(esa_index, EsaLcp());
    // indexRequire(esa_index, EsaChildtab());  // for TopDown iterators
    // indexRequire(esa_index, EsaBwt());       // for (Super-)MaxRepeats iterators
    
    // Save the index.
    if(! save(esa_index, toCString(output_file_name)))
    {
        cerr << "ERROR: Could not save the index." << endl;
        return SAVE_ERROR;
    }

    return SAVE_SUCCESS;
}
// -----------------------------------------------------------------------------------

/**
 * @brief      Copy the content of a file to another.
 *
 * @param      source  The source file
 * @param      dest    The destination file
 *
 * @return     SAVE_SUCCES if the file was copied correctly, otherwise SAVE_ERROR.
 */
SaveFileResult Indexer::copyFile(const string & source, const string & dest)
{
    ifstream  src(source, std::ios::binary);
    ofstream  dst(dest,   std::ios::binary);
    
    if(!src || !dst) {
        cerr << "ERROR: An error occurred while trying to save the profiles table file." << endl;
        
        if(!src)
            cerr << "ERROR: Profiles table file could not be opened." << endl;

        if(!dst)
            cerr << "ERROR: Profiles table file could not be saved." << endl;
        
        return SAVE_ERROR;
    }
    
    dst << src.rdbuf();
    src.close();
    dst.close();
    return SAVE_SUCCESS ;
}

// -----------------------------------------------------------------------------------
int Indexer::runMlstIndexer()
{
    string out_dir;
    string prefix;
    string out_prefix_filename;
    TCharStringSet profiles;
    
    // Process the config file
    processMlstConfigFile();

    // Load loci files
    loadLociFiles();
    
    cout << endl << "Total sequences loaded: " << getSeqsSize() << endl;
    
    cout << endl << "Loading the profiles file... ";
    // ToDO: Implement this block of code as SIndexer method (in a loadLociFiles fashion)
    // Load profile files
    loadProfilesFile();
    
    // Create profile strings to build an ESA index 
    getProfileStrings(profiles);
    
    out_prefix_filename = options.getOutPrefixFilename();
    relativeDirBaseSplit(out_dir, prefix, out_prefix_filename);
    
    if (out_dir == "")
        out_dir = ".";
    
    // Create the output dir if it does not exist
    if (checkDir(out_dir) != DIRECTORY_EXISTS)
    {
        if (mkpath(out_dir.c_str(), 0750) != 0)
        {
            cerr << "ERROR: Could not create the output directory '" << out_dir << "'." << endl;
            return 1;
        }
    }
    
    cout << "Done!" << endl;
    cout << "Creating and saving ESA index from loaded sequences... ";
    
    // Construct and save all the files that make up the database:
    //  - ESA index from the allele sequences
    //  - Loci table (for each loci: id, file path, n. sequences, and index of the last sequence)
    //  - Loci-Alleles index: a vector with the corresponding locus index (form the Loci Table) for each of the alleles
    //  - Profiles file: a copy of the input allelic profiles file
    //  - Allele sequence ids: a vector with all the alleles sequence ids that make up the ESA index
    
    if(createAllelesIndex(out_prefix_filename) == SAVE_SUCCESS && 
       createProfilesIndex(profiles, options.getOutPrefixFilename() + string(".prof_idx")) == SAVE_SUCCESS && 
       saveLociTable(out_prefix_filename + string(LOCI_TABLE_EXT)) == SAVE_SUCCESS &&
       saveLociAllelesIndex(out_prefix_filename + string(ALLELES_LOCI_INDEX_EXT)) == SAVE_SUCCESS &&
       saveProfilesFile(out_prefix_filename + string(ALLELIC_PROFILE_EXT)) == SAVE_SUCCESS && 
       saveAlleleSequenceIds(out_prefix_filename + string(ALLELES_SEQ_IDS_EXT)) == SAVE_SUCCESS)
    {
        cout << "Done!" << endl;
        cout << "Index created successfully with prefix '" << out_prefix_filename << "'" << endl;
    }
    else
    {
        cerr << "ERROR: Could not save index files." << endl;
        return 1;
    }
        
    return 0;
}
// -----------------------------------------------------------------------------------

int Indexer::runGeneDetectorIndexer()
{
    string out_dir;
    string prefix;
    string out_prefix_filename;
    
    // Process the config file
    processGeneDetectorConfigFile();

    // Load loci files
    loadLociFiles();
    
    cout << endl << "Total loaded sequences: " << getSeqsSize() << endl;
    
    cout << endl << "Creating and saving ESA index from loaded sequences..." << endl;
    
    out_prefix_filename = getOptions().getOutPrefixFilename();
    relativeDirBaseSplit(out_dir, prefix, out_prefix_filename);
    
    if (out_dir == "")
        out_dir = ".";
    
    // Create the output dir if it does not exist
    if (checkDir(out_dir) != DIRECTORY_EXISTS)
    {
        // if (makeDir(out_dir) == MAKE_DIR_ERROR)
        if (mkpath(out_dir.c_str(), 0750) != 0)
        {
            cerr << "ERROR: Could not create the output directory '" << out_dir << "'." << endl;
            return 1;
        }
    }
    
    // Construct ans save all the files that make up the database:
    //  - ESA index from the allele sequences
    //  - Loci table (for each loci: id, file path, n. sequences, and index of the last sequence)
    //  - Loci-Alleles index: a vector with the corresponding locus index (form the Loci Table) for each of the alleles
    //  - Allele sequence ids: a vector with all the alleles sequence ids that meke up the ESA index
    
    if(createAllelesIndex(out_prefix_filename) == SAVE_SUCCESS && 
       saveLociTable(out_prefix_filename + string(LOCI_TABLE_EXT)) == SAVE_SUCCESS &&
       saveLociAllelesIndex(out_prefix_filename + string(ALLELES_LOCI_INDEX_EXT)) == SAVE_SUCCESS &&
       saveAlleleSequenceIds(out_prefix_filename + string(ALLELES_SEQ_IDS_EXT)) == SAVE_SUCCESS)
    {
        cout << "Done!" << endl;
        cout << "Index created successfully with prefix '" << out_prefix_filename << "'" << endl;
    }
    else
    {
        cerr << "ERROR: Could not save index files." << endl;
        return 1;
    }
        
    return 0;
}
// -----------------------------------------------------------------------------------
int Indexer::run(int argc, char const ** argv)
{
    setupArgumentParser();

    // Parse command line arguments and options
    ArgumentParser::ParseResult res = parseCommandLine(argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Check if the input file exists
    if (checkFile(options.getConfigFileFilename()) == FILE_NOT_FOUND)
        return 1;
    
    if (options.getMode() == MLST_MODE)
        return runMlstIndexer();
    else if (options.getMode() == GENE_DETECTOR_MODE)
        return runGeneDetectorIndexer();
    
    return 1;
}
// -----------------------------------------------------------------------------------
// Getters and setters
// -----------------------------------------------------------------------------------


/**
 * @brief      Gets the alleles esa index.
 *
 * @return     The alleles esa index.
 */
const TSequenceEsaIndex & Indexer::getAllelesEsaIndex() const
{
    return alleles_esa_index;
}


/**
 * @brief      Sets the alleles esa index.
 *
 * @param[in]  alleles_esa_index  The alleles esa index
 */
void Indexer::setAllelesEsaIndex(const TSequenceEsaIndex & alleles_esa_index)
{
    this->alleles_esa_index = alleles_esa_index;
}


/**
 * @brief      Gets the config map.
 *
 * @return     The config map.
 */
const TConfigMap & Indexer::getConfigMap() const
{
    return config_map;
}


/**
 * @brief      Sets the config map.
 *
 * @param[in]  config_map  The config map
 */
void Indexer::setConfigMap(const TConfigMap & config_map)
{
    this->config_map = config_map;
}


/**
 * @brief      Gets the sequence identifiers of the alleles.
 *
 * @return     The sequence identifiers.
 */
const TCharStringSet & Indexer::getIds() const
{
    return ids;
}


/**
 * @brief      Sets the sequence identifiers of the alleles.
 *
 * @param[in]  ids   The sequence identifiers
 */
void Indexer::setIds(const TCharStringSet & ids)
{
    this->ids = ids;
}


/**
 * @brief      Gets the loci-alleles index.
 *
 * @return     The loci-alleles index.
 */
TUintVector Indexer::getLociAllelesIndex() const
{
    return loci_alleles_index;
}


/**
 * @brief      Sets the loci-alleles index.
 *
 * @param[in]  loci_alleles_index  The loci-alleles index
 */
void Indexer::setLociAllelesIndex(TUintVector loci_alleles_index)
{
    this->loci_alleles_index = loci_alleles_index;
}


/**
 * @brief      Gets the loci table.
 *
 * @return     The loci table.
 */
const TLociTable & Indexer::getLociTable() const
{
    return loci_table;
}


/**
 * @brief      Sets the loci table.
 *
 * @param[in]  loci_table  The loci table
 */
void Indexer::setLociTable(const TLociTable & loci_table)
{
    this->loci_table = loci_table;
}


/**
 * @brief      Gets the Indexer program options.
 *
 * @return     The options.
 */
const IndexerOptions & Indexer::getOptions() const
{
    return options;
}


/**
 * @brief      Sets the Indexer program options.
 *
 * @param[in]  options  The options
 */
void Indexer::setOptions(const IndexerOptions & options)
{
    this->options = options;
}


/**
 * @brief      Gets the command line parser of the Indexer program.
 *
 * @return     The parser.
 */
const ArgumentParser & Indexer::getParser() const
{
    return parser;
}


/**
 * @brief      Sets the command line parser of the Indexer program.
 *
 * @param[in]  parser  The parser
 */
void Indexer::setParser(const ArgumentParser & parser)
{
    this->parser = parser;
}


/**
 * @brief      Gets the profile name.
 *
 * @return     The profile name.
 */
const string & Indexer::getProfileName() const
{
    return profile_name;
}


/**
 * @brief      Sets the profile name.
 *
 * @param[in]  profile_name  The profile name
 */
void Indexer::setProfileName(const string & profile_name)
{
    this->profile_name = profile_name;
}


/**
 * @brief      Gets the profiles ESA index.
 *
 * @return     The profiles ESA index.
 */
const TCharEsaIndex & Indexer::getProfilesEsaIndex() const
{
    return profiles_esa_index;
}


/**
 * @brief      Sets the profiles ESA index.
 *
 * @param[in]  profiles_esa_index  The profiles ESA index
 */
void Indexer::setProfilesEsaIndex(const TCharEsaIndex & profiles_esa_index)
{
    this->profiles_esa_index = profiles_esa_index;
}


/**
 * @brief      Gets the profiles table.
 *
 * @return     The profiles table.
 */
const ProfilesTable & Indexer::getProfilesTable() const
{
    return profiles_table;
}


/**
 * @brief      Sets the profiles table.
 *
 * @param[in]  profiles_table  The profiles table
 */
void Indexer::setProfilesTable(const ProfilesTable & profiles_table)
{
    this->profiles_table = profiles_table;
}


/**
 * @brief      Gets the program name.
 *
 * @return     The program name.
 */
const string & Indexer::getProgramName() const
{
    return program_name;
}


/**
 * @brief      Sets the program name.
 *
 * @param[in]  program_name  The program name
 */
void Indexer::setProgramName(const string & program_name)
{
    this->program_name = program_name;
}


/**
 * @brief      Gets the sequences.
 *
 * @return     The sequences.
 */
const TStringSet & Indexer::getSeqs() const
{
    return seqs;
}


/**
 * @brief      Sets the sequences.
 *
 * @param[in]  seqs  The sequences
 */
void Indexer::setSeqs(const TStringSet & seqs)
{
    this->seqs = seqs;
}


/**
 * @brief      Destroys the object.
 */
Indexer::~Indexer()
{
    
}

// ============================================================================
// Main function
// ============================================================================

/**
 * @brief      Creates an Indexer object and runs the indexer application
 *
 * @param[in]  argc  The number of input arguments
 * @param      argv  The list of input arguments
 *
 * @return     0 if the program ends correctly, 1 otherwise
 */
int main(int argc, char const ** argv)
{
    
    Indexer indexer_app;
    return indexer_app.run(argc, argv);
}
