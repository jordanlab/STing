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


#ifndef TYPES_HPP
#define TYPES_HPP

#include <unordered_map>


#include <seqan/misc/interval_tree.h>

using namespace seqan;
using namespace std;

// ============================================================================
// Forward declarations (cross referencing solving)
// ============================================================================
struct SLocusRecord; 
struct SeqCoverage;
struct SAlleleCount;
class SequenceInfo; 
// ============================================================================
// Types
// #F8B195, #F67280, #C06C84, #6C5B7B, #355C7D
// ============================================================================
typedef String<Dna5>                               TSequence;
// typedef String<Iupac>                              TSequence;
typedef StringSet<TSequence>                       TStringSet;
typedef Index<TStringSet, IndexEsa<> >             TSequenceEsaIndex;
typedef StringSet<CharString>                      TCharStringSet;
typedef Index<TCharStringSet, IndexEsa<> >         TCharEsaIndex;
typedef Finder<TSequenceEsaIndex>                  TSequenceEsaIndexFinder;
typedef Finder<TCharEsaIndex>                      TCharEsaIndexFinder;
typedef Position<TSequenceEsaIndexFinder >::Type   TSequenceEsaIndexPosition;
typedef vector<string>                             TStringVector;
typedef TStringVector::iterator                    TStringVectorIterator;
// typedef FragmentStore<>                            TStore;
typedef pair<uint64_t, uint64_t>                   TSequenceCountPair;
typedef map<string, map<string, string > >         TConfigMap;
typedef map<string, string >                       TLociFilesMap;
typedef vector<SLocusRecord>                       TLociTable;
typedef vector<uint64_t>                           TUintVector;
typedef vector<bool>                               TBoolVector;
typedef vector<double>                             TDoubleVector;
typedef vector<SeqCoverage>                        TSeqCovVector;
typedef vector<SequenceInfo>                       TSeqInfoVector;
typedef vector<SAlleleCount>                       TAlleleCountVector;
typedef IntervalAndCargo<uint64_t, uint64_t>       TInterval;
typedef vector<TInterval>                          TIntervalVector;
typedef unordered_map<uint64_t, uint64_t >         TUMapUintUint;
typedef unordered_map<uint64_t, TUMapUintUint >    TUMapUintUMapUintUint;
typedef unordered_map<uint64_t, TUintVector >      TUMapUintUintVector;
typedef unordered_map<string, TUintVector >        TUMapStringUintVector;
typedef unordered_map<uint64_t, double >           TUMapUintDouble;
typedef unordered_map<uint64_t, TIntervalVector >  TUMapIntervalVector;
typedef unordered_map<string, TAlleleCountVector > TUMapAlleleCountVector;
typedef unordered_map<string, TSeqInfoVector >     TUMapSeqInfoVector;
typedef unordered_map<string, SequenceInfo >       TUMapSequenceInfo;
// typedef unordered_map<uint64_t, TBSTUint* >      TUMapUintBST;
// typedef BST<uint64_t>                            TBSTUint;

typedef unordered_map<uint64_t, TUintVector >::iterator     TUMapUintUintVectorIt;
// ============================================================================
// Enumerations
// ============================================================================
enum SaveFileResult
{
    SAVE_ERROR   = 0,
    SAVE_SUCCESS = 1
};
// ----------------------------------------------------------------------------
enum CheckFileResult
{
    FILE_NOT_FOUND = 0,
    FILE_EXISTS    = 1
};
// ----------------------------------------------------------------------------
enum MakeDirResult
{
    MAKE_DIR_ERROR   = 0,
    MAKE_DIR_SUCCESS = 1
};
// ----------------------------------------------------------------------------
enum CheckDirResult
{
    DIRECTORY_NOT_ACCESSIBLE = 0,
    IS_NOT_DIRECTORY         = 1,
    DIRECTORY_EXISTS         = 2
};
// ----------------------------------------------------------------------------
enum LogLevel
{
    LOG_OFF = 0,
    LOG_ON  = 1
};
// ----------------------------------------------------------------------------
enum IndexerMode
{
    MLST_MODE          = 0,
    GENE_DETECTOR_MODE = 1
};
// ----------------------------------------------------------------------------
// enum StStatus
// {
//     ST_PREDICTED    = 0,
//     NO_KMER_HITS    = 1,
//     NO_ST_IN_TABLE  = 2
// };
// ============================================================================
// Macros
// ============================================================================
// Utilities:
// - Stringfication
#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

// General values
#define FIELD_SEPARATOR          "\t"
#define NA_STRING                "NA"
#define LOW_CONFIDENCE_SYMBOL    "*"
#define CURRENT_DIRECTORY        "."

// Default input param values
#define DEFAULT_KMER_LENGHT      30
#define DEFAULT_TOP_N_ALLELES    2

// Values to check gaps in k-mer depth
#define BASES_MARGIN             10          // from the 10th base to the length of allele-10
#define DEPTH_DIFF_MAX_PEAK      1.414214    // 2^0.5
#define DEPTH_DIFF_MIN_PEAK      0.7071068   // 2^-0.5

// Allele sequence separator
#define ALLELE_SEQ_ID_SEPARATOR   "_"

// Config file related values:
// - Config file comment character
#define COMMENT_CHAR              '#'

// - Config file section headers and keys
#define LOCI_CONFIG_HEADER        "[loci]"
#define LOCI_CONFIG_KEY           "loci"
#define PROFILE_CONFIG_KEY        "profile"
#define PROFILE_CONFIG_HEADER     "[profile]"

// Valid values for the option mode (-m|--mode) in the indexer app
#define MLST_MODE_OPTION          "MLST"
#define GENE_DETECTOR_MODE_OPTION "GDETECT"

// File extensions (index)
#define ESA_INDEX_EXT             ".sa"
#define LOCI_TABLE_EXT            ".dat"
#define ALLELES_LOCI_INDEX_EXT    ".ali"
#define ALLELIC_PROFILE_EXT       ".apr"
#define ALLELES_SEQ_IDS_EXT       ".ids"
#define PROFILE_INDEX_SUFFIX      ".prof_idx"

// App info
#define INDEXER_APP_VERSION       "0.24.2"
#define DETECTOR_APP_VERSION      "0.24.2"
#define TYPER_APP_VERSION         "0.24.2"
#define INDEXER_APP_UPDATE        "01/09/2019"
#define DETECTOR_APP_UPDATE       "03/17/2019"
#define TYPER_APP_UPDATE          "04/03/2019"

#define INDEXER_APP_NAME          "indexer"
#define TYPER_APP_NAME            "typer"
#define DETECTOR_APP_NAME         "detector"

#endif // TYPES_HPP
