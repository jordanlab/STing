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


#ifndef _CUSTOM_UTILS_H_
#define _CUSTOM_UTILS_H_

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

// #include "types.hpp"

using namespace std;
// =========================================================================
// split a string
struct Split
{
  enum EmptiesT { EMPTIES_OK, NO_EMPTIES };
};

template <typename Container>
inline Container& split(
  Container&                            result,
  const typename Container::value_type& s,
  const typename Container::value_type& delimiters,
  Split::EmptiesT                       empties = Split::EMPTIES_OK )
{
  result.clear();
  size_t current;
  size_t next = -1;
  do
  {
    if (empties == Split::NO_EMPTIES)
    {
      next = s.find_first_not_of( delimiters, next + 1 );
      if (next == Container::value_type::npos) break;
      next -= 1;
    }
    current = next + 1;
    next = s.find_first_of( delimiters, current );
    result.push_back( s.substr( current, next - current ) );
  }
  while (next != Container::value_type::npos);
  return result;
}
// =========================================================================
// print a StringSet
// template<typename T>
inline void printStringSet(const TCharStringSet & v, string sep, std::ostream & obj_ostream)
{
    for(uint64_t i = 0; i < length(v); ++i){
        if (i == length(v) - 1 ){
            obj_ostream << v[i];
        }
        else
            obj_ostream << v[i] << sep;
    }
    // obj_ostream << endl;
}
// =========================================================================
// print a std::vector
template<typename T>
inline void printVector(vector<T> & v)
{
    // typename vector<T>::iterator it;
    // for (it = v.begin(); it != v.end(); ++it)
    //     cout << *it << " ";
    string delim = "";
    for (auto x: v)
    {
        cout << delim << x;
        delim = " ";
    }
    
    cout << endl;
}
// =========================================================================
template<typename T>
inline void printVector(vector<T> & v, string sep)
{
    // typename vector<T>::iterator it;
    // for (it = v.begin(); it != v.end(); ++it)
    //     cout << *it << sep;
    string delim = "";
    for (auto x: v)
    {
        cout << delim << x;
        delim = sep;
    }
    
    cout << endl;
}
// =========================================================================
template<typename T>
inline void printVector(vector<T> & v, string sep, std::ostream & obj_ostream)
{
    // typename vector<T>::iterator it;
    // for (it = v.begin(); it != v.end(); ++it)
    //     obj_ostream << *it << sep;
    string delim = "";
    for (auto x: v)
    {
        obj_ostream << delim << x;
        delim = sep;
    }
}
// =========================================================================
template<typename T>
inline void printVector(vector<T> const & v, string sep, std::ostream & obj_ostream)
{
    string delim = "";
    for (auto x: v) {
        obj_ostream << delim << x;
        delim = sep;
    }
}
// =========================================================================
inline string getFileName(const string & s) 
{
    char sep = '/';

#ifdef _WIN32
    sep = '\\';
#endif

    size_t i = s.rfind(sep, s.length());
    if (i != string::npos) {
        return(s.substr(i+1, s.length() - i));
    }

    return(s);
}
// =========================================================================
inline string getPathName(const string & s)
{
    char sep = '/';
    string ret_path = "";

#ifdef _WIN32
    sep = '\\';
#endif
    
    // size_t i = s.rfind(sep, s.length());
    // if (i == 0)  //    /file
    // {
    //     ret_path = "/"; // root in nix systems
    // }
    // else if(i == string::npos) //    file
    // {
    //     ret_path = "."; // current directory in nix systems
    // }
    // else // ./file, /path/to/file, path/to/file
    // {
    //     ret_path = s.substr(0, i);
    // }
    
    
    size_t i = s.rfind(sep, s.length());
    if (i != string::npos) {
       ret_path = s.substr(0, i);
    }

    // return(s);
    
    return(ret_path);
}
// =========================================================================
inline string pathAppend(const string & p1, const string & p2)
{
    char sep   = '/';
    string tmp = p1;
    
#ifdef _WIN32
    sep = '\\';
#endif
    
    regex reg_exp(string(1, sep) + "{2,}");
    // string ret_path = "";
    
    // cout << "p1[p1.length()-1]"<< p1[p1.length()-1] << endl;
    if (p1[p1.length() - 1] != sep) { // Need to add a
        tmp += sep;                // path separator
        // ret_path = tmp + p2;
        return(tmp + p2);
    }
    else
        return(p1 + p2);
        // ret_path = p1 + p2;
    
    // Remove double separator occurrences
    // ret_path = regex_replace(ret_path, reg_exp, string(1, sep));
    // return(ret_path);
}
// // =========================================================================
// inline string getFileName(const string & s) {

//    string sep = "/\\";

//    size_t i = s.rfind(sep, s.length());
//    if (i != string::npos) {
//       return(s.substr(i+1, s.length() - i));
//    }

//    return(s);
// }
// // =========================================================================
// inline string getPathName(const string & s) {

//    string sep = "/\\";

//    size_t i = s.rfind_last_of(sep, s.length());
//    if (i != string::npos) {
//       return(s.substr(0, i));
//    }

//    return(s);
// }
// ==========================================================================
inline TStringVector mergeOptionValues(const TStringVector & option_values)
{
    TStringVector files;
    for (uint64_t i = 0; i < option_values.size(); ++i)
    {
        TStringVector fields;
        split(fields, option_values[i], ",");
        
        for (uint64_t j = 0; j < fields.size(); ++j)
            files.push_back(fields[j]);
    }
    
    return files;
}
// ==========================================================================
inline MakeDirResult makeDir(const string & path)
{
    mode_t n_mode = 0755; // UNIX style permissions
    int n_error = 0;
    
#if defined(_WIN32)
    n_error = _mkdir(path.c_str()); // can be used on Windows
#else 
    n_error = mkdir(path.c_str(), n_mode); // can be used on non-Windows
#endif
    if (n_error != 0)
    {
      cerr << "ERROR: " << strerror(errno) << endl;
      return MAKE_DIR_ERROR;
    }

    return MAKE_DIR_SUCCESS;
}
// ==========================================================================
inline CheckDirResult checkDir(const string & path)
{
    struct stat info;

    if(stat(path.c_str(), &info ) != 0)
        return DIRECTORY_NOT_ACCESSIBLE;
    
    else if(info.st_mode & S_IFDIR)
        return DIRECTORY_EXISTS;
    
    return IS_NOT_DIRECTORY;
}
// ==========================================================================
inline CheckFileResult checkFile(const string & filename)
{
    fstream infile(filename, ios::binary | ios::in);
    if (!infile.good())
    {
        cerr << "ERROR: Could not open the input file '" << filename << "'\n";
        return FILE_NOT_FOUND;
    }
    
    infile.close();
    return FILE_EXISTS;
}
// ==========================================================================
inline CheckFileResult checkFiles(const TStringVector & files)
{
    CheckFileResult result = FILE_EXISTS;
    
    for (uint64_t i = 0; i < files.size(); ++i)
    {
        if(checkFile(files[i]) == FILE_NOT_FOUND)
            result = FILE_NOT_FOUND;
    }
    
    return result;
}
// =========================================================================
// Taken from https://nakkaya.com/2009/11/08/command-line-progress-bar/
// =========================================================================
inline void printProgressBar(uint64_t percent)
{
    std::string bar;

    for(uint64_t i = 0; i < 50; i++){
        if( i < (percent/2))
        {
            bar.replace(i,1,"=");
        }
        else if( i == (percent/2))
        {
          bar.replace(i,1,">");
        }
        else
        {
          bar.replace(i,1," ");
        }
    }

    cout<< "\r" "[" << bar << "] ";
    cout.width( 3 );
    cout<< percent << "%     " << flush;
}
// =========================================================================
// Print the counts map in a descendent sorted way (by count)

typedef std::pair<string, uint64_t> TCountPair;
// -------------------------------------------------------------------------
// Auxiliary function to sort
inline bool sortPairDesc(TCountPair i, TCountPair j) { 
    return i.second > j.second; 
}
// -------------------------------------------------------------------------
// Print function vector of pairs (seq_id, count)
inline void printCountsMap(vector<TCountPair> & counts_vector, 
                           const string &       separator,
                           const bool &         sorted = true)
{   
    // Sort the vector
    if (sorted)
        sort(counts_vector.begin(), counts_vector.end(), sortPairDesc);
    
    // print the sorted vector
    cout << "Seq" << separator << "Count" << endl;
    for (TCountPair i : counts_vector)
        cout << i.first << separator << i.second << endl;
}
// -------------------------------------------------------------------------
// Print function map (seq_id -> count)
inline void printCountsMap(map<string, uint64_t> const & map, 
                           const string &                separator)
{   
    // Wrap the elements from the counts map to a vector of pairs
    std::vector<TCountPair> pairs_vector;
    for (const auto& i : map)
    {
        TCountPair p(i.first, i.second);
        pairs_vector.push_back(p);
    }
    
    // Sort the vector
    sort(pairs_vector.begin(), pairs_vector.end(), sortPairDesc);
    
    // print the sorted vector
    cout << "Seq" << separator << "Count" << endl;
    for (TCountPair i : pairs_vector)
        cout << i.first << separator << i.second << endl;
}

// ==========================================================================
inline TStringVector mergeInValues(const CharString & values)
{
    TStringVector files;
    TStringVector fields;
    
    split(fields, toCString(values), ",");
        
    for (uint64_t j = 0; j < fields.size(); ++j)
        files.push_back(fields[j]);
    
    return files;
}

// =========================================================================
inline void loadLociTable(vector<SLocusRecord> & loci_table, const string & filename)
{
    fstream loci_file(filename, ios::in);
    
    if (! loci_file.good()) 
    {
        cerr << "ERROR: Opening file '" << filename << "' failed.\n";
        exit(1);
    }
    
    string line;
    while(getline(loci_file, line))
    {   
        vector<string> fields;
        if (!line.empty() && line.at(0) != '#') // The current line is not an empty/comment line
        {
          // cout << line << endl;
          split(fields, line, "\t");
          // SLocusRecord loci = {fields[0], fields[1], fields[2], (uint64_t)stoi(fields[3]), (uint64_t)stoi(fields[4])};
          SLocusRecord loci = {fields[0], fields[1], (uint64_t)stoi(fields[2]), (uint64_t)stoi(fields[3])};
          loci_table.push_back(loci);
        }
    }
}
// =========================================================================
inline void loadProfilesFile(ProfilesTable & profiles_table, 
                      const string &  filename)
{
    unordered_map<string, vector<string>> table;
    fstream profiles_file(filename, ios::in);
    
    if (profiles_file.is_open()) 
    {
        string line;
        vector<string> col_names;
        uint64_t current_line = 0;
        while(getline(profiles_file, line))
        {
            vector<string> fields;
            
            if (!line.empty() && line.at(0) != '#') // The current line is not an empty/comment line
            {
                split(fields, line, "\t");
                
                if (++current_line == 1) // get header and create the map
                {   
                    profiles_table.cols = fields;
                    // col_names = fields;
                    // for(uint64_t i = 0; i < fields.size(); ++i) {
                    for(string col : profiles_table.cols) {
                        table[col] = vector<string>();
                        // profiles_table.table[col] = vector<string>();
                        profiles_table.table[col] = table[col];
                    }
                }
                else {
                    for(uint64_t i = 0; i < profiles_table.cols.size(); ++i) {
                        // table[col_name[i]s[i]].push_back(fields[i]);
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
// =========================================================================
inline void loadProfilesFile(ProfilesTable & profiles_table, 
                      TConfigMap &    config_map)
{
    unordered_map<string, vector<string>> table;
    
    string filename;
    
    for(auto const & pofile_item : config_map[PROFILE_CONFIG_KEY]) {
        filename = pofile_item.second;
    }
    
    fstream profiles_file(filename, ios::in);
    
    if (profiles_file.is_open()) 
    {
        string line;
        vector<string> col_names;
        uint64_t current_line = 0;
        while(getline(profiles_file, line))
        {
            vector<string> fields;
            
            if (!line.empty() && line.at(0) != COMMENT_CHAR) // The current line is not an empty/comment line
            {
                split(fields, line, "\t");
                
                if (++current_line == 1) // get header and create the map
                {   
                    profiles_table.cols = fields;
                    // col_names = fields;
                    // for(uint64_t i = 0; i < fields.size(); ++i) {
                    for(string col : profiles_table.cols) {
                        table[col] = vector<string>();
                        // profiles_table.table[col] = vector<string>();
                        profiles_table.table[col] = table[col];
                    }
                }
                else {
                    for(uint64_t i = 0; i < profiles_table.cols.size(); ++i) {
                        // table[col_name[i]s[i]].push_back(fields[i]);
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
// =========================================================================
inline void printProfilesTable(ProfilesTable & profiles_table)
{
    // print header
    for(string col : profiles_table.cols) {
        cout << col << "\t" << endl;
    }
    // print table content
    for (uint64_t i = 0; i < profiles_table.table[profiles_table.cols[0]].size(); ++i)
    {
        for(string col : profiles_table.cols) {
            cout << profiles_table.table[col][i] << "\t";
        }
        cout << endl;
    }
}
// =========================================================================
inline void printSLocusRecord(SLocusRecord & locus_record)
{
    cout << "id" 
         << "\t" 
         << "filename" 
         << "\t" 
         << "n_seqs" 
         << "\t" 
         << "last_seq_idx" 
         << endl;
         
    cout << locus_record.id 
         << "\t" 
         << locus_record.filename
         << "\t" 
         << locus_record.n_seqs
         << "\t" 
         << locus_record.last_seq_idx
         << endl;
}
// =========================================================================
// std::string revComp(std::string seq)
inline void revComp(std::string & seq)
{
    auto lambda = [](const char c) {
        switch (c) {
        case 'A':
            return 'T';
        case 'a':
            return 't';
        case 'T':
            return 'A';
        case 't':
            return 'a';
        case 'G':
            return 'C';
        case 'g':
            return 'c';    
        case 'C':
            return 'G';
        case 'c':
            return 'g';
        case 'N':
            return 'N';
        case 'n':
            return 'n';
        default:
            throw std::domain_error("Invalid nucleotide.");
        }
    };

    std::transform(seq.begin(), seq.end(), seq.begin(), lambda);
    std::reverse(seq.begin(), seq.end());
    // return seq;
}
// =========================================================================

static const char basemap[255] =
    {
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
};

inline void revCompBasemap(std::string & seq)
{
    auto lambda = [](const char c) {
        return basemap[int(c)];
    };

    std::transform(seq.begin(), seq.end(), seq.begin(), lambda);
    std::reverse(seq.begin(), seq.end());
    // return seq;
}

// =========================================================================
inline SaveFileResult copyProfilesFile(const string & source, const string & dest)
{
    std::ifstream  src(source, std::ios::binary);
    std::ofstream  dst(dest,   std::ios::binary);
    
    if(!src || !dst) {
        cerr << "ERROR: An error occurred while triying to save the profiles table file" << endl;
        
        if(!src)
            cerr << "ERROR: Profiles table file could not be opened" << endl;

        if(!dst)
            cerr << "ERROR: Profiles table file could not be saved" << endl;
        
        return SAVE_ERROR;
    }
    
    dst << src.rdbuf();
    return SAVE_SUCCESS ;
}

// =========================================================================
template <typename T>
inline std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << setprecision(n) << a_value;
    return out.str();
}
// =========================================================================
// taken from http://stackoverflow.com/questions/11208971/round-a-float-to-a-given-precision
template <typename T>
inline T pround(T const x, int const precision)
{
    if (x == 0.)
        return x;
    int ex = floor(log10(abs(x))) - precision + 1;
    T div = pow(10, ex);
    return floor(x / div + 0.5) * div;
}

// =========================================================================
template <typename T>
inline std::string getNumberStringWithPrecision(T const value, int const precision = 2) 
{
    stringstream stream;
    stream << fixed << setprecision(precision) << value;
    return stream.str();
}
// =========================================================================
inline bool greaterSAlleleCount(const SAlleleCount & a, const SAlleleCount & b)
{
    return (a.count > b.count);
}

#endif
