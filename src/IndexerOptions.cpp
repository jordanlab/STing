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


#include "IndexerOptions.h"
using namespace std;

// -----------------------------------------------------------------------------------
// Constructors /destructors
// -----------------------------------------------------------------------------------

/**
 * @brief      Constructs the object.
 */
IndexerOptions::IndexerOptions()
{
	config_file_filename = "";
	out_prefix_filename = "";
}

IndexerOptions::IndexerOptions(const IndexerOptions & indexer_options_obj)
{
	config_file_filename = indexer_options_obj.getConfigFileFilename();
	out_prefix_filename  = indexer_options_obj.getOutPrefixFilename();
}

IndexerOptions::~IndexerOptions()
{

}
// -----------------------------------------------------------------------------------
// Class methods (different from accessors)
// -----------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------
// Getters and setters
// -----------------------------------------------------------------------------------
string IndexerOptions::getConfigFileFilename() const
{
	return config_file_filename;
}

void IndexerOptions::setConfigFileFilename(const string & config_file_filename)
{
	this->config_file_filename = config_file_filename;
}

string IndexerOptions::getOutPrefixFilename() const
{
	return out_prefix_filename;
}

void IndexerOptions::setOutPrefixFilename(const string & out_prefix_filename)
{
    this->out_prefix_filename = out_prefix_filename;
}

IndexerMode IndexerOptions::getMode() const
{
    return mode;
}

void IndexerOptions::setMode(const IndexerMode & mode)
{
    this->mode = mode;
}