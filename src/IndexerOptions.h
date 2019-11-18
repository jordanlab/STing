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


#ifndef SINDEXEROPTIONS_H_
#define SINDEXEROPTIONS_H_

#include <string>

#include <seqan/sequence.h>
#include <seqan/index.h>

#include "types.hpp"

using namespace std;

class IndexerOptions
{
		string      config_file_filename;
		string      out_prefix_filename;
		IndexerMode mode;

	public:
		IndexerOptions();
		IndexerOptions(const IndexerOptions & indexer_options_obj);
		virtual ~IndexerOptions();
		string getConfigFileFilename() const;
		void setConfigFileFilename(const string & config_file_filename);
		string getOutPrefixFilename() const;
		void setOutPrefixFilename(const string & out_prefix_filename);
		IndexerMode getMode() const;
		void setMode(const IndexerMode & mode);
};

#endif /* SINDEXEROPTIONS_H_ */
