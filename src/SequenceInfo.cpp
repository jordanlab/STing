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


#include "SequenceInfo.h"

using namespace std;

SequenceInfo::SequenceInfo()
{
    idx              = -1;
    desc             = "";
    allele_number    = "NA";
    length           = 0;
    kmerCount        = 0;
    normKmerCount    = 0.0;
    avgKmerDepth     = 0.0;
    kmerDepthSd      = 0.0;
    score            = 0.0;
    perBaseDepths    = {};
    nBasesCovered    = 0;
    kmerDepthSdSet   = false;
}

SequenceInfo::SequenceInfo(const SequenceInfo& obj) :
    idx(obj.idx),
    desc(obj.desc),
    allele_number(obj.allele_number),
    length(obj.length),
    kmerCount(obj.kmerCount),
    normKmerCount(obj.normKmerCount),
    avgKmerDepth(obj.avgKmerDepth),
    kmerDepthSd(obj.kmerDepthSd),
    score(obj.score),
    perBaseDepths(obj.perBaseDepths),
    nBasesCovered(obj.nBasesCovered),
    kmerDepthSdSet(obj.kmerDepthSdSet)
{
    
}

SequenceInfo::SequenceInfo(const int&      idx, 
                           const string&   desc, 
                           const uint64_t& length)
{
    this->idx           = idx;
    this->desc          = desc;
    this->length        = length;
    this->allele_number = "";
    this->kmerCount     = 0;
    this->normKmerCount = 0.0;
    this->avgKmerDepth  = 0.0;
    this->kmerDepthSd   = 0.0;
    this->score         = 0.0;
    this->perBaseDepths = {};
    this->nBasesCovered = 0;
    
    this->kmerDepthSdSet = false;
}

SequenceInfo::SequenceInfo(const int&      idx, 
                           const char*     desc, 
                           const uint64_t& length)
{
    SequenceInfo(idx, string(desc), length);
}

SequenceInfo::~SequenceInfo()
{

}

bool SequenceInfo::calculateKmerDepthSd()
{
    double sum     = 0.0; 
    double mean    = 0.0;
    double sumDifs = 0.0;
    int dataSize   = perBaseDepths.size();
    
    if (dataSize > 0)
    {
        for (auto const& x: perBaseDepths)
            sum += x;
        
        mean = (double)(sum / dataSize);
        
        for (auto const& x: perBaseDepths)
            sumDifs += pow((x - mean), 2);
        
        kmerDepthSd      = sqrt(sumDifs / (dataSize - 1));
        kmerDepthSdSet = true;
        return true;
    }
    return false;
}

void SequenceInfo::calculateScore()
{
    double sd_trans = 1;
    double cov      = getCoverage();
    if (kmerDepthSdSet)
        sd_trans = 1/kmerDepthSd;
    if (cov == 0.0)
        cov = 1;
    score = normKmerCount * (cov * 0.1) * sd_trans;
}

// double SequenceInfo::getNormKmerCount() const
// {
//     if (length > 0)
//         return (double)((double)kmerCount / (double)length);
//     return 0;
// }

double SequenceInfo::getCoverage() const
{
    if (length > 0)
        return (double)((double)nBasesCovered / (double)length * 100);
    return 0.0;
}

string SequenceInfo::getCoverageStr(const uint64_t & precision) const
{
    stringstream out;
    out << nBasesCovered 
        << "/" 
        << length 
        << "(" 
        << fixed
        << std::setprecision(precision) 
        << getCoverage() 
        << "%)";
    return out.str();
}

bool SequenceInfo::operator<(const SequenceInfo & other)
{
    return (getNormKmerCount() < other.getNormKmerCount());
}

bool SequenceInfo::operator>(const SequenceInfo & other)
{
    return (getNormKmerCount() > other.getNormKmerCount());
}

ostream& operator<<(ostream & out, const SequenceInfo & object)
{
    
    out << endl;
    out << "idx:                " << object.getIdx() << endl;
    out << "desc:               " << object.getDesc() << endl;
    out << "allele_number:      " << object.getAlleleNumber() << endl;
    out << "length:             " << object.getLength() << endl;
    out << "kmerCount:          " << object.getKmerCount() << endl;
    out << "nomrKmerCount:      " << object.getNormKmerCount() << endl;
    out << "nBasesCovered:      " << object.getNBasesCovered() << endl;
    out << "avgKmerDepth:       " << object.getAvgKmerDepth() << endl;
    out << "kmerDepthSd:        " << object.getKmerDepthSd() << endl;
    out << "score:              " << object.getScore() << endl;
    out << "coverage:           " << object.getCoverageStr() << endl;
    out << "perBaseDepths size: " << object.getPerBaseDepths().size() << endl;
    return out;
}

int SequenceInfo::getIdx() const
{
    return idx;
}

string SequenceInfo::getDesc() const
{
    return desc;
}

string SequenceInfo::getAlleleNumber() const
{
    return allele_number;
}

uint64_t SequenceInfo::getLength() const
{
    return length;
}

uint64_t SequenceInfo::getKmerCount() const
{
    return kmerCount;
}

double SequenceInfo::getNormKmerCount() const
{
    return normKmerCount;
}

double SequenceInfo::getAvgKmerDepth() const
{
    return avgKmerDepth;
}

double SequenceInfo::getKmerDepthSd() const
{
    return kmerDepthSd;
}

double SequenceInfo::getScore() const
{
    return score;
}

vector<uint64_t> SequenceInfo::getPerBaseDepths() const
{
    return perBaseDepths;
}

uint64_t SequenceInfo::getNBasesCovered() const
{
    return nBasesCovered;
}

bool SequenceInfo::isKmerDepthSdSet() const
{
    return kmerDepthSdSet;
}

void SequenceInfo::setIdx(const int &  valIdx)
{
    idx = valIdx;
}

void SequenceInfo::setDesc(const string & valDesc)
{
    desc = valDesc;
}

void SequenceInfo::setAlleleNumber(const string & valAlleleNumber)
{
    allele_number = valAlleleNumber;
}

void SequenceInfo::setLength(const uint64_t & valLength)
{
    length = valLength;
}

void SequenceInfo::setKmerCount(const uint64_t & valKmerCount)
{
    kmerCount = valKmerCount;
}

void SequenceInfo::setNormKmerCount(const double & valNormKmerCount)
{
    normKmerCount = valNormKmerCount;
}

void SequenceInfo::setAvgKmerDepth(const double & valAvgKmerDepth)
{
    avgKmerDepth = valAvgKmerDepth;
}

void SequenceInfo::setKmerDepthSd(const double & valKmerDepthSd)
{
    kmerDepthSdSet = true;
    kmerDepthSd      = valKmerDepthSd;
}

void SequenceInfo::setScore(const double & valScore)
{
    score = valScore;
}

void SequenceInfo::setPerBaseDepths(const vector<uint64_t> & valPerBaseDepths)
{
    perBaseDepths = valPerBaseDepths;
}

void SequenceInfo::setNBasesCovered(const uint64_t & valNBasesCovered)
{
    nBasesCovered = valNBasesCovered;
}

