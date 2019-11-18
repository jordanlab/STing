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


#ifndef SEQUENCEINFO_H
#define SEQUENCEINFO_H

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class SequenceInfo
{
    private:
        int idx;
        string desc;
        string allele_number;
        uint64_t length;
        uint64_t kmerCount;
        double normKmerCount;
        double avgKmerDepth;
        double kmerDepthSd;
        double score;
        vector<uint64_t> perBaseDepths;
        uint64_t nBasesCovered;
        // flags
        bool kmerDepthSdSet;

    protected:

    public:
        
        enum Attribute  // Enum to specify the attribute for sorting collections of SequenceInfo objects
        {
            COVERAGE,   
            DEPTH_SD,
            SCORE
        };
        
        
        SequenceInfo();
        
        SequenceInfo(const int&      idx, 
                     const string&   desc, 
                     const uint64_t& length);
        
        SequenceInfo(const int&      idx, 
                     const char*     desc, 
                     const uint64_t& length);
        
        SequenceInfo(const SequenceInfo& obj);
        
        virtual ~SequenceInfo();

        bool calculateKmerDepthSd();
        void calculateScore();
        // virtual double getNormKmerCount() const;
        virtual double getCoverage() const;
        virtual string getCoverageStr(const uint64_t& precision = 1) const;
        
        virtual bool operator<(const SequenceInfo & other);
        virtual bool operator>(const SequenceInfo & other);
        friend ostream& operator<<(ostream& os, const SequenceInfo& object);
        
        // Getters and Setters
        virtual int getIdx() const;
        virtual string getDesc() const;
        virtual uint64_t getLength() const;
        virtual string getAlleleNumber() const;
        virtual uint64_t getKmerCount() const;
        virtual double getNormKmerCount() const;
        virtual double getAvgKmerDepth() const;
        virtual double getKmerDepthSd() const;
        virtual vector<uint64_t> getPerBaseDepths() const;
        virtual uint64_t getNBasesCovered() const;
        virtual bool isKmerDepthSdSet() const;
        virtual void setIdx(const int & valIdx);
        virtual void setDesc(const string & valDesc);
        virtual void setAlleleNumber(const string & valAlleleNumber);
        virtual void setLength(const uint64_t & valLength);
        virtual void setKmerCount(const uint64_t & valKmerCount);
        virtual void setNormKmerCount(const double & valNormKmerCount);
        virtual void setAvgKmerDepth(const double & valAvgKmerDepth);
        virtual void setKmerDepthSd(const double & valKmerDepthSd);
        virtual void setScore(const double & valScore);
        virtual double getScore() const;
        virtual void setPerBaseDepths(const vector<uint64_t> & valPerBaseDepths);
        virtual void setNBasesCovered(const uint64_t & valNBasesCovered);
};

struct SSequenceInfoLess {
  SequenceInfo::Attribute attribute;
  // SSequenceInfoLess(SequenceInfo::Attribute attribute) {this->attribute = attribute;}
  SSequenceInfoLess(SequenceInfo::Attribute attribute):attribute(attribute) {}
  bool operator()(const SequenceInfo& a, const SequenceInfo& b) const {
      if (attribute == SequenceInfo::COVERAGE)
          return a.getCoverage() < b.getCoverage();
      else if (attribute == SequenceInfo::DEPTH_SD)
          return a.getKmerDepthSd() < b.getKmerDepthSd();
      else 
          return a.getScore() < b.getScore();
  }
};

struct SSequenceInfoGreater {
  SequenceInfo::Attribute attribute;
  // SSequenceInfoGreater(SequenceInfo::Attribute attribute) {this->attribute = attribute;}
  SSequenceInfoGreater(SequenceInfo::Attribute attribute):attribute(attribute) {}
  bool operator()(const SequenceInfo& a, const SequenceInfo& b) const {
      if(attribute == SequenceInfo::COVERAGE)
          return a.getCoverage() > b.getCoverage();
      else if (attribute == SequenceInfo::DEPTH_SD)
          return a.getKmerDepthSd() > b.getKmerDepthSd();
      else
        return a.getScore() > b.getScore();
  }
};

struct SSequenceInfoSort {
  // SequenceInfo::Attribute attribute;
  // SSequenceInfoGreater(SequenceInfo::Attribute attribute) {this->attribute = attribute;}
  SSequenceInfoSort(){}
  bool operator()(const SequenceInfo& a, const SequenceInfo& b) const {
      
    // return ((a.getNormKmerCount() > b.getNormKmerCount()) ||
    //         ((a.getNormKmerCount() == b.getNormKmerCount()) && a.getCoverage() > b.getCoverage()) ||
    //         ((a.getNormKmerCount() == b.getNormKmerCount()) && (a.getCoverage() == b.getCoverage()) && (a.getKmerDepthSd() > b.getKmerDepthSd()))
    //        );
    return ((a.getCoverage() > b.getCoverage()) ||
            ((a.getCoverage() == b.getCoverage()) && a.getNormKmerCount() > b.getNormKmerCount()) ||
            ((a.getCoverage() == b.getCoverage()) && (a.getNormKmerCount() > b.getNormKmerCount()) && (a.getKmerDepthSd() < b.getKmerDepthSd()))
           );
  }
};

#endif
