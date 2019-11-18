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


#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <exception>
#include <string>
#include <iostream>
// #include <ctime>
// #include <ratio>

using namespace std::chrono;
using namespace std;

class Timer {
    private:
        steady_clock::time_point t1;
        steady_clock::time_point t2;
        bool started;
        bool ended;
        
    public:
        Timer();
        void start();
        void end();
        void reset();
        duration<double> getElapsedTime();
};

class NoStartException : public exception
{
        string outMessage;
    public:
        // NoStartException() : 
        //     outMessage(string("No start time defined. Call start() method before calling end() or getElapsedTime()"))
        // {}
        virtual const char* what() const throw()
        {
            // return outMessage.c_str();
            return "No start time defined. Call start() method before calling end() or getElapsedTime()";
        }
};

class NoEndException : public exception
{
    public:
        // NoStartException(){}
        virtual const char* what() const throw()
        {
            // return outMessage.c_str();
            return "No end time defined. Call end() method before calling getElapsedTime()";
        }
};

class NoStartEndException : public exception
{
        string outMessage;
    public:
        // NoStartException() : 
        //     outMessage(string("No start and end times defined. Call start() and end() methods before obtaining elapsed time"))
        // {}
        virtual const char* what() const throw()
        {
            // return outMessage.c_str();
            return "No start and end times defined. Call start() and end() methods before obtaining elapsed time";
        }
};

#endif /* TIMER_H_ */