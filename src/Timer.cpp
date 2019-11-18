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


#include "Timer.h"

using namespace std::chrono;
using namespace std;

Timer::Timer()
{
    // time_point<high_resolution_clock> default_value = time_point<high_resolution_clock>();
    // t1      = default_value;
    // t2      = default_value;
    started = false;
    ended   = false;
}

void Timer::start()
{
    t1 = steady_clock::now();
    started = true;
}

void Timer::end()
{
    if (!started)
        throw NoStartException();
    
    t2 = steady_clock::now();
    ended = true;
}

void Timer::reset()
{
    Timer();
    // t1 = duration_cast<duration<double>>(0);
    // t2 = duration_cast<duration<double>>(0);
    // started = false;
    // ended   = false;
}

duration<double> Timer::getElapsedTime()
{
    if (!started && !ended) {
        throw NoStartEndException();
    } 
    else if (!started){
        throw NoStartException();
    }
    else if (!ended)
        throw NoEndException();
    else {        
        return duration_cast<duration<double>>(t2 - t1);
    }
    
    duration<double> default_value(0.0); 
    
    return default_value;
}


// int main(int argc, char const **argv)
// {
//     using namespace std;
    
//     Timer tmr1;
//     Timer tmr2;
    
//     // try {
//     //     cout << "Elapsed time: " << tmr1.getElapsedTime().count() << endl;
//     // }
//     // catch(const std::exception& e) {
//     //     std::cerr << e.what() << '\n';
//     // }
    
//     // tmr1.end();
//     // try {
//     //     cout << "Elapsed time: " << tmr1.getElapsedTime().count() << endl;
//     // }
//     // catch(const std::exception& e) {
//     //     std::cerr << e.what() << '\n';
//     // }
    
//     // tmr1.start();
//     // try {
//     //     cout << "Elapsed time: " << tmr1.getElapsedTime().count() << endl;
//     // }
//     // catch(const std::exception& e) {
//     //     std::cerr << e.what() << '\n';
//     // }
    
//     tmr2.start();
//     for (int i=0; i<50000; ++i) std::cout << "*";
//     cout << endl;
//     tmr2.end();
//     cout << "Elapsed time: " << tmr2.getElapsedTime().count() << endl;
//     tmr2.reset();
    
//     tmr2.start();
//     for (int i=0; i<10000; ++i) std::cout << "*";
//     cout << endl;
//     tmr2.end();
//     cout << "Elapsed time: " << tmr2.getElapsedTime().count() << endl;
    
//     return 0;
// }
