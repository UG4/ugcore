/*
 *  stopwatch.h
 *
 *  Created by Martin Rupp on 03.03.10.
 *  C++11 high resolution added by Torbjörn Klatt on 2012-05-07
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__LIB_ALGEBRA__STOPWATCH_H__
#define __H__UG__LIB_ALGEBRA__STOPWATCH_H__

#ifdef UG_CXX11
#include <chrono>
#else
#include <ctime>
#endif

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////////////
//!
//! stopwatch class for quickly taking times
//! seems to be ok for measuring times > 100 ms
class stopwatch
{
  public:
    stopwatch() {
      // you cant be really sure when constructor is called
#ifdef UG_CXX11
      begin = chrono::high_resolution_clock::now();
      end = chrono::high_resolution_clock::now() - begin;
#else
      beg = end = clock();
#endif
      bRunning = false;
    }
    void start() {
      cout.flush();
#ifdef UG_CXX11
      begin = chrono::high_resolution_clock::now();
#else
      beg = clock();
#endif
      bRunning = true;
    }
    void stop() {
#ifdef UG_CXX11
      end = chrono::high_resolution_clock::now() - begin;
#else
      end = clock();
#endif
      bRunning = false;
    }
    friend std::ostream &operator << ( std::ostream &out, stopwatch &s ) {
      out << s.ms() << " ms";
      return out;
    }

    double ms() {
#ifdef UG_CXX11
      if ( bRunning ) end = chrono::high_resolution_clock::now() - begin;
      return end.count() / 100.0;
#else
      if( bRunning ) end = clock();
      return ( end - beg ) / ( ( double )0.001 * CLOCKS_PER_SEC );
#endif
    }

  private:
#ifdef UG_CXX11
    chrono::high_resolution_clock::time_point begin;
    chrono::microseconds end;
#else
    clock_t beg, end;
#endif
    bool bRunning;
};


} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__STOPWATCH_H__
