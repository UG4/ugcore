/**
 * \file common/stopwatch.h
 * \author Martin Rupp
 * \date 2010-03-03
 * \author Torbjörn Klatt
 * \date 2012-05-07
 * \copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 * \brief stopwatch class for quickly taking times
 */

#ifndef __H__UG__STOPWATCH_H__
#define __H__UG__STOPWATCH_H__

#ifdef UG_CXX11
#include <chrono>
#else
#include <ctime>
#endif

using namespace std;

namespace ug
{

/**
 * \brief stopwatch class for quickly taking times
 *
 * Depending on \c CXX11 flag, two different versions are compiled.
 * If \c CXX11=ON, <tt>std::chrono</tt> from C++11's STL is used providing high
 * resolution (microseconds) time measuring.
 * Otherwise <tt>std::ctime</tt> is used providing millisecond resolution.
 *
 * \note If \c CXX11=OFF timings shorter than 100ms seem to be rather
 * inaccurate.
 */
class stopwatch
{
  public:
    /**
     * \brief Default constructor for the stopwatch
     */
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
    
    /**
     * \brief Starts the stopwatch
     */
    void start() {
      cout.flush();
#ifdef UG_CXX11
      begin = chrono::high_resolution_clock::now();
#else
      beg = clock();
#endif
      bRunning = true;
    }

    /**
     * \brief Stops the stopwatch
     */
    void stop() {
#ifdef UG_CXX11
      end = chrono::high_resolution_clock::now() - begin;
#else
      end = clock();
#endif
      bRunning = false;
    }

    /**
     * \brief Prints number of milliseconds since call of start() to ostream
     *
     * Pretty prints the amount of milliseconds passed between calls of
     * stopwatch::start() and stopwatch::stop() or this function call
     * to the specified std::ostream.
     *
     * \param[out]  out std::ostream to print number of milliseconds to
     * \param[in]   sw  a stopwatch instance (usualy 'this')
     */
    friend std::ostream &operator << ( std::ostream &out, stopwatch &s ) {
      out << s.ms() << " ms";
      return out;
    }

    /**
     * \brief Returns milliseconds since call of start
     *
     * Returns the amount of milliseconds passed between calls of
     * stopwatch::start() and stopwatch::stop() or this function call.
     *
     * \note If compiled with \c CXX11=ON returned milliseconds have microsecond
     * resolution.
     * 
     * \return milliseconds
     */
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
    /// Time point of the start of stopwatch
    chrono::high_resolution_clock::time_point begin;
    /// Number of microseconds since \c begin
    chrono::microseconds end;
#else
    /// Time point of the start of stopwatch
    clock_t beg;
    /// Time point of the end of stopwatch
    clock_t end;
#endif
    /// Flag indicating state of stopwatch
    bool bRunning;
};


} // namespace ug


#endif // __H__UG__STOPWATCH_H__
