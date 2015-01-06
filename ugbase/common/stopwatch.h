/**
 * \file common/stopwatch.h
 * \author Martin Rupp
 * \date 2010-03-03
 * \author Torbjoern Klatt
 * \date 2012-05-07
 * \copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 * \brief stopwatch class for quickly taking times
 */

#ifndef __H__UG__STOPWATCH_H__
#define __H__UG__STOPWATCH_H__

#include <iostream>

#ifdef UG_CXX11
#include <chrono>
#else
#include <ctime>
#endif

#ifdef UG_POSIX
  #include <sys/time.h>
#endif

namespace ug
{



class CuckooClock{

 public:
 CuckooClock() : m_tlaunch(0), m_tdone(0), m_ttotal(0) {}

  void tic(){m_tlaunch = clock();}

  //! returns time since last tic
  double toc()
  {
    m_tdone = clock();
    //std::cerr << m_tdone;
    clock_t delta= (m_tdone-m_tlaunch);
    m_ttotal += delta;
    return 1.0*delta/CLOCKS_PER_SEC;
  }

  //! returns total time (in seconds)
  double cuckoo() {return 1.0*m_ttotal/ CLOCKS_PER_SEC;}

protected:

    clock_t m_tlaunch;
    clock_t m_tdone;
    clock_t m_ttotal;
};


#ifdef UG_POSIX
	inline double get_clock_s()
  {
    timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + time.tv_usec/1000000.0;
  }
#else
  inline double get_clock_s()
  {
    return clock() / ((double)CLOCKS_PER_SEC);
  }
#endif

/// \addtogroup ugbase_common
/// \{

/**
 * \brief Stopwatch class for quickly taking times
 *
 * Depending on \c CXX11 flag, two different versions are compiled.
 * If \c CXX11=ON, <tt>std::chrono</tt> from C++11's STL is used providing high
 * resolution (microseconds) time measuring.
 * Otherwise <tt>std::ctime</tt> is used providing millisecond resolution.
 *
 * \note If \c CXX11=OFF timings shorter than 100ms seem to be rather
 * inaccurate.
 */
class Stopwatch
{
  public:
    /**
     * \brief Default constructor for the Stopwatch
     */
    Stopwatch() {
      // you cant be really sure when constructor is called
#ifdef UG_CXX11
      begin = std::chrono::high_resolution_clock::now();
      end = std::chrono::high_resolution_clock::now() - begin;
#else
      beg = end = get_clock_s();
#endif
      bRunning = false;
    }
    
    /**
     * \brief Starts the Stopwatch
     */
    void start() {
      std::cout.flush();
#ifdef UG_CXX11
      begin = std::chrono::high_resolution_clock::now();
#else
      beg = get_clock_s();
#endif
      bRunning = true;
    }

    /**
     * \brief Stops the Stopwatch
     */
    void stop() {
#ifdef UG_CXX11
      end = std::chrono::high_resolution_clock::now() - begin;
#else
      end = get_clock_s();
#endif
      bRunning = false;
    }

    /**
     * \brief Prints number of milliseconds since call of start() to ostream
     *
     * Pretty prints the amount of milliseconds passed between calls of
     * Stopwatch::start() and Stopwatch::stop() or this function call
     * to the specified std::ostream.
     *
     * \param[out]  out std::ostream to print number of milliseconds to
     * \param[in]   sw  a Stopwatch instance (usualy 'this')
     */
    friend std::ostream &operator << ( std::ostream &out, Stopwatch &s ) {
      out << s.ms() << " ms";
      return out;
    }

    /**
     * \brief Returns milliseconds since call of start
     *
     * Returns the amount of milliseconds passed between calls of
     * Stopwatch::start() and Stopwatch::stop() or this function call.
     *
     * \note If compiled with \c CXX11=ON returned milliseconds have microsecond
     * resolution.
     * 
     * \return milliseconds
     */
    double ms() {
#ifdef UG_CXX11
      if ( bRunning ) end = std::chrono::high_resolution_clock::now();
      std::chrono::microseconds ellapsed = duration_cast<duration<double>> (end-begin);
      return ellapsed.count() / 1000.0;
#else
      if( bRunning ) end = get_clock_s();
      return ( end - beg ) * 1000.0;
#endif
    }

  private:
#ifdef UG_CXX11
    /// Time point of the start of Stopwatch
    std::chrono::high_resolution_clock::time_point begin;
    /// Number of microseconds since \c begin

    std::chrono::high_resolution_clock::time_point end;
#else
    /// Time point of the start of Stopwatch
    double beg;
    /// Time point of the end of Stopwatch
    double end;
#endif
    /// Flag indicating state of Stopwatch
    bool bRunning;
};

// end group ugbase_common
/// \}

} // namespace ug


#endif // __H__UG__STOPWATCH_H__
