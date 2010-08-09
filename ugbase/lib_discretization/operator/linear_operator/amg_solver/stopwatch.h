/*
 *  stopwatch.h
 *
 *  Created by Martin Rupp on 03.03.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__LIB_DISCRETIZATION__STOPWATCH_H__
#define __H__LIB_DISCRETIZATION__STOPWATCH_H__

namespace ug{
////////////////////////////////////////////////////////////////////////////////
//!
//! stopwatch class for quickly taking times
//! seems to be ok for measuring times > 100 ms
class stopwatch
{
public:
	stopwatch() 
	{
		// you cant be really sure when constructor is called
		beg = clock();
		bRunning = false;
	}
	void start()
	{
		cout.flush();
		beg = clock();
		bRunning = true;
	}
	void stop()
	{
		end = clock();
		bRunning = false;
	}
	friend ostream &operator << (ostream &out, stopwatch &s)
	{
		out << s.getTimeDiffMS() << " ms";
		return out;
	}
	
	double getTimeDiffMS()
	{
		if(bRunning) end = clock();
		return (clock()-beg)/((double)0.001*CLOCKS_PER_SEC);
	}
	
	void printTimeDiff()
	{
		cout << "took " << getTimeDiffMS() << " ms" << endl;
		cout.flush();
	}
private:
	clock_t beg, end;
	bool bRunning;
};


} // namespace ug


#endif // __H__LIB_DISCRETIZATION__STOPWATCH_H__
