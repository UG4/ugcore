/*
 * process.h
 *
 *  Created on: 29.04.2013
 *      Author: mrupp
 */

#ifndef PROCESS_H_
#define PROCESS_H_

/*
 *  progress.h
 *  amgframework
 *
 *  Created by Martin Rupp on 05.10.12.
 *  Copyright 2012 . All rights reserved.
 *
 */

#include <iostream>
#include "common/util/string_util.h"

namespace ug {
class progress
{
public:
	progress()
	{
		m_length=70;
	}
	void set_length(int l)
	{
		m_length = l;
	}
	inline void start(double total)
	{
		m_total = total;
		m_now = 0;
		std::cout << "\n." << repeat('_', m_length) << ".\n";
		std::cout << "[";
	}
	inline void set(double now)
	{
		int i=(int)(m_length*m_now/m_total);
		int i2=(int)(m_length*now/m_total);
		for(; i<i2; i++) { std::cout << "-"; std::cout.flush(); }
		m_now = now;
	}
	inline void stop()
	{
		int i=(int)(10*m_now/m_total);
		for(; i<10; i++) std::cout << "-";
		std::cout <<"]\n";
	}
private:
	double m_total;
	double m_now;
	int m_length;
};

}
#endif /* PROCESS_H_ */
