/*
 * histogramm.cpp
 *
 *  Created on: 20.09.2013
 *      Author: mrupp
 */

#include "histogramm.h"
#include <algorithm>    // std::sort
#include <sstream>
#include <iomanip> // setprecision


namespace ug
{

// minimum width is 8
std::string scientificStr(double d, int width)
{
	std::stringstream s;
	if(d < 0)
		std::setprecision(width-7);
	else
		s << std::setprecision(width-6);
	s << std::scientific << std::setw(width) << d;
	return s.str();
}


std::string HistogrammString(std::vector<double> values)
{
	if(values.size() == 0) return "";
	std::stringstream ss;
	std::sort(values.begin(), values.end());

	double vmin = values[0];
	double vmax = values[values.size()-1];

	const size_t N = 80;
	const double Nd=80.0;
	std::vector<size_t> hist;
	hist.resize(N, 0);

	size_t histmax=0;
	for(size_t i=0; i<values.size(); i++)
	{
		size_t idx = (size_t)((N-1)*(values[i]-vmin)/(vmax-vmin));
		hist[idx] ++;
		if(histmax < hist[idx]) histmax = hist[idx];
	}

	ss << "histmax = " << histmax << " vmin = " << vmin << " vmax = " << vmax << "\n";
	const int height = 10;
	for(int h=height; h>=0; h--)
	{
		ss << (histmax*h)/((double)height) << "\t";
		for(size_t i=0; i<N; i++)
		{
			if(hist[i]*height > histmax*h) ss << "#";
			else ss << " ";
		}
		ss <<"\n";
	}

	ss << vmin << "\t";
	for(size_t i =0; i<N; i+=9)
		ss << "    |    ";

	ss << "\n" << vmin << "\t";
	for(size_t i =0; i<N; i+=9)
		ss << scientificStr((i+4)*(vmax-vmin)/(Nd-1)+vmin, 8) << " ";
	ss << vmax << "\n";

	for(size_t i=0; i<N; i++)
		ss << "-";
	ss << "\n";
	return ss.str();
}

}
