
#include "histogramm.h"
#include <algorithm>    // std::sort
#include <sstream>
#include <iomanip> // setprecision


namespace ug
{

// minimum width is 8
std::string scientificStr(double d, size_t width)
{
	std::stringstream s;
	if(d < 0)
		std::setprecision(width-7);
	else
		s << std::setprecision(width-6);
	s << std::left << std::scientific << std::setw(width) << d;
	return s.str();
}

std::string cutString(double d, size_t width)
{
	std::stringstream s;
	s << std::left << std::setw(width) << d;
	if(s.str().length() == width)
		return s.str();
	return scientificStr(d, width);
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


std::string DistributionPercentage(std::vector<double> values)
{
	if(values.size() == 0) return "";
	std::stringstream ss;
	std::sort(values.begin(), values.end());

	double vmin = values[0];
	double vmax = values[values.size()-1];
	double vdiff = vmax-vmin;

	if(vdiff < 1e-12) return "all equal";

	const size_t M = 90;
	size_t N = values.size();

	const int height = 10;
	for(int h=height; h>=0; h--)
	{
		double hv = (vdiff*h)/((double)height);
		ss << cutString(hv+vmin, 8) << " ";
		for(size_t i=0; i<M; i++)
		{
			if(values[(i*(N-1))/M]-vmin >= hv) ss << ".";
			else ss << "#";

		}
		ss <<"\n";
	}

	ss << "         ";
	for(size_t i =0; i<M; i+=9)
		ss << "|        ";

	ss << "\n%smaller ";
	for(size_t i =0; i<M; i+=9)
		ss << std::setw(9) << std::left << i*100.0/M;

	ss << "\n#smaller ";
	for(size_t i =0; i<M; i+=9)
		ss << std::setw(9) << std::left << (i*N)/M;


	ss << "\n%bigger  ";
	for(size_t i =0; i<M; i+=9)
		ss << std::setw(9) << std::left << 100-i*100.0/M;
	ss << "\n#bigger  ";
	for(size_t i =0; i<M; i+=9)
		ss << std::setw(9) << std::left << N-(i*N)/M;


	ss << "\nvalue    ";
	for(size_t i =0; i<M; i+=9)
		ss << cutString(values[(i*(N-1))/M], 8) << " ";

	ss << "\n% of max ";
	for(size_t i =0; i<M; i+=9)
		ss << std::left << cutString(values[(i*(N-1))/M]*100.0/vmax, 8) << " ";


	ss << "\n---------";
	for(size_t i=0; i<M; i++)
		ss << "-";
	ss << "\n";
	return ss.str();
}

}
