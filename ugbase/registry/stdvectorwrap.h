/* 
 * File:   stdvectorwrap.h
 * Author: mrupp
 *
 * Created on 30. Oktober 2012, 14:41
 */

#ifndef STDVECTORWRAP_H
#define	STDVECTORWRAP_H

#include <sstream>
#include <vector>


namespace ug
{
namespace bridge
{

template<typename T>
class std_vector_wrap : public std::vector<T>
{
	public:
	T index(size_t i) { return std::vector<T>::operator[](i); }
	std::string tostring()
	{
		std::stringstream ss;
		ss << "vector[" << size() << "] = ";
		if(size()==0) ss.str();
		for(size_t i=0; i<size()-1; i++)
			ss << "{ " << index(i) << " }, ";
		ss << "{ " << index(size()-1) << " }";
		return ss.str();
	}
	size_t size() { return std::vector<T>::size(); }	
};

}
}
#endif	/* STDVECTORWRAP_H */

