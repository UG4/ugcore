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
	using std::vector<T>::size;
	using std::vector<T>::operator[];
	
	T index(size_t i) { return operator[](i); }
	std_vector_wrap() : std::vector<T>() {}
	std_vector_wrap(size_t s) : std::vector<T>(s) {}
	std_vector_wrap(size_t s, T init) : std::vector<T>(s, init) {}
	std::string tostring()
	{
		std::stringstream ss;
		ss << "vector[" << size() << "] = ";
		if(size()==0) return ss.str();
		for(size_t i=0; i<size()-1; i++)
			ss << i << " = { " << index(i) << " }, ";
		ss << size()-1 << " = { " << index(size()-1) << " }";
		return ss.str();
	}	
	void set(size_t i, T t) { operator[](i) = t; }
	void resize(size_t s) { std::vector<T>::resize(s); }
	void resize(size_t s, T init) { std::vector<T>::resize(s, init); }
	void push_back(T t) { std::vector<T>::push_back(t); }
};

}
}
#endif	/* STDVECTORWRAP_H */

