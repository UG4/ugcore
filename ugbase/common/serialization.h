//	created by Sebastian Reiter, Martin Rupp
//	s.b.reiter@googlemail.com
//	y10 m12 d2

#ifndef __H__UG__SERIALIZATION__
#define __H__UG__SERIALIZATION__

#include <iostream>
#include <vector>
#include "log.h"
namespace ug
{

///	writes data to a binary stream
/**	This template method is used in ug when it comes to writing data to
 * a binary stream.
 * In its default implementation, it simply copies the data to the stream.
 * If you want to overload this method for your types, you can simply do
 * so by using template specialization.
 */
template <class T>
void Serialize(std::ostream& buf, const T& val)
{
	buf.write((char*)&val, sizeof(T));
}

///	reads data from a binary stream
/**	This template method is used in ug when it comes to reading data from
 * a binary stream.
 * In its default implementation, it simply copies the data from the stream.
 * If you want to overload this method for your types, you can simply do
 * so by using template specialization.
 */
template <class T>
void Deserialize(std::istream& buf, T& valOut)
{
	buf.read((char*)&valOut, sizeof(T));
}

///	writes data in a vector to a binary stream
/**	This template method is used in ug when it comes to writing data
 * from a vector into a binary stream.
 * In its default implementation, it first writes the size of the vector
 * and then serializes the entries.
 */
template <class T>
void Serialize(std::ostream& buf, const std::vector<T>& vec)
{
	size_t size = vec.size();
	buf.write((char*)&size, sizeof(size_t));
	for(size_t i = 0; i < size; ++i){
		Serialize(buf, vec[i]);
	}
}

///	writes data in a vector to a binary stream
/**	This template method is used in ug when it comes to writing data
 * from a vector into a binary stream.
 * In its default implementation, it first writes the size of the vector
 * and then serializes the entries.
 */
template <class T>
void Deserialize(std::istream& buf, std::vector<T>& vec)
{
	size_t size = 0;
	buf.read((char*)&size, sizeof(size_t));
	vec.resize(size);
	for(size_t i = 0; i < size; ++i){
		Deserialize(buf, vec[i]);
	}
}

}//	end of namespace

#endif
