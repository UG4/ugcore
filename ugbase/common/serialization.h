//	created by Sebastian Reiter, Martin Rupp
//	s.b.reiter@googlemail.com
//	y10 m12 d2

#ifndef __H__UG__SERIALIZATION__
#define __H__UG__SERIALIZATION__

#include <iostream>
#include <vector>
#include <set>
#include <map>
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

///	deserializes data from a binary stream into a vector
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

/**	This method is used in ug when it comes to writing data
 * from a vector<bool> into a binary stream.
 * This function is to avoid surprises with vector<bool>
 * because of vector<bool>::reference.
 * Note: You could also define Serialize(., vector<bool>::reference)
 */
inline void Serialize(std::ostream& buf, const std::vector<bool>::reference& boolRef)
{
	char b = ((bool)boolRef) ? 1 : 0;
	buf.write(&b, sizeof(char));
}

///	deserializes data from a binary stream into a vector<bool>
// * This function is to avoid surprises with vector<bool>
// note: boolRef is not &boolRef.
inline void Deserialize(std::istream& buf, std::vector<bool>::reference boolRef)
{
	char b;
	buf.read(&b, sizeof(char));
	boolRef = (bool)(b == 1);
}


///	writes data from a map to a binary stream
/**	This template method is used in ug when it comes to writing data
 * from a map into a binary stream.
 * In its default implementation, it first writes the size of the map
 * and then serializes the entries.
 */
template <class Key, class T>
void Serialize(std::ostream& buf, const std::map<Key, T>& m)
{
	size_t size = m.size();
	buf.write((char*)&size, sizeof(size_t));
	for(typename std::map<Key, T>::const_iterator it = m.begin(); it != m.end(); ++it)
	{
		Serialize(buf, it->first);
		Serialize(buf, it->second);
	}
}

///	deserializes data from a binary stream into a map
template <class Key, class T>
void Deserialize(std::istream& buf, std::map<Key, T>& m)
{
	m.clear();
	size_t size = 0;
	buf.read((char*)&size, sizeof(size_t));
	for(size_t i = 0; i < size; ++i)
	{
		Key k;
		Deserialize(buf, k);
		Deserialize(buf, m[k]);
	}
}

///	writes data from a set to a binary stream
/**	This template method is used in ug when it comes to writing data
 * from a map into a binary stream.
 * In its default implementation, it first writes the size of the map
 * and then serializes the entries.
 */
template <class T>
void Serialize(std::ostream& buf, const std::set<T>& m)
{
	size_t size = m.size();
	buf.write((char*)&size, sizeof(size_t));
	for(typename std::set<T>::const_iterator it = m.begin(); it != m.end(); ++it)
		Serialize(buf, *it);
}

///	deserializes data from a binary stream into a set
template <class T>
void Deserialize(std::istream& buf, std::set<T>& myset)
{
	myset.clear();

	size_t size = 0;
	buf.read((char*)&size, sizeof(size_t));
	for(size_t i = 0; i < size; ++i)
	{
		T t;
		Deserialize(buf, t);
		// using myset.end() because data t is sorted.
		myset.insert (myset.end(), t);
	}
}
	
}//	end of namespace

#endif
