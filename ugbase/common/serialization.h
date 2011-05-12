//	created by Sebastian Reiter, Martin Rupp
//	s.b.reiter@googlemail.com
//	y10 m12 d2

#ifndef __H__UG__SERIALIZATION__
#define __H__UG__SERIALIZATION__

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cassert>
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
template <class T, class TOStream>
void Serialize(TOStream& buf, const T& val)
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
template <class T, class TIStream>
void Deserialize(TIStream& buf, T& valOut)
{
	assert(!buf.eof() && "End of buf reached.");
	buf.read((char*)&valOut, sizeof(T));
}

///	writes data in a vector to a binary stream
/**	This template method is used in ug when it comes to writing data
 * from a vector into a binary stream.
 * In its default implementation, it first writes the size of the vector
 * and then serializes the entries.
 */
template <class T, class TOStream>
void Serialize(TOStream& buf, const std::vector<T>& vec)
{
	size_t size = vec.size();
	buf.write((char*)&size, sizeof(size_t));
	for(size_t i = 0; i < size; ++i){
		Serialize(buf, vec[i]);
	}
}

///	deserializes data from a binary stream into a vector
template <class T, class TIStream>
void Deserialize(TIStream& buf, std::vector<T>& vec)
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
template <class TOStream>
inline void Serialize(TOStream& buf, const std::vector<bool>::reference& boolRef)
{
	char b = ((bool)boolRef) ? 1 : 0;
	buf.write(&b, sizeof(char));
}

///	deserializes data from a binary stream into a vector<bool>
// * This function is to avoid surprises with vector<bool>
// note: boolRef is not &boolRef.
template <class TIStream>
inline void Deserialize(TIStream& buf, std::vector<bool>::reference boolRef)
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
template <class Key, class T, class TOStream>
void Serialize(TOStream& buf, const std::map<Key, T>& m)
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
template <class Key, class T, class TIStream>
void Deserialize(TIStream& buf, std::map<Key, T>& m)
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
template <class T, class TOStream>
void Serialize(TOStream& buf, const std::set<T>& m)
{
	size_t size = m.size();
	buf.write((char*)&size, sizeof(size_t));
	for(typename std::set<T>::const_iterator it = m.begin(); it != m.end(); ++it)
		Serialize(buf, *it);
}

///	deserializes data from a binary stream into a set
template <class T, class TIStream>
void Deserialize(TIStream& buf, std::set<T>& myset)
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

///	Writes a string to a binary stream
/**	First the length of the string is written, then its content.*/
template <class TOStream>
void Serialize(TOStream& buf, const std::string& str)
{
	size_t len = str.length();

	buf.write((char*)&len, sizeof(size_t));
	if(len > 0)
		buf.write(str.c_str(), sizeof(char) * len);
}

///	deserializes data from a binary stream into a string
template <class TIStream>
void Deserialize(TIStream& buf, std::string& str)
{
//	the buffers allow us to read small strings fast.
//	for bigger ones we have to temporarily reserve memory.
	char staticBuf[64];
	char* flexBuf = NULL;
	char* tBuf = staticBuf;

	size_t len = 0;
	buf.read((char*)&len, sizeof(size_t));

//	check whether we have to allocate memory
//	don't forget that we have to append a zero at the end
	if(len >= 63){
		flexBuf = new char[len + 1];
		tBuf = flexBuf;
	}

	if(len > 0)
		buf.read(tBuf, sizeof(char) * len);
	tBuf[len] = 0;

//	assign data to the out-string
	str = tBuf;

//	clean up
	if(flexBuf)
		delete[] flexBuf;
}

}//	end of namespace


#endif
