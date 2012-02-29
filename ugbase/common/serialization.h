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
#include "util/variant.h"
#include "error.h"
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

/// method returning value directly
template<typename T, class TIStream>
T Deserialize(TIStream &stream)
{
	T t;
	Deserialize(stream, t);
	return t;
}


////////////////////////////////////////////////////////////////////////////////
//	All specializations should be pre-declared here!
//	This is important, so that the different methods know of each other.
template <class T, class TOStream>
void Serialize(TOStream& buf, const std::set<T>& m);

template <class T, class TIStream>
void Deserialize(TIStream& buf, std::set<T>& myset);

template <class TOStream>
void Serialize(TOStream& buf, const std::string& str);

template <class TIStream>
void Deserialize(TIStream& buf, std::string& str);

template <class TOStream>
void Serialize(TOStream& buf, const Variant& v);

template <class TIStream>
void Deserialize(TIStream& buf, Variant& v);

template <class T, class TOStream>
void Serialize(TOStream& buf, const std::vector<T>& vec);

template <class T, class TIStream>
void Deserialize(TIStream& buf, std::vector<T>& vec);

template <class TOStream>
inline void Serialize(TOStream& buf, const std::vector<bool>::reference& boolRef);

template <class TIStream>
inline void Deserialize(TIStream& buf, std::vector<bool>::reference boolRef);

template <class Key, class T, class TOStream>
void Serialize(TOStream& buf, const std::map<Key, T>& m);

template <class Key, class T, class TIStream>
void Deserialize(TIStream& buf, std::map<Key, T>& m);
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
///	writes data from a set to a binary stream
/**	This template method is used in ug when it comes to writing data
 * from a map into a binary stream.
 * In its default implementation, it first writes the size of the map
 * and then serializes the entries.
 */
template <class T, class TOStream>
void Serialize(TOStream& buf, const std::set<T>& m)
{
	Serialize<size_t>(buf, m.size());
	for(typename std::set<T>::const_iterator it = m.begin(); it != m.end(); ++it)
		Serialize<T>(buf, *it);
}

///	deserializes data from a binary stream into a set
template <class T, class TIStream>
void Deserialize(TIStream& buf, std::set<T>& myset)
{
	myset.clear();
	size_t size = Deserialize<size_t>(buf);
	T t;
	for(size_t i = 0; i < size; ++i)
	{
		Deserialize<T>(buf, t);
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
	Serialize<size_t>(buf, len);
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

	size_t len = Deserialize<size_t>(buf);

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


///	serializes a variant
/**	Note that pointers can't be serialized in a meaningful way. We thus simply
 * do not serialize them. During deserialization pointers will be set to NULL.
 *
 * Note that c-strings (const char*) are converted to std::strings before
 * serialization. This means that they will be deserialized as std::string
 */
template <class TOStream>
void Serialize(TOStream& buf, const Variant& v)
{
	Serialize<uint32, TOStream>(buf, uint32(v.type()));
	switch(v.type()){
		case Variant::VT_INVALID:	break;
		case Variant::VT_BOOL:		Serialize(buf, v.to_bool()); break;
		case Variant::VT_INT:		Serialize(buf, v.to_int()); break;
		case Variant::VT_FLOAT:		Serialize(buf, v.to_float()); break;
		case Variant::VT_DOUBLE:	Serialize(buf, v.to_double()); break;
		case Variant::VT_CSTRING:	Serialize(buf, std::string(v.to_c_string())); break;
		case Variant::VT_STDSTRING:	Serialize(buf, v.to_std_string()); break;
		case Variant::VT_POINTER:	break;

		default:
			UG_THROW("Unknown variant type in Serialize:" << v.type());
			break;
	}
}

///	deserializes data from a binary stream into a variant
/** Note that pointers can't be serialized in a meaningful way. We thus simply
 * do not serialize them. During deserialization pointers will be set to NULL.
 */
template <class TIStream>
void Deserialize(TIStream& buf, Variant& v)
{
	uint32 type = Deserialize<uint32, TIStream>(buf);
	switch(type){
		case Variant::VT_INVALID:	v = Variant(); break;
		case Variant::VT_BOOL:		v = Variant(Deserialize<bool>(buf)); break;
		case Variant::VT_INT:		v = Variant(Deserialize<int>(buf)); break;
		case Variant::VT_FLOAT:		v = Variant(Deserialize<float>(buf)); break;
		case Variant::VT_DOUBLE:	v = Variant(Deserialize<double>(buf)); break;
		case Variant::VT_CSTRING:	// fallthrough
		case Variant::VT_STDSTRING: v = Variant(Deserialize<std::string>(buf)); break;

		case Variant::VT_POINTER:{
			void* val = NULL;
			v = Variant(val);
		} break;

		default:
			UG_THROW("Unknown variant type in Deserialize: " << type);
			break;
	}
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
	Serialize<size_t>(buf, size);
	for(size_t i = 0; i < size; ++i){
		Serialize(buf, vec[i]);
	}
}

///	deserializes data from a binary stream into a vector
template <class T, class TIStream>
void Deserialize(TIStream& buf, std::vector<T>& vec)
{
	vec.clear();
	size_t size = Deserialize<size_t>(buf);
	vec.resize(size);
	for(size_t i = 0; i < size; ++i){
		Deserialize(buf, vec[i]);
	}
}


template<class TIStream>
void Serialize(TIStream &buf, const std::vector<bool> &vec)
{
	size_t size=vec.size();
	Serialize<size_t>(buf, size);
	int j=0;
	char a=0;
	for(size_t i = 0; i < size; ++i)
	{
		if(vec[i]) a |= (1 << j);
		if(++j == 8)
		{
			Serialize<char>(buf, a);
			a = 0;
			j = 0;
		}
	}
	if(j) Serialize<char>(buf, a);
}

template<class TIStream>
void Deserialize(TIStream &buf, std::vector<bool> &vec)
{
	vec.clear();
	size_t size = Deserialize<size_t>(buf);
	vec.resize(size);
	int j=8;
	char a=0;
	for(size_t i = 0; i < size; ++i, ++j)
	{
		if(j==8)
		{
			Deserialize<char>(buf, a);
			j=0;
		}
		vec[i] = a & (1 << j);
		j++;
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
	Serialize<size_t>(buf, m.size());

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
	size_t size = Deserialize<size_t>(buf);
	for(size_t i = 0; i < size; ++i)
	{
		Key k;
		Deserialize(buf, k);
		Deserialize(buf, m[k]);
	}
}

}//	end of namespace


#endif
