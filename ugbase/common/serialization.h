//	created by Sebastian Reiter, Martin Rupp
//	s.b.reiter@googlemail.com
//	y10 m12 d2

#ifndef __H__UG__SERIALIZATION__
#define __H__UG__SERIALIZATION__

#include <iostream>

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

}//	end of namespace

#endif
