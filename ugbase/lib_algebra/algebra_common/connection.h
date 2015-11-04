/*
 * connection.h
 *
 *  Created on: 10.05.2013
 *      Author: mrupp
 */

#ifndef CONNECTION_H_
#define CONNECTION_H_

namespace ug{

template<typename T>
class AlgebraicConnection
{
public:
	size_t iIndex;		// index to
	T dValue; // smallmatrix value;

	AlgebraicConnection() {}
	AlgebraicConnection(size_t i, const T &v)
	: iIndex(i), dValue(v) {}

	void print(){std::cout << *this;}
	friend std::ostream &operator<<(std::ostream &output, const AlgebraicConnection&c)
	{
		output << "(" << c.iIndex << "-> ";
		output << c.dValue;
		output << ")";
		return output;
	}

	void operator = (const AlgebraicConnection &other)
	{
		iIndex = other.iIndex;
		dValue = other.dValue;
	}

	int operator < (const AlgebraicConnection &c) const
	{
		return iIndex < c.iIndex;
	}

	size_t &index()
	{
		return iIndex;
	}
	T &value()
	{
		return dValue;
	}
};

}
#endif /* CONNECTION_H_ */
