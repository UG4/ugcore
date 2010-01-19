#pragma once

template <typename t> class matrix_trait;
template <typename t> class vec_traits;
template<typename entry_type, typename vec_type> struct Mult_Traits;
#include "blockDenseMatrix.h"
#include "blockVector.h"



/*
template<int n>
inline double mnorm2(const  blockVector<n> &v)
{
	return v.norm2();
}*/

template<typename TYPE>
inline double mnorm2(const TYPE &v)
{
	return v.norm2();
}

template<typename TYPE>
inline double mnorm(const TYPE &v)
{
	return sqrt(v.norm2());
}

template <>
inline double mnorm(const double &a)
{
	return dabs(a);
}

template <>
inline double mnorm2(const double &a)
{
	return a*a;
}

//////////////////////////////////////////////////////

template<typename M> inline double getAt(const M &m, int i)
{
	return m(i);
}
template<> inline double getAt(const double &m, int i) { return m; }

//////////////////////////////////////////////////////

template<typename M> inline double getAt(const M &m, int i, int j)
{
	return m(i, j);
}
template<> inline double getAt(const double &m, int i, int j) { return m; }

//////////////////////////////////////////////////////

template<typename M> inline double setAt(M &m, int i, double a)
{
	return m(i) = a;
}
template<> inline double setAt(double &m, int i, double a) { m = a; return a; }

//////////////////////////////////////////////////////

template<typename M> inline double setAt(M &m, int i, int j, double a)
{
	return m(i, j) = a;
}
template<> inline double setAt(double &m, int i, int j, double a) { m = a; return a; }




//////////////////////////////////////////////////////
// wrapper for using doubles
//setSize(t, a, b)
template<typename T>
inline void setSize(T &t, int a, int b)
{
	t.setSize(a, b);
}

template<>
inline void setSize(double &d, int a, int b)
{
	return;
}
//////////////////////////////////////////////////////
//setSize(t, a)
template<typename T>
inline void setSize(T &t, int a)
{
	t.setSize(a);
}

template<>
inline void setSize(double &d, int a)
{
	return;
}
//////////////////////////////////////////////////////
// getSize
template<typename T>
inline int getSize(T &t)
{
	return t.getSize();
}

template<>
inline int getSize(double &t)
{
	return 1;
}

//////////////////////////////////////////////////////
//getRows
template<typename T>
inline int getRows(const T &t)
{
	return t.getRows();
}

template<typename T>
inline int getCols(const T &t)
{
	return t.getCols();
}

//////////////////////////////////////////////////////
template<>
inline int getRows(const double &t)
{
	return 1;
}

template<>
inline int getCols(const double &t)
{
	return 1;
}


///////////////////////////////////////////////////////////////////




template<>
struct matrix_trait<double>
{
	typedef double vec_type;
	typedef double inverse_type;
	enum { nrOfUnknowns = 1 } ;
};


////////////////////////////////////////////////////////////////////


template<>
struct vec_traits<double>
{
	enum { nrOfUnknowns = 1 };
};



////////////////////////////////////////////////////////////////////

template<> struct Mult_Traits<double, double>
{
	typedef double ReturnType;
};





