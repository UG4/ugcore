#pragma once

template <typename t> class matrix_trait;
template <typename t> class vec_traits;
template<typename mat_type, typename vec_type> struct Mult_Traits;
#include "fixedMatrix.h"


/*
template<int n>
inline double mnorm2(const  fixedVector<n> &v)
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

template<typename M> inline double getAt(const M &m, int i)
{
	return m(i);
}
template<> inline double getAt(const double &m, int i) { return m; }

template<typename M> inline double getAt(const M &m, int i, int j)
{
	return m(i, j);
}
template<> inline double getAt(const double &m, int i, int j) { return m; }


template<typename M> inline double setAt(M &m, int i, double a)
{
	return m(i) = a;
}
template<> inline double setAt(double &m, int i, double a) { m = a; return a; }

template<typename M> inline double setAt(M &m, int i, int j, double a)
{
	return m(i, j) = a;
}
template<> inline double setAt(double &m, int i, int j, double a) { m = a; return a; }


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





