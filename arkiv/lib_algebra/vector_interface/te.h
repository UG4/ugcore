/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#ifndef __H__UG__LIB_ALGEBRA__TEMPLATE_EXPRESSIONS2__
#define __H__UG__LIB_ALGEBRA__TEMPLATE_EXPRESSIONS2__

//#include "blockMatrix.h"


namespace ug{

////////////////////////////////////////////////////////////////////////////////
//! AlphaVec_Expression
//! class for Template Expressions of the form
//! double * Vector.
//! \attention to case x = alpha1*y + alpha1*x
template<typename T>
class TE_AlphaVec
{
public:
	const T& cast() const {return static_cast<const T&>(*this); }
	double scaling() const { return cast().scaling(); }
};

template<typename T>
class TE_VecScale : public TE_AlphaVec<TE_VecScale<T>  >
{
public:
	using vector_type = T ;
	double a;
	const T& v;

	double scaling() const { return a;}
	const T &vec() const { return v; }

	inline TE_VecScale(double a_, const T & v_) : a(a_), v(v_)
	{  }

};

template<typename T>
class TE_Vector : public TE_AlphaVec<TE_Vector<T>  >
{
public:
	using vector_type = T ;
	double scaling() const { return 1.0;}
	const T& vec() const {return static_cast<const T&>(*this); }
};



////////////////////////
// alpha*v

template<typename T>
TE_VecScale <typename T::vector_type>  operator * (double alpha, const TE_AlphaVec<T> &r)
{
	return TE_VecScale<typename T::vector_type> (r.scaling()*alpha, r.cast().vec());
}

template<typename T>
TE_AlphaVec <typename T::vector_type>  operator * (const TE_AlphaVec<T> &l, double alpha)
{
	return TE_VecScale<typename T::vector_type> (l.scaling()*alpha, l.cast().vec());
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! AlphaMat_Expression
//! use this class to do double * TE_MAT


template<typename T>
class TE_VecAdd2
{
public:
	double a1;
	const T& v1;
	double a2;
	const T& v2;
	inline TE_VecAdd2(double a1_, const T& v1_, double a2_, const T& v2_)
		: a1(a1_), v1(v1_), a2(a2_), v2(v2_) {}
};


template<typename T>
class TE_VecAdd3
{
public:
	double a1;
	const T& v1;
	double a2;
	const T& v2;
	double a3;
	const T& v3;
	inline TE_VecAdd3(double a1_, const T& v1_, double a2_, const T& v2_, double a3_, const T& v3_)
		: a1(a1_), v1(v1_), a2(a2_), v2(v2_), a3(a3_), v3(v3_) {}
};




// ////////////////////////
// v1 + v2

//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename L, typename R>
TE_VecAdd2<typename L::vector_type> operator + (const TE_AlphaVec<L> &l, const TE_AlphaVec<R> &r)
{
	return TE_VecAdd2<typename L::vector_type> (l.scaling(), l.cast().vec(), r.scaling(), r.cast().vec());
}
template<typename L, typename R>
TE_VecAdd2<typename L::vector_type> operator - (const TE_AlphaVec<L> &l, const TE_AlphaVec<R> &r)
{
	return TE_VecAdd2<typename L::vector_type> (l.scaling(), l.cast().vec(), -r.scaling(), r.cast().vec());
}

////////////////////////
// (v1+v2)+v3
//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename L, typename R>
TE_VecAdd3<typename L::vector_type> operator + (const TE_VecAdd2<L> &l, const TE_AlphaVec<R> &r)
{
	return TE_VecAdd3<typename L::vector_type> (l.a1, l.v1, l.a2, l.v2, r.scaling(), r.cast().vec());
}
template<typename L, typename R>
TE_VecAdd3<typename L::vector_type> operator + (const TE_AlphaVec<R> &r, const TE_VecAdd2<L> &l)
{
	return TE_VecAdd3<typename L::vector_type> (r.scaling(), r.cast().vec(), l.a1, l.v1, l.a2, l.v2);
}

template<typename L, typename R>
TE_VecAdd3<typename L::vector_type> operator - (const TE_VecAdd2<L> &l, const TE_AlphaVec<R> &r)
{
	return TE_VecAdd3<typename L::vector_type> (l.a1, l.v1, l.a2, l.v2, -r.scaling(), r.cast().vec());
}
template<typename L, typename R>
TE_VecAdd3<typename L::vector_type> operator - (const TE_AlphaVec<R> &r, const TE_VecAdd2<L> &l)
{
	return TE_VecAdd3<typename L::vector_type> (r.scaling(), r.cast().vec(), -l.a1, l.v1, -l.a2, l.v2);
}

////////////////////////
// alpha*(v1+v2)

//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename T>
TE_VecAdd2<T> operator * (const TE_VecAdd2<T> &t, double alpha)
{
	return TE_VecAdd2<T> (t.a1*alpha, t.v1, t.a2*alpha, t.v2);
}

//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename T>
TE_VecAdd2<T> operator * (double alpha, const TE_VecAdd2<T> &t)
{
	return TE_VecAdd2<T> (t.a1*alpha, t.v1, t.a2*alpha, t.v2);
}


//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename T>
TE_VecAdd3<T> operator * (const TE_VecAdd3<T> &t, double alpha)
{
	return TE_VecAdd2<T> (t.a1*alpha, t.v1, t.a2*alpha, t.v2, t.a3*alpha, t.v3);
}

//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename T>
TE_VecAdd3<T> operator * (double alpha, const TE_VecAdd2<T> &t)
{
	return TE_VecAdd3<T> (t.a1*alpha, t.v1, t.a2*alpha, t.v2, t.a3*alpha, t.v3);
}

}


#endif
