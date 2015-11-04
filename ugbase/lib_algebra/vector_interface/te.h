
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
	typedef T vector_type;
	double a;
	const T& v;

	double scaling() const { return a;}
	const T &vec() const { return v; }

	inline TE_VecScale(double a_, const T & v_) : a(a_), v(v_)
	{  }

};

template<class T>
class TE_Vector : public TE_AlphaVec<TE_Vector<T>  >
{
public:
	typedef T vector_type;
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
