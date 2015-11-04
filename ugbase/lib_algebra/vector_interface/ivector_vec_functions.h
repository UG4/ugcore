
#ifndef IVECTOR_VEC_FUNCTIONS_H_
#define IVECTOR_VEC_FUNCTIONS_H_


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// T, IVector functions
template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const IVector &v2)
{
	VecAdd(a1, v1, a2, v2.downcast<T>());
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const IVector &v2, double a3, const IVector &v3)
{
	VecAdd(a1, v1, a2, v2.downcast<T>(), a3, v3.downcast<T>());
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const IVector &v2, double a3, const IVector &v3, double a4, const IVector &v4)
{
	VecAdd(a1, v1, a2, v2.downcast<T>(), a3, v3.downcast<T>(), a4, v4.downcast<T>());
}

template<typename T>
inline double VecProd(T &v1, IVector &v2)
{
	return VecProd(v1, v2.downcast<T>());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IVector, IVector functions
void VecAdd(double a1, IVector &v1, double a2, const IVector &v2)
{
	v1.vec_add(a1, a2, v2);
}

void VecAdd(double a1, IVector &v1, double a2, const IVector &v2, double a3, const IVector &v3)
{
	v1.vec_add(a1, a2, v2, a3, v3);
}

void VecAdd(double a1, IVector &v1, double a2, const IVector &v2, double a3, const IVector &v3, double a4, const IVector &v4)
{
	v1.vec_add(a1, a2, v2, a3, v3, a4, v4);
}

double VecProd(const IVector &v1, const IVector &v2)
{
	return v1.vec_prod(v2);
}

double VecNorm2(const IVector &v1)
{
	return v1.norm2();
}

double VecNorm(const IVector &v1)
{
	return v1.norm();
}


#endif /* IVECTOR_VEC_FUNCTIONS_H_ */
