
#ifndef __H__UG__MARTIN_ALGEBRA__BLOCKS__
#define __H__UG__MARTIN_ALGEBRA__BLOCKS__

namespace ug{
	
inline double dabs(double a) { return a > 0 ? a : -a; }


template <typename t> struct block_matrix_traits;
template <typename t> struct block_vector_traits;
template<typename entry_type, typename vec_type> struct block_multiply_traits;

template<typename M>
inline void GetInverse(typename block_matrix_traits<M>::inverse_type &inv, const M &m);

//////////////////////////////////////////////////////

template<typename TYPE>
inline double BlockNorm2(const TYPE &v)
{
	return v.norm2();
}

template<typename TYPE>
inline double BlockNorm(const TYPE &v)
{
	return sqrt(BlockNorm2(v));
}


//////////////////////////////////////////////////////

// get/set vector
template<typename M> inline double &BlockRef(M &m, size_t i)
{
	return m(i);
}

template<typename M> inline const double &BlockRef(const M &m, size_t i)
{
	return m(i);
}

// get/set matrix
template<typename M> inline double &BlockRef(M &m, size_t i, size_t j)
{
	return m(i, j);
}

template<typename M> inline const double &BlockRef(const M &m, size_t i, size_t j)
{
	return m(i, j);
}

//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables 

	
// MATRICES

// todo: replace add_mult etc. with template expressions
// dest = b*vec
template<typename A, typename B, typename C> inline void AssignMult(A &dest, const B &b, const C &vec);
// dest += b*vec
template<typename A, typename B, typename C> inline void AddMult(A &dest, const B &b, const C &vec);
// dest -= b*vec
template<typename A, typename B, typename C> inline void SubMult(A &dest, const B &b, const C &vec);

// VECTORs



//////////////////////////////////////////////////////
//setSize(t, a, b) for matrices
template<typename T>
inline void SetSize(T &t, size_t a, size_t b);

//setSize(t, a) for vectors
template<typename T>
inline void SetSize(T &t, size_t a);

// getSize
template<typename T>
inline size_t GetSize(T &t);

//getRows
template<typename T>
inline size_t GetRows(const T &t);

// getRows
template<typename T>
inline size_t GetCols(const T &t);


template<typename M>
inline bool GetInverse(typename block_matrix_traits<M>::inverse_type &inv, const M &m);

template<typename M>
inline bool Invert(M &m);


} // namespace ug

#include "double.h"
#include "small_matrix/densevector.h"
#include "small_matrix/densematrix.h"
#include "small_matrix/block_dense.h"
#include "storage/storage.h"
#endif
