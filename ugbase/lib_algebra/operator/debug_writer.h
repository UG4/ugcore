/*
 * debug_writer.h
 *
 *  Created on: 18.01.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__DEBUG_WRITER__
#define __H__LIB_ALGEBRA__OPERATOR__DEBUG_WRITER__

namespace ug{

/// base class for all debug writer
/**
 * This is the base class for debug output of algebraic vectors and matrices.
 */
template <typename TAlgebra>
class IDebugWriter
{
	public:
	///	type of algebra
		typedef TAlgebra algebra_type;

	///	type of vector
		typedef typename TAlgebra::vector_type vector_type;

	///	type of matrix
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	write vector
		virtual bool write_vector(const vector_type& vec,
		                          const char* name) = 0;

	///	write matrix
		virtual bool write_matrix(const matrix_type& mat,
		                          const char* name) = 0;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__DEBUG_WRITER__ */
