/*
 * parallel_matrix.h
 *
 *  Created on: 19.10.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__

#include "pcl/pcl.h"
#include "parallel_index_layout.h"
#include "parallelization_util.h"
#include "parallel_storage_type.h"
#include "algebra_layouts.h"
#include "lib_algebra/common/operations.h"

namespace ug
{

///\ingroup lib_algebra_parallelization

///\brief Wrapper for sequential matrices to handle them in parallel
/**
 * A ParallelMatrix is a wrapper around a sequential vector to make it usable
 * in parallel. It has all the function a sequential matrix supports, since it
 * is derived from it. Furthermore the ParallelStorageType is remembered.
 * Currently only additive storage type is implemented. In addition some
 * functions of the sequential matrix are overwritten to adapted the functionality
 * to parallel (e.g. set)
 *
 * Please Note:
 * The Implementation is not yet finished.
 *
 *\tparam		TMatrix		Sequential Matrix Type
 */
template <typename TMatrix>
class ParallelMatrix : public TMatrix
{
	public:
		enum {rows_sorted=TMatrix::rows_sorted};


	private:
	// 	disallow copy constructor
		ParallelMatrix(const ParallelMatrix&);

	public:
	///	own type
		typedef ParallelMatrix<TMatrix> this_type;

	public:
	///	Default Constructor
		ParallelMatrix()
			: TMatrix(), m_type(PST_UNDEFINED), m_spAlgebraLayouts(new HorizontalAlgebraLayouts)
		{}

	///	Constructor setting the layouts
		ParallelMatrix(SmartPtr<HorizontalAlgebraLayouts> layouts)
			: TMatrix(), m_type(PST_UNDEFINED), m_spAlgebraLayouts(layouts)
		{}

		/////////////////////////
		// Storage type handling
		/////////////////////////

	///	returns the algebra layouts
		ConstSmartPtr<HorizontalAlgebraLayouts> layouts() const {return m_spAlgebraLayouts;}

	///	sets the algebra layouts
		void set_layouts(ConstSmartPtr<HorizontalAlgebraLayouts> layouts) {m_spAlgebraLayouts = layouts;}

	/// sets the storage type
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		void set_storage_type(uint type) {m_type = type;}

	/// adds a storage type
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		void add_storage_type(uint type) {m_type |= type;}

	/// removes a storage type
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		void remove_storage_type(uint type) {m_type &= ~type;}

	/// changes to the requested storage type if possible
		bool change_storage_type(ParallelStorageType type);

	/// returns if the current storage type has a given representation
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		bool has_storage_type(uint type) const
			{return type == PST_UNDEFINED ? m_type == PST_UNDEFINED : (m_type & type) == type;}

	/// returns storage type mask
		uint get_storage_mask() const { return m_type; }

		/////////////////////////
		// OverWritten functions
		/////////////////////////

	/// calculate res = A x
		template<typename TPVector>
		bool apply(TPVector &res, const TPVector &x) const;

	/// calculate res = A.T x
		template<typename TPVector>
		bool apply_transposed(TPVector &res, const TPVector &x) const;

	/// calculate res -= A x
		template<typename TPVector>
		bool matmul_minus(TPVector &res, const TPVector &x) const;

	///	assignment
		this_type &operator =(const this_type &M);

	private:
	/// type of storage  (i.e. consistent, additiv, additiv unique)
		uint m_type;

	/// algebra layouts and communicators
		ConstSmartPtr<HorizontalAlgebraLayouts> m_spAlgebraLayouts;
};

//	predaclaration.
//	this type may already be declared somewhere else, which shouldn't hurt.
template<typename T>
struct matrix_algebra_type_traits;

template<typename T>
struct matrix_algebra_type_traits<ParallelMatrix<T> >
{
	static const int type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

} // end namespace ug

#include "parallel_matrix_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__ */
