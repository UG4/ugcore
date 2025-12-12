/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__

//#include "pcl/pcl.h"
#include "parallel_index_layout.h"
//#include "parallelization_util.h"
#include "parallel_storage_type.h"
#include "algebra_layouts.h"
#include "lib_algebra/common/operations.h"
#include "parallel_vector.h"


namespace ug {

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
		ParallelMatrix(const ParallelMatrix&) = delete;

	public:
	///	own type
		using this_type = ParallelMatrix;

	public:
	///	Default Constructor
		ParallelMatrix()
			: TMatrix(), m_type(PST_UNDEFINED), m_spAlgebraLayouts(new AlgebraLayouts)
		{}

	///	Constructor setting the layouts
		explicit ParallelMatrix(SmartPtr<AlgebraLayouts> layouts)
			: TMatrix(), m_type(PST_UNDEFINED), m_spAlgebraLayouts(layouts)
		{}

		/////////////////////////
		// Storage type handling
		/////////////////////////

	///	returns the algebra layouts
		[[nodiscard]] ConstSmartPtr<AlgebraLayouts> layouts() const {return m_spAlgebraLayouts;}

	///	sets the algebra layouts
		void set_layouts(ConstSmartPtr<AlgebraLayouts> layouts) {m_spAlgebraLayouts = layouts;}

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
		[[nodiscard]] bool has_storage_type(uint type) const
			{return type == PST_UNDEFINED ? m_type == PST_UNDEFINED : (m_type & type) == type;}

	/// returns storage type mask
		[[nodiscard]] uint get_storage_mask() const { return m_type; }
		[[nodiscard]] ParallelStorageType get_storage_type() const {
			return static_cast<ParallelStorageType>(m_type);
		}

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
		this_type &operator = (const this_type &M);

	private:
	/// type of storage  (i.e. consistent, additiv, additiv unique)
		uint m_type;

	/// algebra layouts and communicators
		ConstSmartPtr<AlgebraLayouts> m_spAlgebraLayouts;
};

//	predaclaration.
//	this type may already be declared somewhere else, which shouldn't hurt.
template<typename T>
struct matrix_algebra_type_traits;

template<typename T>
struct matrix_algebra_type_traits<ParallelMatrix<T> >
{
	enum{
		type=MATRIX_USE_GLOBAL_FUNCTIONS
	};
};

} // end namespace ug

#include "parallel_matrix_impl.h"

#endif