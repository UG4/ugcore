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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__



#include "pcl/pcl.h"
#include "parallel_index_layout.h"
#include "parallelization_util.h"
#include "parallel_storage_type.h"
#include "algebra_layouts.h"
#include "common/assert.h"
#include "lib_algebra/cpu_algebra/vector.h"

namespace ug
{

///\ingroup lib_algebra_parallelization

/// Parallelization Wrapper to make a Vector usable in parallel
/**
 * A ParallelVector is a wrapper around a sequential vector to make it usable
 * in parallel. It has all the function a sequential vector supports, since it
 * is derived from it. Furthermore the ParallelStorageType is remembered and
 * can be switched. In addition some functions of the sequential vector are
 * overwritten to adapted the functionality to parallel (e.g. norm, set).
 *
 * \tparam 	TVector		Sequential Vector type
 */
template <typename TVector = Vector<double>>
class ParallelVector : public TVector
{
	public:
		using value_type = typename TVector::value_type;
		using size_type = size_t;
		using vector_type = typename TVector::vector_type;
		///	own type
		using this_type = ParallelVector;

	private:
	// 	disallow copy constructor
		//ParallelVector(const ParallelVector&);

	/// catch all other assignments so we don't use copy constructor here
		template<typename T> this_type &operator =(T t);


	public:
	///	Default constructor
		ParallelVector()
			: TVector(), m_type(PST_UNDEFINED), m_spAlgebraLayouts(new AlgebraLayouts)
		{}

	/// Resizing constructor
		ParallelVector(size_t length)
			: TVector(length), m_type(PST_UNDEFINED), m_spAlgebraLayouts(new AlgebraLayouts)
		{}

	///	Constructor setting the Layouts
		ParallelVector(ConstSmartPtr<AlgebraLayouts> layouts)
			: TVector(), m_type(PST_UNDEFINED), m_spAlgebraLayouts(layouts)
		{}

		~ParallelVector() override = default;

		/////////////////////////
		// Storage type handling
		/////////////////////////

	///	returns the algebra layouts
		ConstSmartPtr<AlgebraLayouts> layouts() const {return m_spAlgebraLayouts;}

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
		bool has_storage_type(uint type) const
			{return type == PST_UNDEFINED ? m_type == PST_UNDEFINED : (m_type & type) == type;}

	/// returns storage type mask
		uint get_storage_mask() const { return m_type; }
		ParallelStorageType get_storage_type() const { return (ParallelStorageType) m_type; }

	///	checks correctness of storage type
		void check_storage_type() const;

	///	sets storage type to consistent
		void enforce_consistent_type();

		//////////////////////////////////////////////////
		// overwritten functions of sequential vector
		//////////////////////////////////////////////////

	/// two norm (overwrites TVector::norm())
	/**
	 * Returns the two norm of the Vector. First, the two norm of each process
	 * is computed using the norm() method of the sequential vector. Then,
	 * the norms are summarized over all processes.
	 */
		number norm() const;

	/// max norm (overwrites TVector::maxnorm())
		number maxnorm() const;

	/// dotprod (overwrites TVector::dotprod())
	/**
	 * Returns the dot product of the vector. First, the dot prod of each process
	 * is computed using the dotprod() method of the sequential vector. Then,
	 * the results are summarized over all processes.
	 */
		inline number dotprod(const this_type& v);

	/// assign number to whole Vector
		number operator = (number d);

	/// set all entries to value and the storage type
		void set(number w, ParallelStorageType type);

	///	sets all entries to a value and the storage type to consistent
		void set(number w){return set(w, PST_CONSISTENT);}

	/// set all entries to a random number (overwrites TVector::set_random(number from, number to))
		void set_random(number from, number to, ParallelStorageType type);
		void set_random(number from, number to){return set_random(from, to, PST_CONSISTENT);}

	///	assignment
		this_type &operator =(const this_type &v);


	///	subtract a vector
		this_type &operator -=(const this_type &v);

	///	add a vector
		this_type &operator +=(const this_type &v);

	/// clones the vector (deep-copy) including values
		SmartPtr<this_type> clone() const;

	/// clones the vector (deep-copy) excluding values
		SmartPtr<this_type> clone_without_values() const;

	protected:
		/// virtual clone using covariant return type
		this_type* virtual_clone() const override;

		/// virtual clone using covariant return type excluding values
		this_type* virtual_clone_without_values() const override;

	private:
	// 	type of storage  (i.e. consistent, additive, additive unique)
	//	holds or-combiation of constants enumerated in ug::ParallelStorageType.
		uint m_type;

	/// algebra layouts and communicators
		ConstSmartPtr<AlgebraLayouts> m_spAlgebraLayouts;
};

} // end namespace ug

#include "parallel_vector_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__ */
