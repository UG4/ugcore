/*
 * parallel_vector.h
 *
 *  Created on: 3.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__

#include "pcl/pcl.h"
#include "parallel_index_layout.h"
#include "parallelization_util.h"
#include "parallel_storage_type.h"
#include "algebra_layouts.h"
#include "common/assert.h"

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
template <typename TVector>
class ParallelVector : public TVector
{
	public:
		typedef typename TVector::value_type value_type;
		typedef typename TVector::vector_type vector_type;

	private:
	// 	disallow copy constructor
		//ParallelVector(const ParallelVector&);

	public:
	///	own type
		typedef ParallelVector<TVector> this_type;

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

	/// dotprod (overwrites TVector::dotprod())
	/**
	 * Returns the dot product of the vector. First, the dot prod of each process
	 * is computed using the dotprod() method of the sequential vector. Then,
	 * the results are summarized over all processes.
	 */
		inline number dotprod(const this_type& v);

	/// set all entries to value and the storage type
		bool set(number w, ParallelStorageType type);

	///	sets all entries to a value and the storage type to consistent
		bool set(number w){return set(w, PST_CONSISTENT);}

	/// set all entries to a random number (overwrites TVector::set_random(number from, number to))
		bool set_random(number from, number to, ParallelStorageType type);
		bool set_random(number from, number to){return set_random(from, to, PST_CONSISTENT);}

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
		virtual this_type* virtual_clone() const;

		/// virtual clone using covariant return type excluding values
		virtual this_type* virtual_clone_without_values() const;

	private:
	// 	type of storage  (i.e. consistent, additiv, additiv unique)
	//	holds or-combiation of constants enumerated in ug::ParallelStorageType.
		uint m_type;

	/// algebra layouts and communicators
		ConstSmartPtr<AlgebraLayouts> m_spAlgebraLayouts;
};

} // end namespace ug

#include "parallel_vector_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__ */
