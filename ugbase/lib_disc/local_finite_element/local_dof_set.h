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

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_DOF_SET__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_DOF_SET__

// #include <vector>
// #include <map>

// #include "local_finite_element_id.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_grid/grid/grid_base_objects.h"

namespace ug {

/// \ingroup lib_disc_local_finite_elements
/// @{

/**
 * This class is used to store for a single degree of freedom (DoF) the location
 * within an element. For continuous finite elements the DoFs are usually
 * associated with a sub-geometric object of the element itself (e.g. a vertex).
 * This can be requested from this class, which stores the dimension of the
 * sub-element the DoF is located on, the id of the sub-element (w.r.t. to the
 * numbering in the reference elements) and an offset > 0 if there are more than
 * one DoFs associated with the same sub-element.
 */
class LocalDoF
{
	public:
	///	default constructor
		LocalDoF() : m_dim(-1), m_id(0), m_offset(0) {}

	///	constructor
	/**
	 * Create a pair describing the position of a DoF within the reference element.
	 *
	 * \param[in]	dim		dimension of sub-geometric object
	 * \param[in]	id		number of sub-geometric object (in the numbering
	 * 						used by the reference element)
	 * \param[in]	offset	if several DoFs are associated with the same
	 * 						sub-geometric object the offset specifies the number
	 * 						within all DoFs on that geometric object
	 */
		LocalDoF(int dim, size_t id, size_t offset)
			: m_dim(dim), m_id(id), m_offset(offset)
		{}

	///	sets the values
		void set(int dim, size_t id, size_t offset)
		{
			m_dim = dim; m_id = id; m_offset = offset;
		}

	///	returns the dimension of associated geometric object
		inline int dim() const {return m_dim;}

	///	returns the index for the geometric object (w.r.t reference element numbering)
		inline size_t id() const {return m_id;}

	///	returns the offset for the geometric object
		inline size_t offset() const {return m_offset;}

	///	equality check
		bool operator == (const LocalDoF& v) const{
			return dim() == v.dim() && id() == v.id() && offset() == v.offset();
		}

	///	inequality check
		bool operator != (const LocalDoF& v) const {return !((*this)==v);}

	protected:
	///	dimension of sub-geometric object
		int m_dim;

	///	id of sub-geometric object in counting of reference element
		size_t m_id;

	///	offset if several DoFs associated to the same geometric object
		size_t m_offset;
};

/// writes to the output stream
std::ostream& operator << (std::ostream& out,	const LocalDoF& v);

/**
 * This class provides the interface for the storage of degrees of freedom
 * on a finite element.
 */
class LocalDoFSet
{
	public:
	///	returns the reference dimension
		int dim() const;

	///	returns the Reference object id of the corresponding grid object
		virtual ReferenceObjectID_t roid() const = 0;

	///	returns the total number of dofs on the finite element
	/// \{
				size_t num_dof() const {return num_sh();}
		virtual size_t num_sh() const;
	/// \}

	///	returns the number of DoFs on a sub-geometric object type
		virtual size_t num_dof(ReferenceObjectID_t roid) const = 0;

	///	returns the DoFs storage
		virtual const LocalDoF& local_dof(size_t dof) const = 0;

	///	returns the number of DoFs on a sub-geometric object of dim and id
		size_t num_dof(int d, size_t id) const;

	///	equality check
		bool operator == (const LocalDoFSet& v) const;

	///	inequality check
		bool operator != (const LocalDoFSet& v) const {return !((*this)==v);}

	///	virtual destructor
		virtual ~LocalDoFSet() = default;
};

/**
 * Local DoF Set also providing to local position of the dofs (iff available)
 */
template <int TDim>
class DimLocalDoFSet : public LocalDoFSet
{
	public:
	/// constructor
		DimLocalDoFSet() = default;

	///	returns if the local dof position are exact
		[[nodiscard]] virtual bool exact_position_available() const = 0;

	///	local position of DoF i
	/**
	 * This function returns the local position of a DoF if possible.
	 * \param[in] 	i		number of DoF
	 * \param[out]	pos		Position of DoF
	 * \retval		true 	if position exists
	 * \retval		false 	if no meaningful position available
	 */
		virtual bool position(size_t i, MathVector<TDim>& pos) const = 0;

	///	equality check
		bool operator == (const DimLocalDoFSet<TDim>& v) const;

	///	inequality check
		bool operator != (const DimLocalDoFSet<TDim>& v) const {return !((*this)==v);}
};

/// @}

/// writes to the output stream
std::ostream& operator << (std::ostream& out, const LocalDoFSet& v);
/// writes to the output stream
template <int dim>
std::ostream& operator << (std::ostream& out, const DimLocalDoFSet<dim>& v);

/**
 * Intersection of local dof sets
 */
class CommonLocalDoFSet
{
	public:
	///	indicate not set value
		enum{NOT_SPECIFIED = -1};

	///	constructor
		CommonLocalDoFSet() {clear();}

	///	reset all numbers of dofs to not set
		void clear();

	///	add a local dof set to the intersection
		void add(const LocalDoFSet& set);

	///	number of dofs on a reference element type
		[[nodiscard]] int num_dof(ReferenceObjectID_t roid) const {return m_vNumDoF[roid];}

	protected:
		int m_vNumDoF[NUM_REFERENCE_OBJECTS];
};

/// writes to the output stream
std::ostream& operator << (std::ostream& out, const CommonLocalDoFSet& v);


} // end namespace ug

#endif