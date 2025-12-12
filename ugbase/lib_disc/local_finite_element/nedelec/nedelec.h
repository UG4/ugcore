/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/*
 * This file contains implementations of the local shape function set for
 * the so-called Nedelec (or Whitney-1) elements.
 */
#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__

#include "common/util/provider.h"
#include "../local_dof_set.h"
// #include "lib_disc/common/multi_index.h"

namespace ug {

///////////////////////////////////////////////////////////////////////////////
// Nedelec Set
///////////////////////////////////////////////////////////////////////////////

/// Nedelec, i.e. the edge local dof set
template <typename TRefElem>
class NedelecLDS
{
	protected:
	///	dimension of reference element
		static constexpr int refDim = TRefElem::dim;

	public:
	///	constructor
		NedelecLDS()
		{
			if(refDim < 2)
			{
			//	No dofs if the dimension is less than 2:
				nsh = 0;
				m_vLocalDoF.clear();
				return;
			}

			const TRefElem& rRefElem = Provider<TRefElem>::get();

			nsh = rRefElem.num(1); // number of the edges
		//	set local DoFs (all located at the edges)
			m_vLocalDoF.resize(nsh);
			for(size_t i = 0; i < nsh; ++i)
				m_vLocalDoF[i] = LocalDoF(1, i, 0);
		}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID type) const
		{
			if(ReferenceElementDimension(type) == 1) return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		bool exact_position_available() const {return true;};

	protected:
	///	number of shapes (== number of edges)
		size_t nsh;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

/**
 * Nedelec (or Whitney-1) base function set for a general element:
 * Not implemented, so this class implements error messages only.
 * For the Nedelec base functions for triangles and tetrahedra cf. the
 * specializations below.
 */
template <typename TRefElement>
class NedelecLSFS;

/// Nedelec (or Whitney-1) base function set for triangles
template <>
class NedelecLSFS<ReferenceTriangle>
: public NedelecLDS<ReferenceTriangle>,
  public
	BaseLSFS
		<
			NedelecLSFS<ReferenceTriangle>,
			ReferenceTriangle::dim, ///< dimensionality of the element
			MathVector<ReferenceTriangle::dim>, ///< return type of the shape functions
			MathMatrix<ReferenceTriangle::dim, ReferenceTriangle::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
		using reference_element_type = ReferenceTriangle;

	///	Order of Shape functions
		static constexpr size_t order = 1;

	///	Dimension, where shape functions are defined
		static constexpr int dim = reference_element_type::dim;

	private:
	///	Base class
	using base_type = BaseLSFS;

	public:
	///	Shape type
	using shape_type = shape_type;

	///	Gradient type
	using grad_type = grad_type;

	protected:
	///	number of shapes
		static constexpr size_t nsh = reference_element_type::numEdges;

	public:
	///	Constructor
		NedelecLSFS() {};

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.0;
						return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.5;
						return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 0.5;
						return true;
				default: UG_THROW("NedelecLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline MathVector<dim> shape(const size_t i, const MathVector<dim>& x) const
		{
			UG_THROW ("NedelecLSFS: Nedelec shapes cannot be computed in the reference space.");
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathMatrix<dim,dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			UG_THROW ("NedelecLSFS: Gradients of the Nedelec shapes cannot be computed in the reference space.");
		}
};

/// Nedelec (or Whitney-1) base function set for tetrahedra
template <>
class NedelecLSFS<ReferenceTetrahedron>
: public NedelecLDS<ReferenceTetrahedron>,
  public
	BaseLSFS
		<
			NedelecLSFS<ReferenceTetrahedron>,
			ReferenceTetrahedron::dim, ///< dimensionality of the element
			MathVector<ReferenceTetrahedron::dim>, ///< return type of the shape functions
			MathMatrix<ReferenceTetrahedron::dim, ReferenceTetrahedron::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
	using reference_element_type = ReferenceTetrahedron;

	///	Order of Shape functions
		static constexpr size_t order = 1;

	///	Dimension, where shape functions are defined
		static constexpr int dim = reference_element_type::dim;

	private:
	///	Base class
	using base_type = BaseLSFS;

	public:
	///	Shape type
	using shape_type = shape_type;

	///	Gradient type
	using grad_type = grad_type;

	protected:
	///	number of shapes
		static constexpr size_t nsh = reference_element_type::numEdges;

	public:
	///	Constructor
		NedelecLSFS() = default;

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.0;
						pos[2] = 0.0;
						return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.5;
						pos[2] = 0.0;
						return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 0.5;
						pos[2] = 0.0;
						return true;
				case 3:	pos[0] = 0.0;
						pos[1] = 0.0;
						pos[2] = 0.5;
						return true;
				case 4:	pos[0] = 0.5;
						pos[1] = 0.0;
						pos[2] = 0.5;
						return true;
				case 5:	pos[0] = 0.0;
						pos[1] = 0.5;
						pos[2] = 0.5;
						return true;
				default: UG_THROW("NedelecLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline MathVector<dim> shape(const size_t i, const MathVector<dim>& x) const
		{
			UG_THROW ("NedelecLSFS: Nedelec shapes cannot be computed in the reference space.");
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathMatrix<dim,dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			UG_THROW ("NedelecLSFS: Gradients of the Nedelec shapes cannot be computed in the reference space.");
		}
};

} // namespace ug

#endif