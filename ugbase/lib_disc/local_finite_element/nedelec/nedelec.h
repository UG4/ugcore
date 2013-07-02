/*
 * nedelec.h
 *
 * This file contains implementations of the local shape function set for
 * the so-called Nedelec (or Whitney-1) elements.
 *
 * Created on: 08.02.2013
 * Author: Dmitry Logashenko
 */
#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__

#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Nedelec Set
///////////////////////////////////////////////////////////////////////////////

/// Nedelec, i.e. the edge local dof set
template <typename TRefElem>
class NedelecLDS
{
	protected:
	///	dimension of reference element
		static const int refDim = TRefElem::dim;

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
	BaseLocalShapeFunctionSet
		<
			NedelecLSFS<ReferenceTriangle>,
			ReferenceTriangle::dim, ///< dimensionality of the element
			MathVector<ReferenceTriangle::dim>, ///< return type of the shape functions
			MathMatrix<ReferenceTriangle::dim, ReferenceTriangle::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
		typedef ReferenceTriangle reference_element_type;

	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	private:
	///	Base class
		typedef BaseLocalShapeFunctionSet<NedelecLSFS<reference_element_type>, dim, MathVector<dim>,  MathMatrix<dim,dim> > base_type;

	public:
	///	Shape type
		typedef base_type::shape_type shape_type;

	///	Gradient type
		typedef base_type::grad_type grad_type;

	protected:
	///	number of shapes
		static const size_t nsh = reference_element_type::numEdges;

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
	BaseLocalShapeFunctionSet
		<
			NedelecLSFS<ReferenceTetrahedron>,
			ReferenceTetrahedron::dim, ///< dimensionality of the element
			MathVector<ReferenceTetrahedron::dim>, ///< return type of the shape functions
			MathMatrix<ReferenceTetrahedron::dim, ReferenceTetrahedron::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
		typedef ReferenceTetrahedron reference_element_type;

	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	private:
	///	Base class
		typedef BaseLocalShapeFunctionSet<NedelecLSFS<reference_element_type>, dim, MathVector<dim>,  MathMatrix<dim,dim> > base_type;

	public:
	///	Shape type
		typedef base_type::shape_type shape_type;

	///	Gradient type
		typedef base_type::grad_type grad_type;

	protected:
	///	number of shapes
		static const size_t nsh = reference_element_type::numEdges;

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

#endif // __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__

/* End of File */
