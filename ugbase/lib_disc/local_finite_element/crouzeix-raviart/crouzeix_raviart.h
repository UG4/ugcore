/*
 * crouzeix_raviart.h
 *
 * Created on: 18.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__

#include "lib_grid/grid/grid_base_objects.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"

namespace ug{


/// Crouzeix - Raviart Set
template <typename TRefElem>
class CrouzeixRaviartBase
{
	public:
	///	dimension of reference element
		static const int dim = TRefElem::dim;

	/// Number of shape functions
		static const size_t nsh = TRefElem::numSides;

	public:
	///	constructor
		CrouzeixRaviartBase()
		{
			for(size_t i = 0; i < nsh; ++i)
				m_vLocalDoF[i] = LocalDoF(dim-1, i, 0);
		}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};
		size_t num_sh() const {return nsh;}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID type) const
		{
			if(ReferenceElementDimension(type) == dim-1) return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		bool exact_position_available() const {return true;};

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		bool continuous() const {return false;}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];
};

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem>
class CrouzeixRaviartLSFS;

////////////////////////////////////////////////////////////////////////////////
//	ReferenceTriangle
//
//  function space span {1,x,y}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferenceTriangle>
: public CrouzeixRaviartBase<ReferenceTriangle>,
  public BaseLSFS<CrouzeixRaviartLSFS<ReferenceTriangle>, 2>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = 2;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.0; return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.5; return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 0.5; return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	return 1-2*x[1];
				case 1:	return -1+2*x[0]+2*x[1];
				case 2:	return 1-2*x[0];
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	g[0] = 0.0;
						g[1] = -2.0; return;
				case 1:	g[0] = 2.0;
						g[1] = 2.0; return;
				case 2:	g[0] = -2.0;
						g[1] = 0.0; return;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


////////////////////////////////////////////////////////////////////////////////
//	ReferenceQuadrilateral
//
//  function space span {1,x,y,x^2-y^2}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferenceQuadrilateral>
: public CrouzeixRaviartBase<ReferenceQuadrilateral>,
  public BaseLSFS<CrouzeixRaviartLSFS<ReferenceQuadrilateral>, 2>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = 2;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.0; return true;
				case 1:	pos[0] = 1.0;
						pos[1] = 0.5; return true;
				case 2:	pos[0] = 0.5;
						pos[1] = 1.0; return true;
				case 3: pos[0] = 0.0;
						pos[1] = 0.5; return true;
						return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	return 0.75+x[0]-2*x[1]-x[0]*x[0]+x[1]*x[1];
				case 1:	return -0.25+x[1]+x[0]*x[0]-x[1]*x[1];
				case 2:	return -0.25+x[0]-x[0]*x[0]+x[1]*x[1];
				case 3: return 0.75-2*x[0]+x[1]+x[0]*x[0]-x[1]*x[1];
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	g[0] = 1-2*x[0];
						g[1] = -2.0+2*x[1]; return;
				case 1:	g[0] = 2.0*x[0];
						g[1] = 1-2.0*x[1]; return;
				case 2:	g[0] = 1-2.0*x[0];
						g[1] = 2*x[1]; return;
				case 3: g[0] = -2+2*x[0];
				 	 	g[1] = 1-2*x[1]; return;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

////////////////////////////////////////////////////////////////////////////////
//	ReferenceTetrahedron
//
//  function space span {1,x,y,z}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferenceTetrahedron>
: public CrouzeixRaviartBase<ReferenceTetrahedron>,
  public BaseLSFS<CrouzeixRaviartLSFS<ReferenceTetrahedron>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceTetrahedron::dim;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 1.0/3.0;
						pos[1] = 1.0/3.0;
						pos[2] = 0.0;		return true;
				case 1:	pos[0] = 1.0/3.0;
						pos[1] = 1.0/3.0;
						pos[2] = 1.0/3.0; return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 1.0/3.0;
						pos[2] = 1.0/3.0; return true;
				case 3: pos[0] = 1.0/3.0;
						pos[1] = 0.0;
						pos[2] = 1.0/3.0; return true;
						return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	return 1 - 3*x[2];;
				case 1:	return -2 + 3*x[0]+3*x[1]+3*x[2];
				case 2:	return 1 - 3*x[0];
				case 3: return 1 - 3*x[1];
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	g[0] = 0.0;
						g[1] = 0.0;
						g[2] = -3.0;return;
				case 1:	g[0] = 3.0;
						g[1] = 3.0;
						g[2] = 3.0;return;
				case 2:	g[0] = -3.0;
						g[1] = 0.0;
						g[2] = 0.0;return;
				case 3: g[0] = 0.0;
						g[1] = -3;
						g[2] = 0.0;return;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


////////////////////////////////////////////////////////////////////////////////
//	ReferenceHexahedron
//
//  function space span {1,x,y,z,x^2-y^2,y^2-z^2}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferenceHexahedron>
: public CrouzeixRaviartBase<ReferenceHexahedron>,
  public BaseLSFS<CrouzeixRaviartLSFS<ReferenceHexahedron>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceHexahedron::dim;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.5;
						pos[2] = 0.0;return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.0;
						pos[2] = 0.5; return true;
				case 2:	pos[0] = 1.0;
						pos[1] = 0.5;
						pos[2] = 0.5; return true;
				case 3: pos[0] = 0.5;
						pos[1] = 1.0;
						pos[2] = 0.5; return true;
				case 4: pos[0] = 0.0;
						pos[1] = 0.5;
						pos[2] = 0.5; return true;
				case 5: pos[0] = 0.5;
						pos[1] = 0.5;
						pos[2] = 1.0; return true;
						return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	return 1.0/3.0*(2 + 2*x[0] + 2*x[1] - 7*x[2]-2*x[0]*x[0]-2*x[1]*x[1]+4*x[2]*x[2]);
				case 1:	return 1.0/3.0*(2+2*x[0]-7*x[1]+2*x[2]-2*x[0]*x[0]+4*x[1]*x[1]-2*x[2]*x[2]);
				case 2:	return 1.0/3.0*(-1-x[0]+2*x[1]+2*x[2]+4*x[0]*x[0]-2*x[1]*x[1]-2*x[2]*x[2]);
				case 3: return 1.0/3.0*(-1+2*x[0]-x[1]+2*x[2]-2*x[0]*x[0]+4*x[1]*x[1]-2*x[2]*x[2]);
				case 4: return 1.0/3.0*(2-7*x[0]+2*x[1]+2*x[2]+4*x[0]*x[0]-2*x[1]*x[1]-2*x[2]*x[2]);
				case 5: return 1.0/3.0*(-1+2*x[0]+2*x[1]-x[2]-2*x[0]*x[0]-2*x[1]*x[1]+4*x[2]*x[2]);
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	g[0] = 1.0/3.0*(2-4*x[0]);
						g[1] = 1.0/3.0*(2-4*x[1]);
						g[2] = 1.0/3.0*(-7+8*x[2]);return;
				case 1:	g[0] = 1.0/3.0*(2-4*x[0]);
						g[1] = 1.0/3.0*(-7+8*x[1]);
						g[2] = 1.0/3.0*(2-4*x[2]);return;
				case 2:	g[0] = 1.0/3.0*(-1+8*x[0]);
						g[1] = 1.0/3.0*(2-4*x[1]);
						g[2] = 1.0/3.0*(2-4*x[2]);return;
				case 3: g[0] = 1.0/3.0*(2-4*x[0]);
						g[1] = 1.0/3.0*(-1+8*x[1]);
						g[2] = 1.0/3.0*(2-4*x[2]);return;
				case 4: g[0] = 1.0/3.0*(-7+8*x[0]);
						g[1] = 1.0/3.0*(2-4*x[1]);
						g[2] = 1.0/3.0*(2-4*x[2]);return;
				case 5: g[0] = 1.0/3.0*(2-4*x[0]);
						g[1] = 1.0/3.0*(2-4*x[1]);
						g[2] = 1.0/3.0*(-1+8*x[2]);return;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

////////////////////////////////////////////////////////////////////////////////
//	ReferencePrism
//
//  function space span {1,x,y,z,x^2-y^2-z^2}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferencePrism>
: public CrouzeixRaviartBase<ReferencePrism>,
  public BaseLSFS<CrouzeixRaviartLSFS<ReferencePrism>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = 3;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 1.0/3.0;
						pos[1] = 1.0/3.0;
						pos[2] = 0.0;return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.0;
						pos[2] = 0.5; return true;
				case 2:	pos[0] = 0.5;
						pos[1] = 0.5;
						pos[2] = 0.5; return true;
				case 3: pos[0] = 0.0;
						pos[1] = 0.5;
						pos[2] = 0.5; return true;
				case 4: pos[0] = 1.0/3.0;
						pos[1] = 1.0/3.0;
						pos[2] = 1.0; return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	return 1+x[0]-x[1]-3*x[2]-2*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
				case 1:	return 1.0/3.0*(2-2*x[0]-4*x[1]+4*x[2]+4*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]));
				case 2:	return 1.0/3.0*(-4+4*x[0]+8*x[1]+4*x[2]+4*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]));
				case 3:	return 1.0/3.0*(2-8*x[0]+2*x[1]+4*x[2]+4*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]));
				case 4:	return x[0]-x[1]-x[2]-2*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	g[0] = 1-4*x[0];
						g[1] = -1+4*x[1];
						g[2] = -3+4*x[2];return;
				case 1:	g[0] = 1.0/3.0*(-2+8*x[0]);
						g[1] = 1.0/3.0*(-4-8*x[1]);
						g[2] = 1.0/3.0*(4-8*x[2]);return;
				case 2:	g[0] = 1.0/3.0*(4+8*x[0]);
						g[1] = 1.0/3.0*(8-8*x[1]);
						g[2] = 1.0/3.0*(4-8*x[2]);return;
				case 3: g[0] = 1.0/3.0*(-8+8*x[0]);
						g[1] = 1.0/3.0*(2-8*x[1]);
						g[2] = 1.0/3.0*(4-8*x[2]);return;
				case 4: g[0] = 1-4*x[0];
						g[1] = -1+4*x[1];
						g[2] = -1+4*x[2];return;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


////////////////////////////////////////////////////////////////////////////////
//	ReferencePyramid
//
//  function space span {1,x,y,z,x^2-y^2-z^2}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferencePyramid>
: public CrouzeixRaviartBase<ReferencePyramid>,
  public BaseLSFS<CrouzeixRaviartLSFS<ReferencePyramid>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = 3;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.5;
						pos[2] = 0.0;return true;
				case 1:	pos[0] = 1.0/3.0;
						pos[1] = 0.0;
						pos[2] = 1.0/3.0; return true;
				case 2:	pos[0] = 2.0/3.0;
						pos[1] = 1.0/3.0;
						pos[2] = 1.0/3.0; return true;
				case 3: pos[0] = 1.0/3.0;
						pos[1] = 2.0/3.0;
						pos[2] = 1.0/3.0; return true;
				case 4: pos[0] = 0.0;
						pos[1] = 1.0/3.0;
						pos[2] = 1.0/3.0; return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	return 1-3*x[2];
				case 1:	return 0.75+1.5*x[0]-3*x[1]-0.75*x[2]-2.25*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
				case 2:	return -0.75+1.5*x[1]+2.25*x[2]+2.25*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
				case 3:	return -0.75+1.5*x[0]+0.75*x[2]-2.25*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
				case 4:	return 0.75-3*x[0]+1.5*x[1]+0.75*x[2]+2.25*(x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0:	g[0] = 0.0;
						g[1] = 0.0;
						g[2] = -3.0;return;
				case 1:	g[0] = 1.5-4.5*x[0];
						g[1] = -3+4.5*x[1];
						g[2] = -0.75+4.5*x[2];return;
				case 2:	g[0] = 4.5*x[0];
						g[1] = 1.5-4.5*x[1];
						g[2] = 2.25-4.5*x[2];return;
				case 3: g[0] = 1.5-4.5*x[0];
						g[1] = 4.5*x[1];
						g[2] = 0.75+4.5*x[2];return;
				case 4: g[0] = -3+4.5*x[0];
						g[1] = 1.5-4.5*x[1];
						g[2] = 0.75-4.5*x[2];return;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__ */

