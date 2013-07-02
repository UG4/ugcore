/*
 * crouzeix_raviart.h
 *
 * Created on: 18.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE__

#include "../common/lagrange1d.h"
#include "../local_finite_element_provider.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem>
class MiniBubbleLSFS;


/// MiniBubble Set (2D only!)
template <typename TRefElem>
class MiniBubbleLDS : public LocalDoFSet
{
	protected:
	///	dimension of reference element
		static const int refDim = TRefElem::dim;

	public:
	///	constructor
		MiniBubbleLDS()
		{
			// get _the_ reference element
			const TRefElem& rRefElem = Provider<TRefElem>::get();

			if(refDim >= 2)
			{
				// face (or volume???)
				//	set local DoFs (located at vertices+bubble)
				nsh = rRefElem.num(0)+1;

				m_vLocalDoF.resize(nsh);
				for(size_t i = 0; i< nsh-1; ++i)
					m_vLocalDoF[i] = LocalDoF(0, i, 0);

				m_vLocalDoF[nsh-1] = LocalDoF(refDim, nsh-1, 0); // bubble located at element
			}
			else
			{
				// edge or vertex
				nsh = refDim+1;
				m_vLocalDoF.resize(nsh);
				for(size_t i = 0; i< nsh-1; ++i)
					m_vLocalDoF[i] = LocalDoF(0, i, 0);

			}
		}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID type) const
		{
			const int d = ReferenceElementDimension(type);
			if (d==0) return 1;         // vertices
			if (d == refDim)   return 1;    // element
			return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		bool exact_position_available() const {return true;};

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

////////////////////////////////////////////////////////////////////////////////
//	ReferenceEdge
//
//  function space span {1,x}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceEdge>
: public MiniBubbleLDS<ReferenceEdge>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceEdge>, 1>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 1;

	/// Number of shape functions
		static const size_t nsh = 3;

	public:
	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceEdge reference_element_type;

	public:
	///	Constructor
		MiniBubbleLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0;return true;
				case 1:	pos[0] = 1;return true;
				case 2:	pos[0] = 0.5;return true; // bubble
				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceEdge::check_position(x);

			switch(i)
			{
				case 0:	return 1-x[0];
				case 1:	return x[0];
				case 2:	return x[0]*(1-x[0]);
				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceEdge::check_position(x);

			switch(i)
			{
				case 0:	g[0] = -1.0;
				case 1:	g[0] =  1.0;
				case 2:	g[0] =  1.0-2*x[0];
				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

////////////////////////////////////////////////////////////////////////////////
//	ReferenceTriangle
//
//  function space span {1,x,y}+bubble
//
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceTriangle>
: public MiniBubbleLDS<ReferenceTriangle>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceTriangle>, 2>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 2;

	/// Number of shape functions
		static const size_t nsh = 4;
	public:
	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceTriangle reference_element_type;

	public:
	///	Constructor
		MiniBubbleLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.0;
						pos[1] = 0.0; return true;
				case 1:	pos[0] = 0.1;
						pos[1] = 0.0; return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 0.1; return true;
				case 3: pos[0] = 1.0/3.0;
						pos[1] = 1.0/3.0; return true;
				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTriangle::check_position(x);

			switch(i)
			{
				case 0:	return (1.0-x[0]-x[1]);
				case 1:	return x[0];
				case 2:	return x[1];
				case 3: return x[0]*x[1]*(1.0-x[0]-x[1]); // bubble
				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTriangle::check_position(x);

			switch(i)
			{
				case 0:	g[0] = -1.0;
						g[1] = -1.0; return;
				case 1:	g[0] = 1.0;
						g[1] = 0.0; return;
				case 2:	g[0] = 0.0;
						g[1] = 1.0; return;
				case 3: g[0] = x[1]*(1.0-x[1]-2.0*x[0]);
						g[1] = x[0]*(1.0-x[0]-2.0*x[1]); return;
				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


////////////////////////////////////////////////////////////////////////////////
//	ReferenceQuadrilateral
//
//  function space span {1,x,y,x^2-y^2}+bubble
//
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceQuadrilateral>
: public MiniBubbleLDS<ReferenceQuadrilateral>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceQuadrilateral>, 2>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 2;

	/// Number of shape functions
		static const size_t nsh = 5;

	public:
	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceQuadrilateral reference_element_type;

	public:
	///	Constructor
		MiniBubbleLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.0;
						pos[1] = 0.0; return true;
				case 1:	pos[0] = 1.0;
						pos[1] = 0.0; return true;
				case 2:	pos[0] = 1.0;
						pos[1] = 1.0; return true;
				case 3: pos[0] = 0.0;
						pos[1] = 1.0; return true;
				case 4: pos[0] = 0.5;
						pos[1] = 0.5; return true;

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceQuadrilateral::check_position(x);

			switch(i)
			{
				case 0:	return (1.0-x[0])*(1.0-x[1]);
				case 1:	return x[0]*(1.0-x[1]);
				case 2:	return x[0]*x[1];
				case 3:	return x[1]*(1.0-x[0]);
				case 4:	return x[0]*x[1]*(1.0-x[0])*(1.0-x[1]);

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceQuadrilateral::check_position(x);

			switch(i)
			{
				case 0:	g[0] = -1.0*(1.0-x[1]);
						g[1] = -1.0*(1.0-x[0]); return;
				case 1:	g[0] = (1.0-x[1]);
						g[1] = -x[0]; return;
				case 2:	g[0] = x[1];
						g[1] = x[0]; return;
				case 3: g[0] = -x[1];
						g[1] = 1.0-x[0]; return;
				case 4: g[0] = x[1]*(1.0-x[1])*(1.0-2.0*x[0]);
						g[1] = x[0]*(1.0-x[0])*(1.0-2.0*x[1]); return;

				default: UG_THROW("MiniLSFS: shape function "<<i<<
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
class MiniBubbleLSFS<ReferenceTetrahedron>
: public MiniBubbleLDS<ReferenceTetrahedron>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceTetrahedron>, 3>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 3;

	/// Number of shape functions
		static const size_t nsh = 4+1;

	public:
	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceTetrahedron reference_element_type;

	public:
	///	Constructor
		MiniBubbleLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
					case 0:	pos[0] = 0.0; pos[1] = 0.0; pos[2] = 0.0; return true;
					case 1:	pos[0] = 1.0; pos[1] = 0.0; pos[2] = 0.0; return true;
					case 2:	pos[0] = 0.0; pos[1] = 1.0; pos[2] = 0.0; return true;
					case 3: pos[0] = 0.0; pos[1] = 0.0; pos[2] = 1.0; return true;
					case 4: pos[0] = 0.25; pos[1] = 0.25; pos[2] = 0.25; return true;
				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTetrahedron::check_position(x);
			//number prod = x[0]*x[1]*x[2];
			//number lambda4 = (1.0-x[0]-x[1]-x[2]);

			switch(i)
			{
				case 0:	return (1.0-x[0]-x[1]-x[2]);
				case 1:	return x[0];
				case 2:	return x[1];
				case 3:	return x[2];
				case 4: return x[0]*x[1]*x[2]*(1.0-x[0]-x[1]-x[2]); // bubble
				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTetrahedron::check_position(x);
			number prod = x[0]*x[1]*x[2];
			number lambda4 = (1.0-x[0]-x[1]-x[2]);

			switch(i)
			{
				case 0:	g[0] = -1.0;
						g[1] = -1.0;
						g[1] = -1.0; return;
				case 1:	g[0] =  1.0;
						g[1] =  0.0;
						g[2] =  0.0; return;
				case 2:	g[0] =  0.0;
						g[1] =  1.0;
						g[2] =  0.0; return;
				case 3:	g[0] =  0.0;
						g[1] =  0.0;
						g[2] =  1.0; return;
				case 4:	g[0] = -prod + lambda4*x[1]*x[2];
						g[1] = -prod + lambda4*x[0]*x[2];
						g[2] = -prod + lambda4*x[0]*x[2]; return;

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__ */

