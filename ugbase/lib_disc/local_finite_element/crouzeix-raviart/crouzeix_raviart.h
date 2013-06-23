/*
 * crouzeix_raviart.h
 *
 * Created on: 18.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__

#include "../common/lagrange1d.h"
#include "../local_shape_function_set.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem>
class CrouzeixRaviartLSFS;

////////////////////////////////////////////////////////////////////////////////
//	ReferenceEdge
//
//  function space span {1,x}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferenceEdge>
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferenceEdge>, 1>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 1;

	/// Number of shape functions
		static const size_t nsh = 2;

	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceEdge reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0;return true;
				case 1:	pos[0] = 1;return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
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
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
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
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

////////////////////////////////////////////////////////////////////////////////
//	ReferenceTriangle
//
//  function space span {1,x,y}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class CrouzeixRaviartLSFS<ReferenceTriangle>
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferenceTriangle>, 2>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 2;

	/// Number of shape functions
		static const size_t nsh = 3;
	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceTriangle reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
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
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTriangle::check_position(x);

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
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTriangle::check_position(x);

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
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferenceQuadrilateral>, 2>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 2;

	/// Number of shape functions
		static const size_t nsh = 4;

	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceQuadrilateral reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
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
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceQuadrilateral::check_position(x);

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
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceQuadrilateral::check_position(x);

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
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferenceTetrahedron>, 3>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 3;

	/// Number of shape functions
		static const size_t nsh = 4;

	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceTetrahedron reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
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
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTetrahedron::check_position(x);

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
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTetrahedron::check_position(x);

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
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferenceHexahedron>, 3>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 3;

	/// Number of shape functions
		static const size_t nsh = 6;

	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferenceHexahedron reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
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
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceHexahedron::check_position(x);

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
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceHexahedron::check_position(x);

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
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferencePrism>, 3>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 3;

	/// Number of shape functions
		static const size_t nsh = 5;
	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferencePrism reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
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
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePrism::check_position(x);

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
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePrism::check_position(x);

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
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<ReferencePyramid>, 3>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 3;

	/// Number of shape functions
		static const size_t nsh = 5;

	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Reference Element type
		typedef ReferencePyramid reference_element_type;

	public:
	///	Constructor
		CrouzeixRaviartLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.0;
						pos[1] = 0.0;
						pos[2] = 0.0;return true;
				case 1:	pos[0] = 1.0;
						pos[1] = 0.0;
						pos[2] = 0.0; return true;
				case 2:	pos[0] = 1.0;
						pos[1] = 1.0;
						pos[2] = 0.0; return true;
				case 3: pos[0] = 0.0;
						pos[1] = 1.0;
						pos[2] = 0.0; return true;
				case 4: pos[0] = 0.0;
						pos[1] = 0.0;
						pos[2] = 1.0; return true;
				default: UG_THROW("CrouzeixRaviartLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePyramid::check_position(x);

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
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePyramid::check_position(x);

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

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__ */

