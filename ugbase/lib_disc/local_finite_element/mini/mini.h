/*
 * crouzeix_raviart.h
 *
 * Created on: 18.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE__

#include "../common/lagrange1d.h"
#include "../local_shape_function_set.h"
#include "../local_dof_set.h"
#include "mini_local_dof.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem>
class MiniBubbleLSFS;

////////////////////////////////////////////////////////////////////////////////
//	ReferenceEdge
//
//  function space span {1,x}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceEdge>
: public MiniBubbleLDS<ReferenceEdge>,
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferenceEdge>, 1>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 1;

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
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferenceTriangle>, 2>
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
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferenceQuadrilateral>, 2>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 2;

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
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferenceTetrahedron>, 3>
{
	public:
	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = 3;

	/// Number of shape functions
		static const size_t nsh = 4+1;

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


////////////////////////////////////////////////////////////////////////////////
//	ReferenceHexahedron
//
//  function space span {1,x,y,z,x^2-y^2,y^2-z^2}
//
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceHexahedron>
: public MiniBubbleLDS<ReferenceHexahedron>,
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferenceHexahedron>, 3>
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

				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceHexahedron::check_position(x);

			switch(i)
			{

				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceHexahedron::check_position(x);

			switch(i)
			{

				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
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
class MiniBubbleLSFS<ReferencePrism>
: public MiniBubbleLDS<ReferencePrism>,
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferencePrism>, 3>
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

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePrism::check_position(x);

			switch(i)
			{

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePrism::check_position(x);

			switch(i)
			{

				default: UG_THROW("MiniLSFS: shape function "<<i<<
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
class MiniBubbleLSFS<ReferencePyramid>
: public MiniBubbleLDS<ReferencePyramid>,
  public BaseLocalShapeFunctionSet<MiniBubbleLSFS<ReferencePyramid>, 3>
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
		MiniBubbleLSFS(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline LFEID type() const {return LFEID(LFEID::MINI, dim, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{

				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePyramid::check_position(x);

			switch(i)
			{

				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			// ReferencePyramid::check_position(x);

			switch(i)
			{

				default: UG_THROW("MiniBubbleLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__ */

