/*
 * Copyright (c) 2012-2016:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE__

#include "../common/lagrange1d.h"
#include "../local_finite_element_provider.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
// #include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/grid_base_objects.h"
// #include "common/util/provider.h"
#include "lib_disc/reference_element/reference_element_util.h"
// #include "lib_disc/common/multi_index.h"

namespace ug {

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem>
class MiniBubbleLSFS;


/// MiniBubble Set (2D only!)
template <typename TRefElem>
class MiniBubbleLDS : public LocalDoFSet
{
	protected:
	///	dimension of reference element
		static constexpr int refDim = TRefElem::dim;

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
		[[nodiscard]] ReferenceObjectID_t roid() const override {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		[[nodiscard]] size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		[[nodiscard]] size_t num_dof(ReferenceObjectID_t type) const override {
			const int d = ReferenceElementDimension(type);
			if (d==0) return 1;         // vertices
			if (d == refDim)   return 1;    // element
			return 0;
		}

	///	returns the dof storage
		[[nodiscard]] const LocalDoF& local_dof(size_t dof) const override {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		[[nodiscard]] bool exact_position_available() const {return true;};

	protected:
	///	number of shapes
		size_t nsh;

	///	order
	///	size_t p;

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
		static constexpr size_t order = 1;

	///	Dimension, where shape functions are defined
		static constexpr int dim = 1;

	/// Number of shape functions
		static constexpr size_t nsh = 3;

	public:
	///	Shape type
		using shape_type = number;

	///	Gradient type
		using grad_type = MathVector<dim>;

	///	Reference Element type
		using reference_element_type = ReferenceEdge;

	public:
	///	Constructor
		MiniBubbleLSFS()= default;

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		[[nodiscard]] inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		[[nodiscard]] inline size_t num_sh() const override {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		[[nodiscard]] inline bool position(size_t i, MathVector<dim>& pos) const
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
		[[nodiscard]] inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceEdge::check_position(x);

			switch(i)
			{
				case 0:	return 1-x[0];
				case 1:	return x[0];
				case 2:	return 2.0*x[0]*(1-x[0]);
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
				case 2:	g[0] =  2.0-4.0*x[0];
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
		static constexpr size_t order = 2;

	///	Dimension, where shape functions are defined
		static constexpr int dim = 2;

	/// Number of shape functions
		static constexpr size_t nsh = 4;
	public:
	///	Shape type
	using shape_type = number;

	///	Gradient type
	using grad_type = MathVector<dim>;

	///	Reference Element type
	using reference_element_type = ReferenceTriangle;

	public:
	///	Constructor
		MiniBubbleLSFS()= default;

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		[[nodiscard]] inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		[[nodiscard]] inline size_t num_sh() const override {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		[[nodiscard]] inline bool position(size_t i, MathVector<dim>& pos) const
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
		[[nodiscard]] inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceTriangle::check_position(x);

			switch(i)
			{
				case 0:	return (1.0-x[0]-x[1]);
				case 1:	return x[0];
				case 2:	return x[1];
				case 3: return 27.0*x[0]*x[1]*(1.0-x[0]-x[1]); // bubble
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
				case 3: g[0] = 27*x[1]*(1.0-x[1]-2.0*x[0]);
						g[1] = 27*x[0]*(1.0-x[0]-2.0*x[1]);
						std::cout << "MiniGrad: " << g[0] << ","<< g[1] << "at"<< x[0]<< ","<< x[1]<< std::endl;
						return;
				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};


////////////////////////////////////////////////////////////////////////////////
//	ReferenceQuadrilateral
//
//  function space: span {1,x,y,x^2-y^2} + 2 bubbles
//  ref: Wen Bai, CMAME 143 (1997), 41-47
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceQuadrilateral>
: public MiniBubbleLDS<ReferenceQuadrilateral>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceQuadrilateral>, 2>
{
	protected:
		static const double SQRT_FIVE;
		static const double SQRT_FIVTH;

	public:
	///	Order of Shape functions
		static constexpr size_t order = 2;

	///	Dimension, where shape functions are defined
		static constexpr int dim = 2;

	/// Number of shape functions
		static constexpr size_t nsh = 5;

	public:
	///	Shape type
		using shape_type = number;

	///	Gradient type
		using grad_type = MathVector<dim>;

	///	Reference Element type
		using reference_element_type = ReferenceQuadrilateral;

	public:
	///	Constructor
		MiniBubbleLSFS()= default;

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		[[nodiscard]] inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		[[nodiscard]] inline size_t num_sh() const override {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		[[nodiscard]] inline bool position(size_t i, MathVector<dim>& pos) const
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
				case 5: pos[0] = SQRT_FIVTH;
						pos[1] = SQRT_FIVTH; return true;

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		[[nodiscard]] inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			// ReferenceQuadrilateral::check_position(x);

			switch(i)
			{
				case 0:	return (1.0-x[0])*(1.0-x[1]);
				case 1:	return x[0]*(1.0-x[1]);
				case 2:	return x[0]*x[1];
				case 3:	return x[1]*(1.0-x[0]);

				// two bubbles (condensable)
				case 4:	return x[0]*(1.0-x[0])*x[1]*(1.0-x[1]);
				case 5:	return (1.0-x[0]*x[0])*(1.0-x[1]*x[1])*(x[0]+x[1])/(0.8*0.8*2.0)*SQRT_FIVTH; // max: sqrt(1/5) = sqrt(5)/5

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
				case 0:	g[0] = -(1.0-x[1]);
						g[1] = -(1.0-x[0]); return;
				case 1:	g[0] = (1.0-x[1]);
						g[1] = -(0.0+x[0]); return;
				case 2:	g[0] = (0.0+x[1]);
						g[1] = (0.0+x[0]); return;
				case 3: g[0] = -(0.0+x[1]);
						g[1] = (1.0-x[0]); return;
				// bubble 1
				case 4: g[0] = (1.0-2.0*x[0])*x[1]*(1.0-x[1]);
						g[1] = x[0]*(1.0-x[0])*(1.0-2.0*x[1]); return;
				// bubble 2, grad = 3*x^2*(y^2-1) + 2xy*(y^2-1) -(y^2-1) = (y^2-1)*(3*x^2+2xy-1)
				case 5: g[0] =(x[1]*x[1]-1.0)*(3.0*x[0]*x[0]+2.0*x[0]*x[1]-1.0)/(0.8*0.8*2.0)*SQRT_FIVTH;
						g[1] =(x[0]*x[0]-1.0)*(3.0*x[1]*x[1]+2.0*x[1]*x[0]-1.0)/(0.8*0.8*2.0)*SQRT_FIVTH; return;

				default: UG_THROW("MiniLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

////////////////////////////////////////////////////////////////////////////////
//	ReferenceTetrahedron
//
//  function space span {1,x,y,z} + bubble
//
////////////////////////////////////////////////////////////////////////////////

template <>
class MiniBubbleLSFS<ReferenceTetrahedron>
: public MiniBubbleLDS<ReferenceTetrahedron>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceTetrahedron>, 3>
{
	public:
	///	Order of Shape functions
		static constexpr size_t order = 1;

	///	Dimension, where shape functions are defined
		static constexpr int dim = 3;

	/// Number of shape functions
		static constexpr size_t nsh = 4+1;

	public:
	///	Shape type
	using shape_type = number;

	///	Gradient type
	using grad_type = MathVector<dim>;

	///	Reference Element type
	using reference_element_type = ReferenceTetrahedron;

	public:
	///	Constructor
		MiniBubbleLSFS()= default;

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		[[nodiscard]] inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		[[nodiscard]] inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		[[nodiscard]] inline bool position(size_t i, MathVector<dim>& pos) const
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
		[[nodiscard]] inline number shape(const size_t i, const MathVector<dim>& x) const
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
//  function space span {1,x} \times {1,y} \times {1,z} + bubble
//
////////////////////////////////////////////////////////////////////////////////
template <>
class MiniBubbleLSFS<ReferenceHexahedron>
: public MiniBubbleLDS<ReferenceHexahedron>,
  public BaseLSFS<MiniBubbleLSFS<ReferenceHexahedron>, 3>
{
	public:
	///	Order of Shape functions
		static constexpr size_t order = 6;

	///	Dimension, where shape functions are defined
		static constexpr int dim = 3;

	/// Number of shape functions
		static constexpr size_t nsh = 8+1;

	public:
	///	Shape type
	using shape_type = number;

	///	Gradient type
	using grad_type = MathVector<dim>;

	///	Reference Element type
	using reference_element_type = ReferenceHexahedron;

	public:
	///	Constructor
		MiniBubbleLSFS()= default;

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		[[nodiscard]] inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		[[nodiscard]] inline size_t num_sh() const override {return nsh;}


	///	\copydoc ug::LocalShapeFunctionSet::position()
		[[nodiscard]] inline bool position(size_t i, MathVector<dim>& pos) const
		{
			static const DimReferenceElement<3>& refElem = ReferenceElementProvider::get<3>(ROID_HEXAHEDRON);

			if (i<=7)
			{ 	// corner
				pos = refElem.corner(i); return true;
			} else if (i==8)
			{	// bubble
				pos[0] = 0.5; pos[1] = 0.5; pos[2] = 0.5; return true;
			} else
			{
				UG_THROW("MiniLSFS: shape function "<<i<< " not found. Only "<<nsh<<" shapes present.");
			}
			return false;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		[[nodiscard]] inline number shape(const size_t i, const MathVector<dim>& x) const
		{

			switch(i)
			{
			// corners 0..7
			case 0: return((1.0-x[0])*(1.0-x[1])*(1.0-x[2]));
			case 1: return((x[0])*(1.0-x[1])*(1.0-x[2]));
			case 2: return((x[0])*(x[1])*(1.0-x[2]));
			case 3: return((1.0-x[0])*(x[1])*(1.0-x[2]));
			case 4: return((1.0-x[0])*(1.0-x[1])*(x[2]));
			case 5: return((x[0])*(1.0-x[1])*(x[2]));
			case 6: return((x[0])*(x[1])*(x[2]));
			case 7: return((1.0-x[0])*(x[1])*(x[2]));
			// bubble
			case 8: return 64.0*x[0]*(1-x[0])*x[1]*(1-x[1])*x[2]*(1-x[2]);
			default: UG_THROW("MiniLSFS: shape function "<<i<<
					" not found. Only "<<nsh<<" shapes present.");
			}

		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathVector<dim>& value, const size_t i, const MathVector<dim>& x) const
		{

			switch(i)

			{
			case 0:
				value[0] = -(1.0-x[1])*(1.0-x[2]);
				value[1] = -(1.0-x[0])*(1.0-x[2]);
				value[2] = -(1.0-x[0])*(1.0-x[1]);
				return;
			case 1:
				value[0] = (1.0-x[1])*(1.0-x[2]);
				value[1] = -(x[0])*(1.0-x[2]);
				value[2] = -(x[0])*(1.0-x[1]);
				return;
			case 2:
				value[0] = (x[1])*(1.0-x[2]);
				value[1] = (x[0])*(1.0-x[2]);
				value[2] = -x[0]*x[1];
				return;
			case 3:
				value[0] = -(x[1])*(1.0-x[2]);
				value[1] = (1.0-x[0])*(1.0-x[2]);
				value[2] = -(1.0-x[0])*(x[1]);
				return;
			case 4:
				value[0] = -(1.0-x[1])*(x[2]);
				value[1] = -(1.0-x[0])*(x[2]);
				value[2] = (1.0-x[0])*(1.0-x[1]);
				return;
			case 5:
				value[0] = (1.0-x[1])*x[2];
				value[1] = -(x[0])*x[2];
				value[2] = (x[0])*(1.0-x[1]);
				return;
			case 6:
				value[0] = (x[1])*x[2];
				value[1] = (x[0])*x[2];
				value[2] = x[0]*x[1];
				return;
			case 7:
				value[0] = -(x[1])*x[2];
				value[1] = (1.0-x[0])*x[2];
				value[2] = (1.0-x[0])*x[1];
				return;
			case 8: // bubble
				value[0] = 64.0*(1.0-2.0*x[0])*x[1]*(1-x[1])*x[2]*(1-x[2]);
				value[1] = 64.0*(1.0-2.0*x[1])*x[0]*(1-x[0])*x[2]*(1-x[2]);
				value[2] = 64.0*(1.0-2.0*x[2])*x[0]*(1-x[0])*x[1]*(1-x[1]);
				return;

			default: UG_THROW("MiniLSFS: Invalid shape fct index: "<<i);
			}


		}
};

} //namespace ug

#endif