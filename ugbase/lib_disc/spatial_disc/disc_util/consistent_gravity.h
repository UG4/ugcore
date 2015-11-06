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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONSISTENT_GRAVITY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONSISTENT_GRAVITY__

// other ug4 modules
#include "common/common.h"
#include <vector>

namespace ug{

/// Class for the computation of the standard version ('Voss-Souza-type') of the consistent gravity
/**
 * Density driven flow models involve the convection velocity depending on the
 * gradient of the pressure and the density-dependent gravity force, typically
 * of the form \f$\mathbf{v} = \mathbf{K} (- \nabla p + \rho \mathbf{g}) / \mu\f$,
 * where \f$\mathbf{K}\f$ is a matrix, \f$\mu\f$ a scalar, \f$\mathbf{g}\f$
 * gravity vector (all independent of \f$p\f$ or \f$\rho\f$), \f$p\f$ is
 * pressure, \f$rho\f$ is density (depending on the concentration, ...). (Both
 * the pressure and the concentration are unknown grid functions in the PDEs.)
 * Treating \f$p\f$ and \f$\rho\f$ as grid functions of the same class (like
 * the piecewise linear/bilinear interpolations of nodal values) leads to oscillatory
 * solutions because \f$\nabla p\f$ and \f$\rho \mathbf{g}\f$ belong to different
 * classes of grid functions. In particular, if \f$p\f$ is a piecewise linear
 * function, \f$\nabla p\f$ can never cancel the of the gravitaty force if 
 * \f$\rho\f$ increases linearly from the top down, although this happens in the
 * analytical solution. The idea of the consistent gravity is to consider a
 * vector function \f$\mathbf{h} = (h_x, h_y, \dots)\f$ that is in some sence a
 * primitive function for \f$\rho \mathbf{g}\f$ so that the velocity can be
 * written as \f$\mathbf{v} = - \mathbf{K} (p_x - h_x, p_y - h_y, \dots) / \mu\f$.
 * In the discretization, the function \f$\mathbf{h}\f$ should belong to the
 * same class of the grid functions as \f$p\f$. In the present implementation,
 * \f$\mathbf{h}\f$ is a piecewise linear/bilinear grid function.
 *
 * Class method 'prepare' computes the values of \f$h\f$ at the corners
 * of an element. Using these values, method 'compute' computes the
 * consistent gravity force \f$\rho \mathbf{g}\f$ (not \f$\mathbf{h}\f$!) at any
 * given point.
 *
 * \remark (Consistent gravity in the computation of the Jacobian)
 * The consistent gravity computed by the implemented methods depends linearly
 * on the density. Thus, to compute the derivative of the consistent gravity
 * w.r.t. the density at one of the corners of the element, set the density to
 * 1 at that corner and to 0 at all other corners, then call the functions. Then
 * 'prepare' prepares the nodal values of the derivatives and 'compute' computes
 * the derivative itself.
 *
 * Alternatively, you can set the density at that corner to the derivative
 * of the density w.r.t. the concentration (instead of 1). Then you get
 * directly the derivative of the consistent gravity w.r.t. the concentration
 * at that node.
 *
 * \remark There is an enhanced version of the consistent gravity. \see StdLinConsistentGravityX
 *
 * References:
 * <ul>
 *  <li> P. Frolkovic, P. Knabner, Consistent Velocity Approximations in Finite
 *       Element or Volume Discretizations of Density Driven Flow, In: Computational
 *       Methods in WaterResources XI, Vol. 1 (A.A. Aldama et al., eds.),
 *       Computational Mechanics Publication, Southhampten, 1996, p. 93-100
 *  </li>
 * </ul>
 *
 * \tparam refDim	dimensionality of the reference element (e.g. 2 for triangles, 3 for tetrahedra)
 */
template <int refDim>
class StdLinConsistentGravity
{
private:
//	static constants
	static const size_t _X_ = 0;
	static const size_t _Y_ = 1;
	static const size_t _Z_ = 2;
	
public:

///	constructor (sets the 'not init.' flag)
	StdLinConsistentGravity () : m_nCo (0) {};

///	computation of the primary function for the consistent gravity at corners, cf. the specializations
	template <int dim>
	inline void prepare
	(
		MathVector<refDim>* vConsGravity, ///< where to save the values (n_co vectors)
		const int n_co, ///< number of corners of the element
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
		const number* vDensity, ///< corner density (n_co scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		UG_THROW ("StdLinConsistentGravity: Combination of the world dim " << dim <<
			"and the reference element dim " << refDim << " is not implemented.");
	}
	
///	computation of the consistent gravity at a given point
	template <int dim>
	inline void compute
	(
		MathVector<dim>& ConsistentGravity, ///< where to save the vector
		const MathVector<refDim>& LocalCoord, ///< local coordinates of the point
		const MathMatrix<dim, refDim>& JTInv, ///< inverse transposed Jacobian
		const MathVector<refDim>* vLocalGrad, ///< gradients of the shape functions at the given point
		const MathVector<refDim>* vConsGravity ///< primary function of the consistent gravity at corners
	)
	{
		UG_ASSERT (m_nCo > 0, "StdLinConsistentGravity: Object not initialized.");
		
		MathVector<refDim> LocalGravity;
		VecSet(LocalGravity, 0.0);
		
	//	Loop shape functions
		for(size_t sh = 0; sh < (size_t) m_nCo; sh++)
			for(size_t d = 0; d < refDim; d++)
				LocalGravity[d] += vConsGravity[sh][d] * vLocalGrad[sh][d];

	//	Multiply by JacobianTransposedInverse
		MatVecMult(ConsistentGravity, JTInv, LocalGravity);
	}
	
protected:
	
	int m_nCo; ///< number of corners of the element for which the object is init. (0 if not init)

///	computation of the primary function for the consistent gravity at corners of an edge
	template <int dim>
	inline void prepare_edge
	(
		MathVector<1>* vConsGravity, ///< where to save the values (2 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (2 vectors)
		const number* vDensity, ///< corner density (2 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<1> LocalPoint;
		MathVector<1> LocalGravity;
		MathMatrix<1,dim> JT;
		static ReferenceMapping<ReferenceEdge, dim> EdgeMapping;
		EdgeMapping.update (vCorners);

		/* compute the local gravity */
		VecSet (LocalPoint, 0.0);
		EdgeMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravity, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravity[_X_]*(vDensity[0] + vDensity[1])*0.5;
	}
	
///	computation of the primary function for the consistent gravity at corners of a triangle
	template <int dim>
	inline void prepare_triangle
	(
		MathVector<2>* vConsGravity, ///< where to save the values (3 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (3 vectors)
		const number* vDensity, ///< corner density (3 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<2> LocalPoint;
		MathVector<2> LocalGravity;
		MathMatrix<2,dim> JT;
		static ReferenceMapping<ReferenceTriangle, dim> TriangleMapping;
		TriangleMapping.update (vCorners);

		/* compute the local gravity */
		VecSet (LocalPoint, 0.0);
		TriangleMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravity, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0; vConsGravity[2][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravity[_X_]*(vDensity[0] + vDensity[1])*0.5;

		vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
		vConsGravity[2][_Y_] = LocalGravity[_Y_]*(vDensity[0] + vDensity[2])*0.5;
	}

///	computation of the primary function for the consistent gravity at corners of a quadrilateral
	template <int dim>
	inline void prepare_quadrilateral
	(
		MathVector<2>* vConsGravity, ///< where to save the values (4 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (4 vectors)
		const number* vDensity, ///< corner density (4 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<2> LocalPoint;
		MathVector<2> LocalGravityAt000, LocalGravityAt110;
		MathMatrix<2,dim> JT;
		static ReferenceMapping<ReferenceQuadrilateral, dim> QuadMapping;
		QuadMapping.update (vCorners);

		/* compute the local gravity at local corner (0,0) */
		VecSet (LocalPoint, 0.0);
		QuadMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt000, JT, PhysicalGravity);

		/* compute the local gravity at local corner (1,1) */
		VecSet (LocalPoint, 1.0);
		QuadMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt110, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0; vConsGravity[3][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
		vConsGravity[2][_X_] = LocalGravityAt110[_X_]*(vDensity[2] + vDensity[3])*0.5;

		vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
		vConsGravity[2][_Y_] = LocalGravityAt110[_Y_]*(vDensity[1] + vDensity[2])*0.5;
		vConsGravity[3][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[3])*0.5;
	}
	
///	computation of the primary function for the consistent gravity at corners of a tetrahedron
	template <int dim>
	inline void prepare_tetrahedron
	(
		MathVector<3>* vConsGravity, ///< where to save the values (4 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (4 vectors)
		const number* vDensity, ///< corner density (4 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<3> LocalPoint;
		MathVector<3> LocalGravity;
		MathMatrix<3,dim> JT;
		static ReferenceMapping<ReferenceTetrahedron, dim> TetMapping;
		TetMapping.update (vCorners);

		/* compute the local gravity */
		LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
		TetMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravity, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0; vConsGravity[2][_X_] = 0.0; vConsGravity[3][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravity[_X_]*(vDensity[0] + vDensity[1])*0.5;

		vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0; vConsGravity[3][_Y_] = 0.0;
		vConsGravity[2][_Y_] = LocalGravity[_Y_]*(vDensity[0] + vDensity[2])*0.5;

		vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0; vConsGravity[2][_Z_] = 0.0;
		vConsGravity[3][_Z_] = LocalGravity[_Z_]*(vDensity[0] + vDensity[3])*0.5;
	}
	
///	computation of the primary function for the consistent gravity at corners of a pyramid
/**
 * TODO: Verify this implementation! Cf. UG3.
 */
	template <int dim>
	inline void prepare_pyramid
	(
		MathVector<3>* vConsGravity, ///< where to save the values (5 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (5 vectors)
		const number* vDensity, ///< corner density (5 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<3> LocalPoint;
		MathVector<3> LocalGravityAt000, LocalGravityAt110;
		MathMatrix<3,dim> JT;
		static ReferenceMapping<ReferencePyramid, dim> PyramidMapping;
		PyramidMapping.update (vCorners);

		/* compute the local gravity at local corner (0,0,0) */
		LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
		PyramidMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt000, JT, PhysicalGravity);

		/* compute the local gravity at local corner (1,1,0) */
		LocalPoint.x() = 1.0; LocalPoint.y() = 1.0; LocalPoint.z() = 0.0;
		PyramidMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt110, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0; vConsGravity[3][_X_] = 0.0; vConsGravity[4][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
		vConsGravity[2][_X_] = LocalGravityAt110[_X_]*(vDensity[2] + vDensity[3])*0.5;

		vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0; vConsGravity[4][_Y_] = 0.0;
		vConsGravity[2][_Y_] = LocalGravityAt110[_Y_]*(vDensity[1] + vDensity[2])*0.5;
		vConsGravity[3][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[3])*0.5;

		vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0;
		vConsGravity[2][_Z_] = 0.0; vConsGravity[3][_Z_] = 0.0;
		vConsGravity[4][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[4])*0.5;
	}
	
///	computation of the primary function for the consistent gravity at corners of a prism
	template <int dim>
	inline void prepare_prism
	(
		MathVector<3>* vConsGravity, ///< where to save the values (6 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (6 vectors)
		const number* vDensity, ///< corner density (6 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<3> LocalPoint;
		MathVector<3> LocalGravityAt000, LocalGravityAt101, LocalGravityAt011;
		MathMatrix<3,dim> JT;
		static ReferenceMapping<ReferencePrism, dim> PrismMapping;
		PrismMapping.update (vCorners);

		/* compute the local gravity at local corner (0,0,0) */
		LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
		PrismMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt000, JT, PhysicalGravity);

		/* compute the local gravity at local corner (1,0,1) */
		LocalPoint.x() = 1.0; LocalPoint.y() = 0.0; LocalPoint.z() = 1.0;
		PrismMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt101, JT, PhysicalGravity);

		/* compute the local gravity at local corner (0,1,1) */
		LocalPoint.x() = 0.0; LocalPoint.y() = 1.0; LocalPoint.z() = 1.0;
		PrismMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt011, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0; vConsGravity[2][_X_] = 0.0;
		vConsGravity[3][_X_] = 0.0; vConsGravity[5][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
		vConsGravity[4][_X_] = LocalGravityAt011[_X_]*(vDensity[3] + vDensity[4])*0.5;

		vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
		vConsGravity[3][_Y_] = 0.0; vConsGravity[4][_Y_] = 0.0;
		vConsGravity[2][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[2])*0.5;
		vConsGravity[5][_Y_] = LocalGravityAt011[_Y_]*(vDensity[3] + vDensity[5])*0.5;

		vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0; vConsGravity[2][_Z_] = 0.0;
		vConsGravity[3][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[3])*0.5;
		vConsGravity[4][_Z_] = LocalGravityAt101[_Z_]*(vDensity[1] + vDensity[4])*0.5;
		vConsGravity[5][_Z_] = LocalGravityAt011[_Z_]*(vDensity[2] + vDensity[5])*0.5;
	}
	
///	computation of the primary function for the consistent gravity at corners of a hexahedron
	template <int dim>
	inline void prepare_hexahedron
	(
		MathVector<3>* vConsGravity, ///< where to save the values (8 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (8 vectors)
		const number* vDensity, ///< corner density (8 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		MathVector<3> LocalPoint;
		MathVector<3> LocalGravityAt000, LocalGravityAt110, LocalGravityAt101, LocalGravityAt011;
		MathMatrix<3,dim> JT;
		static ReferenceMapping<ReferenceHexahedron, dim> HexMapping;
		HexMapping.update (vCorners);

		/* compute the local gravity at local corner (0,0,0) */
		LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
		HexMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt000, JT, PhysicalGravity);

		/* compute the local gravity at local corner (1,1,0) */
		LocalPoint.x() = 1.0; LocalPoint.y() = 1.0; LocalPoint.z() = 0.0;
		HexMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt110, JT, PhysicalGravity);

		/* compute the local gravity at local corner (1,0,1) */
		LocalPoint.x() = 1.0; LocalPoint.y() = 0.0; LocalPoint.z() = 1.0;
		HexMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt101, JT, PhysicalGravity);

		/* compute the local gravity at local corner (0,1,1) */
		LocalPoint.x() = 0.0; LocalPoint.y() = 1.0; LocalPoint.z() = 1.0;
		HexMapping.jacobian_transposed (JT, LocalPoint);
		MatVecMult (LocalGravityAt011, JT, PhysicalGravity);

		vConsGravity[0][_X_] = 0.0; vConsGravity[3][_X_] = 0.0;
		vConsGravity[4][_X_] = 0.0; vConsGravity[7][_X_] = 0.0;
		vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
		vConsGravity[2][_X_] = LocalGravityAt110[_X_]*(vDensity[2] + vDensity[3])*0.5;
		vConsGravity[5][_X_] = LocalGravityAt101[_X_]*(vDensity[4] + vDensity[5])*0.5;
		vConsGravity[6][_X_] = LocalGravityAt011[_X_]*(vDensity[6] + vDensity[7])*0.5;

		vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
		vConsGravity[4][_Y_] = 0.0; vConsGravity[5][_Y_] = 0.0;
		vConsGravity[2][_Y_] = LocalGravityAt110[_Y_]*(vDensity[1] + vDensity[2])*0.5;
		vConsGravity[3][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[3])*0.5;
		vConsGravity[6][_Y_] = LocalGravityAt101[_Y_]*(vDensity[5] + vDensity[6])*0.5;
		vConsGravity[7][_Y_] = LocalGravityAt011[_Y_]*(vDensity[4] + vDensity[7])*0.5;

		vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0;
		vConsGravity[2][_Z_] = 0.0; vConsGravity[3][_Z_] = 0.0;
		vConsGravity[4][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[4])*0.5;
		vConsGravity[5][_Z_] = LocalGravityAt101[_Z_]*(vDensity[1] + vDensity[5])*0.5;
		vConsGravity[6][_Z_] = LocalGravityAt110[_Z_]*(vDensity[2] + vDensity[6])*0.5;
		vConsGravity[7][_Z_] = LocalGravityAt011[_Z_]*(vDensity[3] + vDensity[7])*0.5;
	}
};

/// spacialization of the method for edges (reference dimension 1)
template <>
template <int dim>
void StdLinConsistentGravity<1>::prepare
(
	MathVector<1>* vConsGravity, ///< where to save the values (n_co vectors)
	const int n_co, ///< number of corners of the element (should be 2)
	const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
	const number* vDensity, ///< corner density (n_co scalar values)
	const MathVector<dim>& PhysicalGravity ///< the gravity vector
)
{
	UG_ASSERT (n_co == 2, "StdLinConsistentGravity: Illegal number of corners of an edge.");
	m_nCo = n_co;
	this->template prepare_edge<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
}

/// spacialization of the method for faces (reference dimension 2)
template <>
template <int dim>
void StdLinConsistentGravity<2>::prepare
(
	MathVector<2>* vConsGravity, ///< where to save the values (n_co vectors)
	const int n_co, ///< number of corners of the element
	const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
	const number* vDensity, ///< corner density (n_co scalar values)
	const MathVector<dim>& PhysicalGravity ///< the gravity vector
)
{
	switch (m_nCo = n_co)
	{
		case 3:
			this->template prepare_triangle<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
			break;
		case 4:
			this->template prepare_quadrilateral<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
			break;
		default:
			UG_THROW ("StdLinConsistentGravity: Illegal number of corners ("
				<< n_co << ") of an element with reference dimension 2.");
	}
}

/// spacialization of the method for volumes (reference dimension 3)
template <>
template <int dim>
void StdLinConsistentGravity<3>::prepare
(
	MathVector<3>* vConsGravity, ///< where to save the values (n_co vectors)
	const int n_co, ///< number of corners of the element
	const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
	const number* vDensity, ///< corner density (n_co scalar values)
	const MathVector<dim>& PhysicalGravity ///< the gravity vector
)
{
	switch (m_nCo = n_co)
	{
		case 4:
			this->template prepare_tetrahedron<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
			break;
		case 5:
			this->template prepare_pyramid<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
			break;
		case 6:
			this->template prepare_prism<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
			break;
		case 8:
			this->template prepare_hexahedron<dim> (vConsGravity, vCorners, vDensity, PhysicalGravity);
			break;
		default:
			UG_THROW ("StdLinConsistentGravity: Illegal number of corners ("
				<< n_co << ") of an element with reference dimension 3.");
	}
}

/// Class for the computation of the enhanced version ('Frolkovic-type') of the consistent gravity
/**
 * Implementation of the enhanced ('Frolkovic-type') version of the consistent
 * gravity for simplices (triangles and tetrahedra) in the case of the gravity
 * force parallel to the z-axis. For all other elements and other gravities,
 * the same method as in StdLinConsistentGravity is used.
 * \see StdLinConsistentGravity
 *
 * References:
 * <ul>
 *  <li> P. Frolkovic, Consistent velocity approximation for density driven
 *       flow and transport, In: Advanced Computational Methods in Engineering,
 *       Part 2: Contributed papers; (R. Van Keer at al., eds.), Shaker Publishing,
 *       Maastricht, 1998, p. 603-611
 *  </li>
 * </ul>
 *
 * \tparam refDim	dimensionality of the reference element (e.g. 2 for triangles, 3 for tetrahedra)
 */
template <int refDim>
class StdLinConsistentGravityX : public StdLinConsistentGravity<refDim>
{
	typedef StdLinConsistentGravity<refDim> base_type;
	
public:

///	constructor
	StdLinConsistentGravityX () {}
	
///	computation of the primary function for the consistent gravity at corners, cf. the specializations
	template <int dim>
	inline void prepare
	(
		MathVector<refDim>* vConsGravity, ///< where to save the values (n_co vectors)
		const int n_co, ///< number of corners of the element
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
		const number* vDensity, ///< corner density (n_co scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector
	)
	{
		base_type::template prepare<dim>
			(vConsGravity, n_co, vCorners, vDensity, PhysicalGravity);
		// C.f. the specializations below for the enhanced version
	}
	
///	computation of the consistent gravity at a given point
	template <int dim>
	inline void compute
	(
		MathVector<dim>& ConsistentGravity, ///< where to save the vector
		const MathVector<refDim>& LocalCoord, ///< local coordinates of the point
		const MathMatrix<dim, refDim>& JTInv, ///< inverse transposed Jacobian
		const MathVector<refDim>* vLocalGrad, ///< gradients of the shape functions at the given point
		const MathVector<refDim>* vConsGravity ///< primary function of the consistent gravity at corners
	)
	{
		if (base_type::m_nCo > 0) // use the standard version
			base_type::template compute<dim> (ConsistentGravity, LocalCoord, JTInv, vLocalGrad, vConsGravity);
		else
		{
		//	special method for triangles and tetrahedra, currently only in the full dimension
			UG_ASSERT (dim == refDim && base_type::m_nCo == -(refDim+1), "StdLinConsistentGravityX: Illegal initialization of the object.");
		
			MathVector<refDim> LocalGravity;
			VecSet (LocalGravity, 0.0);
			
			for (size_t d = 0; d < refDim; d++)
				LocalGravity[d] = vConsGravity[d+1][dim-1];
	
			//	multiply by JacobianTransposedInverse
			MatVecMult(ConsistentGravity, JTInv, LocalGravity);
			
			//	correct the first coordinates
			for (size_t d = 0; d < dim-1; d++)
				for (size_t sh = 1; sh <= refDim; sh++)
					ConsistentGravity[d] -= vConsGravity[sh][d] * LocalCoord[sh-1];
		}
	}
	
protected:
	
///	computation of the extended version of the corner consistent gravity for simplices (only in full dimension!)
	template <typename refElem, int dim>
	inline void prepare_simplex
	(
		MathVector<refDim>* vConsGravity, ///< where to save the values (dim+1 vectors)
		const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (dim+1 vectors)
		const number* vDensity, ///< corner density (dim+1 scalar values)
		const MathVector<dim>& PhysicalGravity ///< the gravity vector (MUST BE (0, ..., 0, g))
	)
	{
		number DensityIP, Diff, gradient;
		MathVector<refDim> LocalPoint;
		MathVector<dim> ShiftedGlobalPoint;
		MathMatrix<refDim,dim> JT;
		static ReferenceMapping<refElem, dim> Mapping;
		Mapping.update (vCorners);
		Mapping.jacobian_transposed_inverse (JT, vCorners[0]); // the 2. argument is dummy

		for (size_t i = 1; i < dim+1; i++)
		{
			VecSubtract (ShiftedGlobalPoint, vCorners[i], vCorners[0]);
			Diff = ShiftedGlobalPoint [dim-1];

			ShiftedGlobalPoint [dim-1] = 0.0;
			TransposedMatVecMult (LocalPoint, JT, ShiftedGlobalPoint);
			
			DensityIP = 1;
			for (size_t j = 0; j < dim; j++) DensityIP -= LocalPoint[j]; //GN(0)=1-local(0)-local(1)-...
			DensityIP *= vDensity[0];
			for (size_t j = 1; j < dim+1; j++)
				DensityIP += vDensity[j] * LocalPoint[j-1]; //GN(j)=local(j-1) for j=1,2,...
			
			for(size_t j = 0; j < dim-1; ++j)
			{
				gradient = 0.0;
				for (size_t k = 0; k < dim; k++)
					gradient += JT[j][k] * (vDensity[k+1] - vDensity[0]);
				vConsGravity[i][j] = Diff * PhysicalGravity[dim-1] * gradient;
			}
			vConsGravity[i][dim-1] = Diff * PhysicalGravity[dim-1] * (vDensity[i] + DensityIP)*0.5;
		}
	}
};

/// spacialization of the method for faces (reference dimension 2)
template <>
template <int dim>
void StdLinConsistentGravityX<2>::prepare
(
	MathVector<2>* vConsGravity, ///< where to save the values (n_co vectors)
	const int n_co, ///< number of corners of the element
	const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
	const number* vDensity, ///< corner density (n_co scalar values)
	const MathVector<dim>& PhysicalGravity ///< the gravity vector
)
{
	if (dim == 2 && n_co == 3 && PhysicalGravity[0] == 0.0) // use the enhanced version
	{
		m_nCo = -3;
		this->template prepare_simplex<ReferenceTriangle, dim>
			(vConsGravity, vCorners, vDensity, PhysicalGravity);
	}
	else // use the standard version
		base_type::template prepare<dim>
			(vConsGravity, n_co, vCorners, vDensity, PhysicalGravity);
}

/// spacialization of the method for volumes (reference dimension 3)
template <>
template <int dim>
void StdLinConsistentGravityX<3>::prepare
(
	MathVector<3>* vConsGravity, ///< where to save the values (n_co vectors)
	const int n_co, ///< number of corners of the element
	const MathVector<dim>* vCorners, ///< (global) coordinates of the corners (n_co vectors)
	const number* vDensity, ///< corner density (n_co scalar values)
	const MathVector<dim>& PhysicalGravity ///< the gravity vector
)
{
	if (dim == 3 && n_co == 4 && PhysicalGravity[0] == 0.0 && PhysicalGravity[1] == 0.0) // use the enhanced version
	{
		m_nCo = -4;
		this->template prepare_simplex<ReferenceTetrahedron, dim>
			(vConsGravity, vCorners, vDensity, PhysicalGravity);
	}
	else // use the standard version
		base_type::template prepare<dim>
			(vConsGravity, n_co, vCorners, vDensity, PhysicalGravity);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONSISTENT_GRAVITY__ */
