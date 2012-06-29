/*
 * reference_element_mapping.h
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__

#include <cassert>
#include <iostream>
#include <sstream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

/**
 * This class describes the mapping from a reference element into the real
 * (physical) world. The mapping is initialized by the physical positions of
 * the vertices of the real world element. The order of those points must be
 * given as indicated by the corresponding reference element.
 *
 * Let \f$R\f$ be the reference element and \f$T\f$ be the element. Then, the
 * reference mapping is a mapping:
 * \f[
 * 	\phi:	R \mapsto T
 * \f]
 *
 * \tparam	TRefElem		reference element
 * \tparam	TWorldDim		world dimension
 */
template <typename TRefElem, int TWorldDim>
class ReferenceMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefElem::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping();

	///	Constructor setting the corners of the element
		ReferenceMapping(const MathVector<worldDim>* vCorner);

	///	refresh mapping for new set of corners
		void update(const MathVector<worldDim>* vCorner);

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
		                     const MathVector<dim> locPos) const;

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                         const MathVector<dim> locPos) const;

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                 const MathVector<dim> locPos) const;

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const;
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Concrete Reference Mappings
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Vertex
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Edge
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferenceEdge, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceEdge::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = true;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCorner)
		{
			m_vCo = vCorner;
			VecSubtract(a10, m_vCo[1], m_vCo[0]);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{
			globPos = m_vCo[0];
			VecScaleAppend(globPos, locPos[0], a10);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
		{
			for(int i = 0; i < worldDim; ++i) JT(0,i) = a10[i];
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			return a10[0];
		}

	private:
		const MathVector<worldDim>* m_vCo;

		MathVector<worldDim> a10;
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Triangle
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferenceTriangle, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceTriangle::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCorner)
		{
			m_vCo = vCorner;
			VecSubtract(a10, m_vCo[1], m_vCo[0]);
			VecSubtract(a20, m_vCo[2], m_vCo[0]);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{
			globPos = m_vCo[0];
			VecScaleAppend(globPos, locPos[0], a10);
			VecScaleAppend(globPos, locPos[1], a20);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
		{
			for(int i = 0; i < worldDim; ++i)
			{
				JT(0, i) = a10[i];
				JT(1, i) = a20[i];
			}
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			return a10[0]*a20[1] - a10[1]*a20[0];
		}

	private:
		const MathVector<worldDim>* m_vCo;

		MathVector<worldDim> a10, a20;

};
///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Quadrilateral
///////////////////////////////////////////////////////////////////////////////


template <int TWorldDim>
class ReferenceMapping<ReferenceQuadrilateral, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceQuadrilateral::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCorner)
		{
			m_vCo = vCorner;
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{
			VecScaleAdd(globPos, 	(1.-locPos[0])*(1.-locPos[1]), m_vCo[0],
									locPos[0]*(1.-locPos[1])     , m_vCo[1],
									locPos[0]*locPos[1]          , m_vCo[2],
									(1.-locPos[0])*locPos[1]     , m_vCo[3]);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
		{
			const number a = 1. - locPos[1];
			const number b = 1. - locPos[0];

			for(int i = 0; i < worldDim; ++i)
			{
				JT(0, i) = a*(m_vCo[1][i] - m_vCo[0][i]) + locPos[1]*(m_vCo[2][i] - m_vCo[3][i]);
				JT(1, i) = b*(m_vCo[3][i] - m_vCo[0][i]) + locPos[0]*(m_vCo[2][i] - m_vCo[1][i]);
			}
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			MathMatrix<dim, worldDim> JT;

			jacobian_transposed(JT, locPos);

			if( (worldDim == 2) && (dim==2) )
			{
				const number det = JT(0, 0)*JT(1, 1) - JT(0, 1)*JT(1, 0);
				return det;
			}

			//TODO: Implement pseudo inverse
			UG_ASSERT(0, "Not implemented");
			return 0.0;
		}

	private:
		const MathVector<worldDim>* m_vCo;

		MathVector<worldDim> a10, a20;

};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferenceTetrahedron, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceTetrahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = true;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornr)
		{
			m_vCo = vCornr;
			VecSubtract(a10, m_vCo[1], m_vCo[0]);
			VecSubtract(a20, m_vCo[2], m_vCo[0]);
			VecSubtract(a30, m_vCo[3], m_vCo[0]);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{
			globPos = m_vCo[0];
			VecScaleAppend(globPos, locPos[0], a10);
			VecScaleAppend(globPos, locPos[1], a20);
			VecScaleAppend(globPos, locPos[2], a30);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
		{
			for(int i = 0; i < worldDim; ++i)
			{
				JT[0][i] = a10[i];
				JT[1][i] = a20[i];
				JT[2][i] = a30[i];
			}
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			MathMatrix<dim, worldDim> JT;

		//	compute jacobian transposed
			jacobian_transposed(JT, locPos);

		//	only in quad case defined
			if((dim==3) && (worldDim==3))
			{
				const number det =
				JT(0,0)*JT(1,1)*JT(2,2)
				+ JT(0,1)*JT(1,2)*JT(2,0)
				+ JT(0,2)*JT(1,0)*JT(2,1)
				- JT(0,0)*JT(1,2)*JT(2,1)
				- JT(0,1)*JT(1,0)*JT(2,2)
				- JT(0,2)*JT(1,1)*JT(2,0);
				return det;
			}

			UG_ASSERT(0, "Not implemented.");
			return 0.0;
		}

	private:
		const MathVector<worldDim>* m_vCo;

		MathVector<worldDim> a10, a20, a30;
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Pyramid
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferencePyramid, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferencePyramid::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCorner)
		{
			m_vCo = vCorner;
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{
			number a,b,a0,a1,a2,a3;
			const MathVector<worldDim>* x = m_vCo;
			a = 1.0 - locPos[0];
			b = 1.0 - locPos[1];
			if (locPos[0] > locPos[1])
			{
				a0 = a * b - locPos[2] * b;
				a1 = locPos[0] * b - locPos[2]*locPos[1];
				a2 = locPos[0] * locPos[1] + locPos[2]*locPos[1];
				a3 = a * locPos[1] - locPos[2] * locPos[1];
				globPos[0] =
						a0*x[0][0]+a1*x[1][0]+a2*x[2][0]+a3*x[3][0]+locPos[2]*x[4][0];
				globPos[1] =
						a0*x[0][1]+a1*x[1][1]+a2*x[2][1]+a3*x[3][1]+locPos[2]*x[4][1];
				globPos[2] =
						a0*x[0][2]+a1*x[1][2]+a2*x[2][2]+a3*x[3][2]+locPos[2]*x[4][2];
			}
			else
			{
				a0 = a * b - locPos[2] * a;
				a1 = locPos[0] * b - locPos[2]*locPos[0];
				a2 = locPos[0] * locPos[1] + locPos[2]*locPos[0];
				a3 = a * locPos[1] - locPos[2] * locPos[0];
				globPos[0] =
						a0*x[0][0]+a1*x[1][0]+a2*x[2][0]+a3*x[3][0]+locPos[2]*x[4][0];
				globPos[1] =
						a0*x[0][1]+a1*x[1][1]+a2*x[2][1]+a3*x[3][1]+locPos[2]*x[4][1];
				globPos[2] =
						a0*x[0][2]+a1*x[1][2]+a2*x[2][2]+a3*x[3][2]+locPos[2]*x[4][2];
			}
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
	   {
			number a,b,c;
			const MathVector<worldDim>* x = m_vCo;

			a = x[0][0]-x[1][0]+x[2][0]-x[3][0];
			b = x[0][1]-x[1][1]+x[2][1]-x[3][1];
			c = x[0][2]-x[1][2]+x[2][2]-x[3][2];

			if (locPos[0] > locPos[1])
			{
				JT(0,0) = x[1][0]-x[0][0]+locPos[1]*a;
				JT(0,1) = x[1][1]-x[0][1]+locPos[1]*b;
				JT(0,2) = x[1][2]-x[0][2]+locPos[1]*c;
				JT(1,0) = x[3][0]-x[0][0]+(locPos[0]+locPos[2])*a;
				JT(1,1) = x[3][1]-x[0][1]+(locPos[0]+locPos[2])*b;
				JT(1,2) = x[3][2]-x[0][2]+(locPos[0]+locPos[2])*c;
				JT(2,0) = x[4][0]-x[0][0]+locPos[1]*a;
				JT(2,1) = x[4][1]-x[0][1]+locPos[1]*b;
				JT(2,2) = x[4][2]-x[0][2]+locPos[1]*c;
			}
			else
			{
				JT(0,0) = x[1][0]-x[0][0]+(locPos[1]+locPos[2])*a;
				JT(0,1) = x[1][1]-x[0][1]+(locPos[1]+locPos[2])*b;
				JT(0,2) = x[1][2]-x[0][2]+(locPos[1]+locPos[2])*c;
				JT(1,0) = x[3][0]-x[0][0]+locPos[0]*a;
				JT(1,1) = x[3][1]-x[0][1]+locPos[0]*b;
				JT(1,2) = x[3][2]-x[0][2]+locPos[0]*c;
				JT(2,0) = x[4][0]-x[0][0]+locPos[0]*a;
				JT(2,1) = x[4][1]-x[0][1]+locPos[0]*b;
				JT(2,2) = x[4][2]-x[0][2]+locPos[0]*c;
			}
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			MathMatrix<dim, worldDim> JT;
			jacobian_transposed(JT, locPos);
			if((dim==3) && (worldDim==3))
			{
				const number det
				= JT(0,0)*JT(1,1)*JT(2,2)
				+ JT(0,1)*JT(1,2)*JT(2,0)
				+ JT(0,2)*JT(1,0)*JT(2,1)
				- JT(0,0)*JT(1,2)*JT(2,1)
				- JT(0,1)*JT(1,0)*JT(2,2)
				- JT(0,2)*JT(1,1)*JT(2,0);
				return det;
			}

			UG_ASSERT(0, "Not implemented");
			return 0.0;
		}

	private:
		const MathVector<worldDim>* m_vCo;
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Prism
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferencePrism, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferencePrism::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCorner)
		{
			m_vCo = vCorner;
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{
			number a,b, a0,a1,a2,a3,a4,a5;
			const MathVector<worldDim>* x = m_vCo;

			a = 1.0 - locPos[0] - locPos[1];
			b = 1.0 - locPos[2];
			a0 = a * b;
			a1 = locPos[0] * b;
			a2 = locPos[1] * b;
			a3 = a * locPos[2];
			a4 = locPos[0] * locPos[2];
			a5 = locPos[1] * locPos[2];
			globPos[0] =
				a0*x[0][0]+a1*x[1][0]+a2*x[2][0]+a3*x[3][0]+a4*x[4][0]+a5*x[5][0];
			globPos[1] =
				a0*x[0][1]+a1*x[1][1]+a2*x[2][1]+a3*x[3][1]+a4*x[4][1]+a5*x[5][1];
			globPos[2] =
				a0*x[0][2]+a1*x[1][2]+a2*x[2][2]+a3*x[3][2]+a4*x[4][2]+a5*x[5][2];
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
	   {
	        number a0,a1,a2,b0,b1,b2;
			const MathVector<worldDim>* x = m_vCo;
	          a0 = x[0][0]-x[1][0]-x[3][0]+x[4][0];
	          a1 = x[0][1]-x[1][1]-x[3][1]+x[4][1];
	          a2 = x[0][2]-x[1][2]-x[3][2]+x[4][2];
	          b0 = x[0][0]-x[2][0]-x[3][0]+x[5][0];
	          b1 = x[0][1]-x[2][1]-x[3][1]+x[5][1];
	          b2 = x[0][2]-x[2][2]-x[3][2]+x[5][2];
	          JT(0,0) = x[1][0]-x[0][0]+locPos[2]*a0;
	          JT(0,1) = x[1][1]-x[0][1]+locPos[2]*a1;
	          JT(0,2) = x[1][2]-x[0][2]+locPos[2]*a2;
	          JT(1,0) = x[2][0]-x[0][0]+locPos[2]*b0;
	          JT(1,1) = x[2][1]-x[0][1]+locPos[2]*b1;
	          JT(1,2) = x[2][2]-x[0][2]+locPos[2]*b2;
	          JT(2,0) = x[3][0]-x[0][0]+locPos[0]*a0+locPos[1]*b0;
	          JT(2,1) = x[3][1]-x[0][1]+locPos[0]*a1+locPos[1]*b1;
	          JT(2,2) = x[3][2]-x[0][2]+locPos[0]*a2+locPos[1]*b2;
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			MathMatrix<dim, worldDim> JT;
			jacobian_transposed(JT, locPos);
			if((dim==3) && (worldDim==3))
			{
				const number det
				= JT(0,0)*JT(1,1)*JT(2,2)
				+ JT(0,1)*JT(1,2)*JT(2,0)
				+ JT(0,2)*JT(1,0)*JT(2,1)
				- JT(0,0)*JT(1,2)*JT(2,1)
				- JT(0,1)*JT(1,0)*JT(2,2)
				- JT(0,2)*JT(1,1)*JT(2,0);
				return det;
			}

			UG_ASSERT(0, "Not implemented");
			return 0.0;
		}

	private:
		const MathVector<worldDim>* m_vCo;
};


///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferenceHexahedron, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceHexahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping() : m_vCo(NULL) {}

	///	Constructor setting the corners
		ReferenceMapping(const MathVector<worldDim>* vCorner) {update(vCorner);}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCorner)
		{
			m_vCo = vCorner;
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim> locPos) const
		{

			number a,b,c,a0,a1,a2,a3,a4,a5,a6,a7;
			const MathVector<worldDim>* x = m_vCo;
			a = 1.0 - locPos[0];
			b = 1.0 - locPos[1];
			c = 1.0 - locPos[2];
			a0 = a * b * c;
			a1 = locPos[0] * b * c;
			a2 = locPos[0] * locPos[1] * c;
			a3 = a * locPos[1] * c;
			a4 = a * b * locPos[2];
			a5 = locPos[0] * b * locPos[2];
			a6 = locPos[0] * locPos[1] * locPos[2];
			a7 = a * locPos[1] * locPos[2];
			globPos[0] =
			        a0*x[0][0]+a1*x[1][0]+a2*x[2][0]+a3*x[3][0]+
			        a4*x[4][0]+a5*x[5][0]+a6*x[6][0]+a7*x[7][0];
			globPos[1] =
			        a0*x[0][1]+a1*x[1][1]+a2*x[2][1]+a3*x[3][1]+
			        a4*x[4][1]+a5*x[5][1]+a6*x[6][1]+a7*x[7][1];
			globPos[2] =
			        a0*x[0][2]+a1*x[1][2]+a2*x[2][2]+a3*x[3][2]+
			        a4*x[4][2]+a5*x[5][2]+a6*x[6][2]+a7*x[7][2];
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim> locPos) const
	   {
			number a,b,c,a0,a1,a2,a3;
			const MathVector<worldDim>* x = m_vCo;
			a = 1.0 - locPos[0];
			b = 1.0 - locPos[1];
			c = 1.0 - locPos[2];
			a0 = b * c;
			a1 = locPos[1] * c;
			a2 = locPos[1] * locPos[2];
			a3 = b * locPos[2];
			JT(0,0) = a0*(x[1][0]-x[0][0])+a1*(x[2][0]-x[3][0])
						+ a2*(x[6][0]-x[7][0])+a3*(x[5][0]-x[4][0]);
	        JT(0,1) = a0*(x[1][1]-x[0][1])+a1*(x[2][1]-x[3][1])
	                    + a2*(x[6][1]-x[7][1])+a3*(x[5][1]-x[4][1]);
	        JT(0,2) = a0*(x[1][2]-x[0][2])+a1*(x[2][2]-x[3][2])
	                    + a2*(x[6][2]-x[7][2])+a3*(x[5][2]-x[4][2]);
	        a0 = a * c;
	        a1 = locPos[0] * c;
	        a2 = locPos[0] * locPos[2];
	        a3 = a * locPos[2];
	        JT(1,0) = a0*(x[3][0]-x[0][0])+a1*(x[2][0]-x[1][0])
	                    + a2*(x[6][0]-x[5][0])+a3*(x[7][0]-x[4][0]);
	        JT(1,1) = a0*(x[3][1]-x[0][1])+a1*(x[2][1]-x[1][1])
	                    + a2*(x[6][1]-x[5][1])+a3*(x[7][1]-x[4][1]);
	        JT(1,2) = a0*(x[3][2]-x[0][2])+a1*(x[2][2]-x[1][2])
	                    + a2*(x[6][2]-x[5][2])+a3*(x[7][2]-x[4][2]);
	        a0 = a * b;
	        a1 = locPos[0] * b;
	        a2 = locPos[0] * locPos[1];
	        a3 = a * locPos[1];
	        JT(2,0) = a0*(x[4][0]-x[0][0])+a1*(x[5][0]-x[1][0])
	                    + a2*(x[6][0]-x[2][0])+a3*(x[7][0]-x[3][0]);
	        JT(2,1) = a0*(x[4][1]-x[0][1])+a1*(x[5][1]-x[1][1])
	                    + a2*(x[6][1]-x[2][1])+a3*(x[7][1]-x[3][1]);
	        JT(2,2) = a0*(x[4][2]-x[0][2])+a1*(x[5][2]-x[1][2])
	                    + a2*(x[6][2]-x[2][2])+a3*(x[7][2]-x[3][2]);
		}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim> locPos) const
		{
		//	temporary matrix for jacobian transposed
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(JT, locPos);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const
		{
			MathMatrix<dim, worldDim> JT;
			jacobian_transposed(JT, locPos);
			if((dim==3) && (worldDim==3))
			{
				const number det
				= JT(0,0)*JT(1,1)*JT(2,2)
				+ JT(0,1)*JT(1,2)*JT(2,0)
				+ JT(0,2)*JT(1,0)*JT(2,1)
				- JT(0,0)*JT(1,2)*JT(2,1)
				- JT(0,1)*JT(1,0)*JT(2,2)
				- JT(0,2)*JT(1,1)*JT(2,0);
				return det;
			}

			UG_ASSERT(0, "Not implemented");
			return 1.0;
		}

	private:
		const MathVector<worldDim>* m_vCo;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__ */
