/*
 * reference_prism.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_PRISM__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_PRISM__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

///	reference element for prism
class ReferencePrism : public DimReferenceElement<3>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PRISM;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int num_corners = 6;

	///	number of eges
		static const int num_edges = 9;

	///	number of faces
		static const int num_faces = 5;

	///	number of volumes
		static const int num_volumes = 1;

	public:
	///	Constructor
		ReferencePrism();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 0.5;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0 &&
					  pos[0]+pos[1] <= 1.0 &&
					  pos[2] >= 0.0 && pos[2] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

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

}

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_PRISM__ */
