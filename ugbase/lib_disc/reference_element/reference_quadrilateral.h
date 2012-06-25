/*
 * reference_quadrilateral.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

class ReferenceQuadrilateral : public DimReferenceElement<2>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;

	///	dimension of reference element
		static const int dim = 2;

	///	number of corners
		static const int num_corners = 4;

	///	number of eges
		static const int num_edges = 4;

	///	number of faces
		static const int num_faces = 1;

	///	number of volumes
		static const int num_volumes = 0;

	public:
		ReferenceQuadrilateral();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

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

}

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__ */
