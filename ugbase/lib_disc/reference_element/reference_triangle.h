/*
 * reference_triangle.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_TRIANGLE__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_TRIANGLE__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

class ReferenceTriangle : public DimReferenceElement<2>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;

	///	dimension of reference element
		static const int dim = 2;

	///	number of corners
		static const int num_corners = 3;

	///	number of eges
		static const int num_edges = 3;

	///	number of faces
		static const int num_faces = 1;

	///	number of volumes
		static const int num_volumes = 0;

	public:
	///	Constructor filling the arrays
		ReferenceTriangle();

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
			          pos[0]+pos[1] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};


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

}

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_TRIANGLE__ */
