/*
 * reference_quadrilateral.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

class ReferenceQuadrilateral{
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

	/// \copydoc ug::ReferenceElement::num(int)
		size_t num(int dim) const	{return m_vNum[dim];}

	/// \copydoc ug::ReferenceElement::num(int, size_t, int)
		size_t num(int dim_i, size_t i, int dim_j) const
			{return m_vSubNum[dim_i][i][dim_j];}

	/// \copydoc ug::ReferenceElement::id()
		int id(int dim_i, size_t i, int dim_j, size_t j) const
			{return m_id[dim_i][i][dim_j][j];}

	/// \copydoc ug::ReferenceElement::num_ref_elem()
		size_t num_ref_elem(ReferenceObjectID type) const {return m_vNumRefElem[type];}

	/// \copydoc ug::ReferenceElement::ref_elem_type()
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const{	return m_vRefElemType[dim_i][i];}

	/// \copydoc ug::DimReferenceElement::corner()
		const MathVector<dim>& corner(size_t i) const {return m_vCorner[i];}

	/// \copydoc ug::DimReferenceElement::corner()
		const MathVector<dim,int>* corner() const {return m_vCoInt;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}

	private:
	/// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2};
		enum{MAXOBJECTS = 4};

	/// number of Geometric Objects of a dimension
	
		size_t m_vNum[dim+1];

	/// number of Geometric Objects contained in a (Sub-)Geometric Object of the Element
		size_t m_vSubNum[dim+1][MAXOBJECTS][dim+1];

	///	coordinates of Reference Corner
		MathVector<dim> m_vCorner[num_corners];
		MathVector<dim, int> m_vCoInt[num_corners];

	/// indices of GeomObjects
		int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

	///	number of reference elements
		size_t m_vNumRefElem[NUM_REFERENCE_OBJECTS];

	///	type of reference elements
		ReferenceObjectID m_vRefElemType[dim+1][MAXOBJECTS];
};

template <>
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

		void local_to_global(	const MathVector<dim>& locPos,
								MathVector<worldDim>& globPos) const
		{
			VecScaleAdd(globPos, 	(1.-locPos[0])*(1.-locPos[1]), m_vCo[0],
									locPos[0]*(1.-locPos[1])     , m_vCo[1],
									locPos[0]*locPos[1]          , m_vCo[2],
									(1.-locPos[0])*locPos[1]     , m_vCo[3]);
		}

		void jacobian_transposed(	const MathVector<dim>& locPos,
									MathMatrix<dim, worldDim>& JT) const
		{
			number a = 1. - locPos[1];

			JT(0, 0) = a*(m_vCo[1][0] - m_vCo[0][0]) + locPos[1]*(m_vCo[2][0] - m_vCo[3][0]);
			JT(0, 1) = a*(m_vCo[1][1] - m_vCo[0][1]) + locPos[1]*(m_vCo[2][1] - m_vCo[3][1]);

			a = 1. - locPos[0];
			JT(1, 0) = a*(m_vCo[3][0] - m_vCo[0][0]) + locPos[0]*(m_vCo[2][0] - m_vCo[1][0]);
			JT(1, 1) = a*(m_vCo[3][1] - m_vCo[0][1]) + locPos[0]*(m_vCo[2][1] - m_vCo[1][1]);
		}

		void jacobian_transposed_inverse(	const MathVector<dim>& locPos,
											MathMatrix<worldDim, dim>& JTInv) const
		{
			MathMatrix<dim, worldDim> JT;

			jacobian_transposed(locPos, JT);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

		number jacobian_det(const MathVector<dim>& locPos) const
		{
			MathMatrix<dim, worldDim> JT;

			jacobian_transposed(locPos, JT);

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

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__ */
