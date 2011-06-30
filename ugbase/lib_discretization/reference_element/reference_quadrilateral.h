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

	/// \copydoc ug::ReferenceElement::num_obj()
		size_t num_obj(int dim) const	{return m_vNum[dim];}

	/// \copydoc ug::ReferenceElement::num_obj_of_obj()
		size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const
			{return m_vSubNum[dim_i][i][dim_j];}

	/// \copydoc ug::ReferenceElement::id()
		int id(int dim_i, size_t i, int dim_j, size_t j) const
			{return m_id[dim_i][i][dim_j][j];}

	/// \copydoc ug::ReferenceElement::num_ref_elem()
		size_t num_ref_elem(ReferenceObjectID type) const {return m_vNumRefElem[type];}

	/// \copydoc ug::ReferenceElement::ref_elem_type()
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const{	return m_vRefElemType[dim_i][i];}

	/// \copydoc ug::DimReferenceElement::corner()
		const MathVector<dim>& corner(int i) const {return m_vCorner[i];}

	private:
	/// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2};
		enum{MAXOBJECTS = 4};

	/// number of Geometric Objects of Reference Element
	//  (m_num_obj[dim] = number of GeomObjects of dimension dim)
		size_t m_vNum[dim+1];

	/// number of Geometric Objects contained in a (Sub-)Geometric Object of the Element
		size_t m_vSubNum[dim+1][MAXOBJECTS][dim+1];

	///	coordinates of Reference Corner
		MathVector<dim> m_vCorner[num_corners];

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
		static const int world_dim = TWorldDim;
		static const int dim = ReferenceQuadrilateral::dim;

	public:
		ReferenceMapping() : m_corners(NULL)
		{}

		void update(const MathVector<world_dim>* corners)
		{
			m_corners = corners;
		}

		bool local_to_global(	const MathVector<dim>& loc_pos,
								MathVector<world_dim>& glob_pos) const
		{
			VecScaleAdd(glob_pos, 	(1.-loc_pos[0])*(1.-loc_pos[1]), m_corners[0],
									loc_pos[0]*(1.-loc_pos[1])     , m_corners[1],
									loc_pos[0]*loc_pos[1]          , m_corners[2],
									(1.-loc_pos[0])*loc_pos[1]     , m_corners[3]);
			return true;
		}

		bool jacobian_transposed(	const MathVector<dim>& loc_pos,
									MathMatrix<dim, world_dim>& JT) const
		{
			number a = 1. - loc_pos[1];

			JT(0, 0) = a*(m_corners[1][0] - m_corners[0][0]) + loc_pos[1]*(m_corners[2][0] - m_corners[3][0]);
			JT(0, 1) = a*(m_corners[1][1] - m_corners[0][1]) + loc_pos[1]*(m_corners[2][1] - m_corners[3][1]);

			a = 1. - loc_pos[0];
			JT(1, 0) = a*(m_corners[3][0] - m_corners[0][0]) + loc_pos[0]*(m_corners[2][0] - m_corners[1][0]);
			JT(1, 1) = a*(m_corners[3][1] - m_corners[0][1]) + loc_pos[0]*(m_corners[2][1] - m_corners[1][1]);
			return true;
		}

		bool jacobian_transposed_inverse(	const MathVector<dim>& loc_pos,
											MathMatrix<world_dim, dim>& JTInv) const
		{
			MathMatrix<dim, world_dim> JT;

			if(!jacobian_transposed(loc_pos, JT)) return false;

			// compute right inverse
			RightInverse(JTInv, JT);

			return true;
		}

		bool jacobian_det(const MathVector<dim>& loc_pos, number& det) const
		{
			MathMatrix<dim, world_dim> JT;

			if(!jacobian_transposed(loc_pos, JT)) return false;

			if( (world_dim == 2) && (dim==2) )
			{
				det = JT(0, 0)*JT(1, 1) - JT(0, 1)*JT(1, 0);
				return true;
			}
			//TODO: Implement pseudo inverse
			return false;
		}

	private:
		const MathVector<world_dim>* m_corners;

		MathVector<world_dim> a10, a20;

};

}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__ */
