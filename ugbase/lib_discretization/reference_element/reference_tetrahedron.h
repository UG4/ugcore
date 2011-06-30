/*
 * reference_tetrahedron.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

class ReferenceTetrahedron{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TETRAHEDRON;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int num_corners = 4;

	///	number of eges
		static const int num_edges = 6;

	///	number of faces
		static const int num_faces = 4;

	///	number of volumes
		static const int num_volumes = 1;

	public:
	///	Constructor
		ReferenceTetrahedron();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0/6.0;}

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

	private:
	/// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME= 3};
		enum{MAXOBJECTS = 6};

	/// number of Geometric Objects of a dimension
	
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
class ReferenceMapping<ReferenceTetrahedron, TWorldDim>
{
	public:
		static const int world_dim = TWorldDim;
		static const int dim = ReferenceTetrahedron::dim;

	public:
		ReferenceMapping() : m_corners(NULL)
		{}

		void update(const MathVector<world_dim>* corners)
		{
			m_corners = corners;
			VecSubtract(a10, m_corners[1], m_corners[0]);
			VecSubtract(a20, m_corners[2], m_corners[0]);
			VecSubtract(a30, m_corners[3], m_corners[0]);
		}

		bool local_to_global(	const MathVector<dim>& loc_pos,
								MathVector<world_dim>& glob_pos) const
		{
			glob_pos = m_corners[0];
			VecScaleAppend(glob_pos, loc_pos[0], a10);
			VecScaleAppend(glob_pos, loc_pos[1], a20);
			VecScaleAppend(glob_pos, loc_pos[1], a30);
			return true;
		}

		bool jacobian_transposed(	const MathVector<dim>& loc_pos,
									MathMatrix<dim, world_dim>& JT) const
		{
			for(int i = 0; i < world_dim; ++i)
			{
				JT[0][i] = a10[i];
				JT[1][i] = a20[i];
				JT[2][i] = a30[i];
			}
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
			if((dim==3) && (world_dim==3))
			{
				det = JT(0,0)*JT(1,1)*JT(2,2)
				+ JT(0,1)*JT(1,2)*JT(2,0)
				+ JT(0,2)*JT(1,0)*JT(2,1)
				- JT(0,0)*JT(1,2)*JT(2,1)
				- JT(0,1)*JT(1,0)*JT(2,2)
				- JT(0,2)*JT(1,1)*JT(2,0);
				return true;
			}
			return false;
		}

	private:
		const MathVector<world_dim>* m_corners;

		MathVector<world_dim> a10, a20, a30;
};


}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__ */
