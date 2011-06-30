/*
 * reference_prism.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_PRISM__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_PRISM__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

///	reference element for prism
class ReferencePrism{
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
		enum{MAXOBJECTS = 10};

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
class ReferenceMapping<ReferencePrism, TWorldDim>
{
	public:
		static const int world_dim = TWorldDim;
		static const int dim = ReferencePrism::dim;

	public:
		ReferenceMapping() : m_corners(NULL)
		{}

		void update(const MathVector<world_dim>* corners)
		{
			m_corners = corners;
		}

		bool local_to_global(	const MathVector<dim>& local,
								MathVector<world_dim>& global) const
		{
			number a,b, a0,a1,a2,a3,a4,a5;
			const MathVector<world_dim>* x = m_corners;

			a = 1.0 - (local)[0] - (local)[1];
			b = 1.0 - (local)[2];
			a0 = a * b;
			a1 = (local)[0] * b;
			a2 = (local)[1] * b;
			a3 = a * (local)[2];
			a4 = (local)[0] * (local)[2];
			a5 = (local)[1] * (local)[2];
			(global)[0] =
				a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+
				a4*(x)[4][0]+a5*(x)[5][0];
			(global)[1] =
					a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+
					a4*(x)[4][1]+a5*(x)[5][1];
			(global)[2] =
					a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+
					a4*(x)[4][2]+a5*(x)[5][2];
			return true;
		}

		bool jacobian_transposed(	const MathVector<dim>& local,
									MathMatrix<dim, world_dim>& JT) const
	   {
	        number a0,a1,a2,b0,b1,b2;
			const MathVector<world_dim>* x = m_corners;
	          a0 = (x)[0][0]-(x)[1][0]-(x)[3][0]+(x)[4][0];
	          a1 = (x)[0][1]-(x)[1][1]-(x)[3][1]+(x)[4][1];
	          a2 = (x)[0][2]-(x)[1][2]-(x)[3][2]+(x)[4][2];
	          b0 = (x)[0][0]-(x)[2][0]-(x)[3][0]+(x)[5][0];
	          b1 = (x)[0][1]-(x)[2][1]-(x)[3][1]+(x)[5][1];
	          b2 = (x)[0][2]-(x)[2][2]-(x)[3][2]+(x)[5][2];
	          JT(0,0) = (x)[1][0]-(x)[0][0]+(local)[2]*a0;
	          JT(0,1) = (x)[1][1]-(x)[0][1]+(local)[2]*a1;
	          JT(0,2) = (x)[1][2]-(x)[0][2]+(local)[2]*a2;
	          JT(1,0) = (x)[2][0]-(x)[0][0]+(local)[2]*b0;
	          JT(1,1) = (x)[2][1]-(x)[0][1]+(local)[2]*b1;
	          JT(1,2) = (x)[2][2]-(x)[0][2]+(local)[2]*b2;
	          JT(2,0) = (x)[3][0]-(x)[0][0]+(local)[0]*a0+(local)[1]*b0;
	          JT(2,1) = (x)[3][1]-(x)[0][1]+(local)[0]*a1+(local)[1]*b1;
	          JT(2,2) = (x)[3][2]-(x)[0][2]+(local)[0]*a2+(local)[1]*b2;
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
};

}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_PRISM__ */
