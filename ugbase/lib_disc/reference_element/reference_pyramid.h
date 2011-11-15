/*
 * reference_pyramid.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_PYRAMID__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_PYRAMID__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

///	reference element for a pyramid
class ReferencePyramid{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PYRAMID;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int num_corners = 5;

	///	number of eges
		static const int num_edges = 8;

	///	number of faces
		static const int num_faces = 5;

	///	number of volumes
		static const int num_volumes = 1;

	public:
	///	Constructor
		ReferencePyramid();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0/3.0;}

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
			//\todo: add check
		}

	private:
	/// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME= 3};
		enum{MAXOBJECTS = 8};

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

}

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_PYRAMID__ */
