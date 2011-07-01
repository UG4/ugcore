/*
 * reference_edge.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_EDGE__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_EDGE__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

///	reference element for an edge
class ReferenceEdge
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;

	///	dimension of reference element
		static const int dim = 1;

	///	number of corners
		static const int num_corners = 2;

	///	number of eges
		static const int num_edges = 1;

	///	number of faces
		static const int num_faces = 0;

	///	number of volumes
		static const int num_volumes = 0;

	public:
	///	Constructor
		ReferenceEdge();

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

	/// \copydoc ug::DimReferenceElement::corner()
		const MathVector<dim,int>* corner() const {return m_vCoInt;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0,
			          "Local position "<<pos<<" outside Reference Element");
		}

	private:
	/// to make it more readable
		enum{POINT = 0, EDGE = 1};
		enum{MAXOBJECTS = 2};

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
class ReferenceMapping<ReferenceEdge, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceEdge::dim;

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

		void local_to_global(	const MathVector<dim>& locPos,
								MathVector<worldDim>& globPos) const
		{
			globPos = m_vCo[0];
			VecScaleAppend(globPos, locPos[0], a10);
		}

		void jacobian_transposed(	const MathVector<dim>& locPos,
									MathMatrix<dim, worldDim>& JT) const
		{
			for(int i = 0; i < worldDim; ++i) JT(0,i) = a10[i];
		}

		void jacobian_transposed_inverse(	const MathVector<dim>& locPos,
											MathMatrix<worldDim, dim>& JTInv) const
		{
			MathMatrix<dim, worldDim> JT;

		// 	get jacobian transposed
			jacobian_transposed(locPos, JT);

		// 	compute right inverse
			RightInverse(JTInv, JT);
		}

		number jacobian_det(const MathVector<dim>& locPos) const
		{
			return a10[0];
		}

	private:
		const MathVector<worldDim>* m_vCo;

		MathVector<worldDim> a10;
};


} // namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_EDGE__ */
