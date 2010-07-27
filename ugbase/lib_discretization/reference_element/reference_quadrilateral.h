/*
 * reference_quadrilateral.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__

namespace ug{


class ReferenceQuadrilateral{
	public:
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
		static const int dim = 2;
		static const int num_corners = 4;
		static const int num_edges = 4;
		static const int num_faces = 1;
		static const int num_volumes = 0;

	public:
		ReferenceQuadrilateral(){initializeArrays();}

		/// reference object id
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

		/// Dimension where reference element lives
		int dimension() const {return dim;}

		/// size of reference triangle
		number size() const	{return 0.5;}

		/// number of objects of dim
		size_t num_obj(int dim) const	{return m_num_obj[dim];}

		/// number of object of dim
		size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const
			{return m_num_obj_of_obj[dim_i][i][dim_j];}

		/// id of object j in dimension dim_j of obj i in dimension dim_i
		int id(int dim_i, size_t i, int dim_j, size_t j) const
			{return m_id[dim_i][i][dim_j][j];}

		/// number of reference elements this element is contained of
		size_t num_ref_elem(ReferenceObjectID type) const {return m_ref_elem[type];}

		/// reference element type of obj nr i in dimension dim_i */
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const{	return m_ref_elem_type[dim_i][i];}

		/// coordinates of reference corner (i = 0 ... num_obj(0))
		const MathVector<dim>& corner(int i) const {return m_corner[i];}

	private:
		// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2};
		enum{MAXOBJECTS = 4};

		/* number of Geometric Objects of Reference Element
		 * (m_num_obj[dim] = number of GeomObjects of dimension dim) */
		size_t m_num_obj[dim+1];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		size_t m_num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
		/* coordinates of Reference Corner */
		MathVector<dim> m_corner[num_corners];
		// indices of GeomObjects
		int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

		size_t m_ref_elem[NUM_REFERENCE_OBJECTS];
		ReferenceObjectID m_ref_elem_type[dim+1][MAXOBJECTS];

		void initializeArrays()
		{
			//number of Geometric Objects
		 	m_num_obj[POINT] = 4;
		 	m_num_obj[EDGE] = 4;
		 	m_num_obj[FACE] = 1;

			// number of Geometric Objects
		 	m_num_obj_of_obj[FACE][0][POINT] = 4;
		 	m_num_obj_of_obj[FACE][0][EDGE] = 4;
		 	m_num_obj_of_obj[FACE][0][FACE] = 1;
		 	m_ref_elem_type[FACE][0] = ROID_QUADRILATERAL;

		 	for(size_t i = 0; i < m_num_obj[EDGE]; ++i)
		 	{
		 		m_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	m_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	m_num_obj_of_obj[EDGE][i][FACE] = 1;

			 	m_ref_elem_type[EDGE][i] = ROID_EDGE;
		 	}

		 	for(size_t i = 0; i < m_num_obj[EDGE]; ++i)
		 	{
		 		m_num_obj_of_obj[POINT][i][POINT] = 1;
		 		m_num_obj_of_obj[POINT][i][EDGE] = 2;
		 		m_num_obj_of_obj[POINT][i][FACE] = 1;

			 	m_ref_elem_type[POINT][i] = ROID_VERTEX;
		 	}

			//reset m_id to -1
			for(int i=0; i<=dim; ++i)
				for(size_t j=0; j<MAXOBJECTS; ++j)
					for(int k=0; k<=dim; ++k)
						for(size_t l=0; l<MAXOBJECTS; l++)
						{
						 	m_id[i][j][k][l] = -1;
						}

			//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
			for(int d=0; d<=dim; ++d)
				for(size_t j=0; j<m_num_obj[d]; ++j)
				{
				 	m_id[d][j][d][0] = j;
				}

			//Edges <-> Face
			for(size_t i=0; i<m_num_obj[EDGE]; ++i)
			{
			 	m_id[FACE][0][EDGE][i] = i;
			 	m_id[EDGE][i][FACE][0] = 0;
			}

			// Points <-> Face
			for(size_t i=0; i<m_num_obj[POINT]; ++i)
			{
			 	m_id[FACE][0][POINT][i] = i;
			 	m_id[POINT][i][FACE][0] = 0;
			}

			// Points of Edges
			// edge 0 = (0,1)
		 	m_id[EDGE][0][POINT][0] = 0;
		 	m_id[EDGE][0][POINT][1] = 1;
			// edge 1 = (1,2)
		 	m_id[EDGE][1][POINT][0] = 1;
		 	m_id[EDGE][1][POINT][1] = 2;
			// edge 2 = (2,3)
		 	m_id[EDGE][2][POINT][0] = 2;
		 	m_id[EDGE][2][POINT][1] = 3;
			// edge 3 = (3,0)
		 	m_id[EDGE][3][POINT][0] = 3;
		 	m_id[EDGE][3][POINT][1] = 0;

		 	// Edges of Point
		 	m_id[POINT][0][EDGE][0] = 3;
		 	m_id[POINT][0][EDGE][1] = 0;

		 	m_id[POINT][1][EDGE][0] = 0;
		 	m_id[POINT][1][EDGE][1] = 1;

		 	m_id[POINT][2][EDGE][0] = 1;
		 	m_id[POINT][2][EDGE][1] = 2;

		 	m_id[POINT][3][EDGE][0] = 2;
		 	m_id[POINT][3][EDGE][1] = 3;

			// Reference Corners
		 	m_corner[0] = MathVector<dim>(0.0, 0.0);
		 	m_corner[1] = MathVector<dim>(1.0, 0.0);
		 	m_corner[2] = MathVector<dim>(1.0, 1.0);
		 	m_corner[3] = MathVector<dim>(0.0, 1.0);

		 	// Reference Element Types
		 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		 	{
				m_ref_elem[i] = 0;
		 	}
		 	m_ref_elem[ROID_VERTEX] = 4;
		 	m_ref_elem[ROID_EDGE] = 4;
		 	m_ref_elem[ROID_QUADRILATERAL] = 1;
		}

};

template <>
template <int TWorldDim>
class ReferenceMapping<ReferenceQuadrilateral, TWorldDim>
{
	public:
		static const int world_dim = TWorldDim;
		static const int dim = ReferenceTriangle::dim;

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

			JT[0][0] = a*(m_corners[1][0] - m_corners[0][0]) + loc_pos[1]*(m_corners[2][0] - m_corners[3][0]);
			JT[0][1] = a*(m_corners[1][1] - m_corners[0][1]) + loc_pos[1]*(m_corners[2][1] - m_corners[3][1]);

			a = 1. - loc_pos[0];
			JT[1][0] = a*(m_corners[3][0] - m_corners[0][0]) + loc_pos[0]*(m_corners[2][0] - m_corners[1][0]);
			JT[1][1] = a*(m_corners[3][1] - m_corners[0][1]) + loc_pos[0]*(m_corners[2][1] - m_corners[1][1]);
			return true;
		}

		bool jacobian_transposed_inverse(	const MathVector<dim>& loc_pos,
											MathMatrix<world_dim, dim>& JTInv) const
		{
			MathMatrix<dim, world_dim> JT;

			if(!jacobian_transposed(loc_pos, JT)) return false;

			if( (world_dim == 2) && (dim==2) )
			{
				const number det = JT[0][0]*JT[1][1] - JT[0][1]*JT[1][0];
				UG_ASSERT(det != 0.0, "Zero Determinant. Impossible to invert");

				JTInv[0][0] = JT[1][1] / det;
				JTInv[1][0] = -JT[1][0] / det;
				JTInv[0][1] = -JT[0][1] / det;
				JTInv[1][1] = JT[0][0] / det;
				return true;
			}

			//TODO: Implement pseudo inverse
			return false;
		}

		bool jacobian_det(const MathVector<dim>& loc_pos, number& det) const
		{
			MathMatrix<dim, world_dim> JT;

			if(!jacobian_transposed(loc_pos, JT)) return false;

			if( (world_dim == 2) && (dim==2) )
			{
				det = JT[0][0]*JT[1][1] - JT[0][1]*JT[1][0];
				return true;
			}
			//TODO: Implement pseudo inverse
			return false;
		}

	private:
		const MathVector<world_dim>* m_corners;

		MathVector<world_dim> a10, a20;

};


template <>
class reference_element_traits<Quadrilateral>
{
	public:
		typedef ReferenceQuadrilateral reference_element_type;
};



}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_QUADRILATERAL__ */
