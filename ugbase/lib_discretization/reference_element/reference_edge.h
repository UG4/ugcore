/*
 * reference_edge.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_EDGE__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_EDGE__

namespace ug{

class ReferenceEdge
{
	public:
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
		static const int dim = 1;
		static const int num_corners = 2;
		static const int num_edges = 1;
		static const int num_faces = 0;
		static const int num_volumes = 0;

	public:
		ReferenceEdge() {initializeArrays();}

		/// reference object id
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

		/// Dimension where reference element lives
		int dimension() const {return dim;}

		/// size of reference triangle
		number size() const	{return 0.5;}

		/// number of objects of dim
		size_t num_obj(int dim)	const {return m_num_obj[dim];}

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
		enum{POINT = 0, EDGE = 1};
		enum{MAXOBJECTS = 2};

		/* number of Geometric Objects of Reference Element
		 * (m_num_obj[dim] = number of GeomObjects of dimension dim) */
		size_t m_num_obj[dim+1];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		size_t m_num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
		/* coordinates of Reference Corner */
		MathVector<dim> m_corner[2];
		// indices of GeomObjects
		int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

		size_t m_ref_elem[NUM_REFERENCE_OBJECTS];
		ReferenceObjectID m_ref_elem_type[dim+1][MAXOBJECTS];

		void initializeArrays()
		{
			//number of Geometric Objects
			m_num_obj[POINT] = 2;
			m_num_obj[EDGE] = 1;

			// number of Geometric Objects
			m_num_obj_of_obj[EDGE][0][POINT] = 2;
			m_num_obj_of_obj[EDGE][0][EDGE] = 1;

			m_ref_elem_type[EDGE][0] = ROID_EDGE;

			for(size_t i = 0; i < m_num_obj[POINT]; ++i)
			{
				m_num_obj_of_obj[POINT][i][POINT] = 1;
				m_num_obj_of_obj[POINT][i][EDGE] = 1;

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
			for(int i=0; i<=dim; ++i)
				for(size_t j=0; j<m_num_obj[i]; ++j)
				{
					m_id[i][j][i][0] = j;
				}

			// Points <-> Face
			for(size_t i=0; i<m_num_obj[POINT]; ++i)
			{
				m_id[EDGE][0][POINT][i] = i;
				m_id[POINT][i][EDGE][0] = 0;
			}

			// Reference Corners
			m_corner[0][0] = 0.0;
			m_corner[1][0] = 1.0;

			// Reference Element Types
			for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
			{
				m_ref_elem[i] = 0;
			}
			m_ref_elem[ROID_VERTEX] = 2;
			m_ref_elem[ROID_EDGE] = 1;
		}
};


template <>
template <int TWorldDim>
class ReferenceMapping<ReferenceEdge, TWorldDim>
{
	public:
		static const int world_dim = TWorldDim;
		static const int dim = ReferenceEdge::dim;

	public:
		ReferenceMapping() : m_corners(NULL)
		{}

		void update(const MathVector<world_dim>* corners)
		{
			m_corners = corners;
			VecSubtract(a10, m_corners[1], m_corners[0]);
		}

		bool local_to_global(	const MathVector<dim>& loc_pos,
								MathVector<world_dim>& glob_pos) const
		{
			glob_pos = m_corners[0];
			VecScaleAppend(glob_pos, loc_pos[0], a10);
			return true;
		}

		bool jacobian_transposed(	const MathVector<dim>& loc_pos,
									MathMatrix<dim, world_dim>& JT) const
		{
			for(int i = 0; i < world_dim; ++i)
			{
				JT(0,i) = a10[i];
			}
			return true;
		}

		bool jacobian_transposed_inverse(	const MathVector<dim>& loc_pos,
											MathMatrix<world_dim, dim>& JTInv) const
		{
			MathMatrix<dim, world_dim> JT;

			// get jacobian transposed
			if(!jacobian_transposed(loc_pos, JT))
				{UG_LOG("Cannot get jacobian transposed.\n");return false;}

			// compute right inverse
			RightInverse(JTInv, JT);

			return true;
		}

		bool jacobian_det(const MathVector<dim>& loc_pos, number& det) const
		{
			det = a10[0];
			return true;
		}

	private:
		const MathVector<world_dim>* m_corners;

		MathVector<world_dim> a10;

};


template <>
class reference_element_traits<Edge>
{
	public:
		typedef ReferenceEdge reference_element_type;
};




} // namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_EDGE__ */
