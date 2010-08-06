/*
 * reference_pyramid.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_PYRAMID__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_PYRAMID__

namespace ug{


class ReferencePyramid{
	public:
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PYRAMID;
		static const int dim = 3;
		static const int num_corners = 5;
		static const int num_edges = 8;
		static const int num_faces = 5;
		static const int num_volumes = 1;

	public:
		ReferencePyramid(){initializeArrays();}

		/// reference object id
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

		/// Dimension where reference element lives
		int dimension() const {return dim;}

		/// size of reference triangle
		number size() const	{return 1.0/3.0;}

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
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME= 3};
		enum{MAXOBJECTS = 8};

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
		 	m_num_obj[POINT] = 5;
		 	m_num_obj[EDGE] = 8;
		 	m_num_obj[FACE] = 5;
		 	m_num_obj[VOLUME] = 1;

			// number of Geometric Objects
		 	m_num_obj_of_obj[VOLUME][0][POINT] = 5;
		 	m_num_obj_of_obj[VOLUME][0][EDGE] = 8;
		 	m_num_obj_of_obj[VOLUME][0][FACE] = 5;
		 	m_num_obj_of_obj[VOLUME][0][VOLUME] = 1;

			m_num_obj_of_obj[FACE][0][POINT] = 4;
			m_num_obj_of_obj[FACE][0][EDGE] = 4;
			m_num_obj_of_obj[FACE][0][FACE] = 1;
			m_num_obj_of_obj[FACE][0][VOLUME] = 1;
			m_ref_elem_type[FACE][0] = ROID_QUADRILATERAL;

			for(size_t i = 1; i < m_num_obj[FACE]; ++i)
		 	{
				m_num_obj_of_obj[FACE][i][POINT] = 3;
				m_num_obj_of_obj[FACE][i][EDGE] = 3;
				m_num_obj_of_obj[FACE][i][FACE] = 1;
				m_num_obj_of_obj[FACE][i][VOLUME] = 1;

				m_ref_elem_type[FACE][i] = ROID_TRIANGLE;
		 	}

		 	for(size_t i = 0; i < m_num_obj[EDGE]; ++i)
		 	{
			 	m_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	m_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	m_num_obj_of_obj[EDGE][i][FACE] = 2;
			 	m_num_obj_of_obj[EDGE][i][VOLUME] = 1;

			 	m_ref_elem_type[EDGE][i] = ROID_EDGE;
		 	}

		 	for(size_t i = 0; i < m_num_obj[POINT] - 1; ++i)
		 	{
			 	m_num_obj_of_obj[POINT][i][POINT] = 1;
			 	m_num_obj_of_obj[POINT][i][EDGE] = 3;
			 	m_num_obj_of_obj[POINT][i][FACE] = 3;
			 	m_num_obj_of_obj[POINT][i][VOLUME] = 1;

			 	m_ref_elem_type[POINT][i] = ROID_VERTEX;
		 	}
		 	m_num_obj_of_obj[POINT][4][POINT] = 1;
		 	m_num_obj_of_obj[POINT][4][EDGE] = 4;
		 	m_num_obj_of_obj[POINT][4][FACE] = 4;
		 	m_num_obj_of_obj[POINT][4][VOLUME] = 1;

		 	m_ref_elem_type[POINT][4] = ROID_VERTEX;

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

			// Face <-> Volume
			for(size_t i=0; i<m_num_obj[VOLUME]; ++i)
			{
			 	m_id[VOLUME][0][FACE][i] = i;
			 	m_id[FACE][i][VOLUME][0] = 0;
			}

			// Edge <-> Volume
			for(size_t i=0; i<m_num_obj[EDGE]; ++i)
			{
			 	m_id[VOLUME][0][EDGE][i] = i;
			 	m_id[EDGE][i][VOLUME][0] = 0;
			}

			// Point <-> Volume
			for(size_t i=0; i<m_num_obj[POINT]; ++i)
			{
			 	m_id[VOLUME][0][POINT][i] = i;
			 	m_id[POINT][i][VOLUME][0] = 0;
			}

			// Points <-> Faces
		 	m_id[FACE][0][POINT][0] = 0;
		 	m_id[FACE][0][POINT][1] = 3;
		 	m_id[FACE][0][POINT][2] = 2;
		 	m_id[FACE][0][POINT][3] = 1;

		 	m_id[FACE][1][POINT][0] = 0;
		 	m_id[FACE][1][POINT][1] = 1;
		 	m_id[FACE][1][POINT][2] = 4;

		 	m_id[FACE][2][POINT][0] = 1;
		 	m_id[FACE][2][POINT][1] = 2;
		 	m_id[FACE][2][POINT][2] = 4;

		 	m_id[FACE][3][POINT][0] = 2;
		 	m_id[FACE][3][POINT][1] = 3;
		 	m_id[FACE][3][POINT][2] = 4;

		 	m_id[FACE][4][POINT][0] = 3;
		 	m_id[FACE][4][POINT][1] = 0;
		 	m_id[FACE][4][POINT][2] = 4;

		 	m_id[POINT][0][FACE][0] = 0;
		 	m_id[POINT][0][FACE][1] = 1;
		 	m_id[POINT][0][FACE][2] = 4;

		 	m_id[POINT][1][FACE][0] = 0;
		 	m_id[POINT][1][FACE][1] = 1;
		 	m_id[POINT][1][FACE][2] = 2;

		 	m_id[POINT][2][FACE][0] = 0;
		 	m_id[POINT][2][FACE][1] = 2;
		 	m_id[POINT][2][FACE][2] = 3;

		 	m_id[POINT][3][FACE][0] = 0;
		 	m_id[POINT][3][FACE][1] = 3;
		 	m_id[POINT][3][FACE][2] = 4;

		 	m_id[POINT][4][FACE][0] = 1;
		 	m_id[POINT][4][FACE][1] = 2;
		 	m_id[POINT][4][FACE][2] = 3;
		 	m_id[POINT][4][FACE][3] = 4;

		 	// Edges <-> Faces
		 	m_id[FACE][0][EDGE][0] = 3;
		 	m_id[FACE][0][EDGE][1] = 2;
		 	m_id[FACE][0][EDGE][2] = 1;
		 	m_id[FACE][0][EDGE][3] = 0;

		 	m_id[FACE][1][EDGE][0] = 0;
		 	m_id[FACE][1][EDGE][1] = 5;
		 	m_id[FACE][1][EDGE][2] = 4;

		 	m_id[FACE][2][EDGE][0] = 5;
		 	m_id[FACE][2][EDGE][1] = 6;
		 	m_id[FACE][2][EDGE][2] = 1;

		 	m_id[FACE][3][EDGE][0] = 2;
		 	m_id[FACE][3][EDGE][1] = 7;
		 	m_id[FACE][3][EDGE][2] = 6;

		 	m_id[FACE][4][EDGE][0] = 3;
		 	m_id[FACE][4][EDGE][1] = 4;
		 	m_id[FACE][4][EDGE][2] = 7;

		 	m_id[EDGE][0][FACE][0] = 0;
		 	m_id[EDGE][0][FACE][1] = 1;

		 	m_id[EDGE][1][FACE][0] = 0;
		 	m_id[EDGE][1][FACE][1] = 2;

		 	m_id[EDGE][2][FACE][0] = 0;
		 	m_id[EDGE][2][FACE][1] = 3;

		 	m_id[EDGE][3][FACE][0] = 0;
		 	m_id[EDGE][3][FACE][1] = 4;

		 	m_id[EDGE][4][FACE][0] = 1;
		 	m_id[EDGE][4][FACE][1] = 4;

		 	m_id[EDGE][5][FACE][0] = 2;
		 	m_id[EDGE][5][FACE][1] = 1;

		 	m_id[EDGE][6][FACE][0] = 3;
		 	m_id[EDGE][6][FACE][1] = 2;

		 	m_id[EDGE][7][FACE][0] = 4;
		 	m_id[EDGE][7][FACE][1] = 3;

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
			// edge 4 = (0,4)
		 	m_id[EDGE][4][POINT][0] = 0;
		 	m_id[EDGE][4][POINT][1] = 4;
			// edge 5 = (1,4)
		 	m_id[EDGE][5][POINT][0] = 1;
		 	m_id[EDGE][5][POINT][1] = 4;
			// edge 6 = (2,4)
		 	m_id[EDGE][6][POINT][0] = 2;
		 	m_id[EDGE][6][POINT][1] = 4;
			// edge 7 = (3,4)
		 	m_id[EDGE][7][POINT][0] = 3;
		 	m_id[EDGE][7][POINT][1] = 4;

		 	// Edges of Point
		 	m_id[POINT][0][EDGE][0] = 3;
		 	m_id[POINT][0][EDGE][1] = 0;
		 	m_id[POINT][0][EDGE][2] = 4;

		 	m_id[POINT][1][EDGE][0] = 0;
		 	m_id[POINT][1][EDGE][1] = 5;
		 	m_id[POINT][1][EDGE][2] = 1;

		 	m_id[POINT][2][EDGE][0] = 1;
		 	m_id[POINT][2][EDGE][1] = 6;
		 	m_id[POINT][2][EDGE][2] = 2;

		 	m_id[POINT][3][EDGE][0] = 2;
		 	m_id[POINT][3][EDGE][1] = 7;
		 	m_id[POINT][3][EDGE][2] = 3;

		 	m_id[POINT][4][EDGE][0] = 4;
		 	m_id[POINT][4][EDGE][1] = 5;
		 	m_id[POINT][4][EDGE][2] = 6;
		 	m_id[POINT][4][EDGE][3] = 7;

		 	// Reference Corners
		 	m_corner[0] = MathVector<dim>(0.0, 0.0, 0.0);
		 	m_corner[1] = MathVector<dim>(1.0, 0.0, 0.0);
		 	m_corner[2] = MathVector<dim>(1.0, 1.0, 0.0);
		 	m_corner[3] = MathVector<dim>(0.0, 1.0, 0.0);
		 	m_corner[4] = MathVector<dim>(0.0, 0.0, 1.0);

		 	// Reference Element Types
		 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		 	{
				m_ref_elem[i] = 0;
		 	}
		 	m_ref_elem[ROID_VERTEX] = 5;
		 	m_ref_elem[ROID_EDGE] = 18;
		 	m_ref_elem[ROID_QUADRILATERAL] = 1;
		 	m_ref_elem[ROID_TRIANGLE] = 4;
		 	m_ref_elem[ROID_PYRAMID] = 1;
		}


};

template <>
template <int TWorldDim>
class ReferenceMapping<ReferencePyramid, TWorldDim>
{
	public:
		static const int world_dim = TWorldDim;
		static const int dim = ReferencePyramid::dim;

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
			number a,b,a0,a1,a2,a3;
			const MathVector<world_dim>* x = m_corners;
			a = 1.0 - (local)[0];
			b = 1.0 - (local)[1];
			if ((local)[0] > (local)[1]) {
			a0 = a * b - (local)[2] * b;
			a1 = (local)[0] * b - (local)[2]*(local)[1];
			a2 = (local)[0] * (local)[1] + (local)[2]*(local)[1];
			a3 = a * (local)[1] - (local)[2] * (local)[1];
			(global)[0] =
			        a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+(local)[2]*(x)[4][0];
			(global)[1] =
			        a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+(local)[2]*(x)[4][1];
			(global)[2] =
			        a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+(local)[2]*(x)[4][2];}
			else {
			a0 = a * b - (local)[2] * a;
			a1 = (local)[0] * b - (local)[2]*(local)[0];
			a2 = (local)[0] * (local)[1] + (local)[2]*(local)[0];
			a3 = a * (local)[1] - (local)[2] * (local)[0];
			(global)[0] =
			        a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+(local)[2]*(x)[4][0];
			(global)[1] =
			        a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+(local)[2]*(x)[4][1];
			(global)[2] =
			        a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+(local)[2]*(x)[4][2];}

			return true;
		}

		bool jacobian_transposed(	const MathVector<dim>& local,
									MathMatrix<dim, world_dim>& JT) const
	   {
			number a,b,c;
			const MathVector<world_dim>* x = m_corners;
			a = (x)[0][0]-(x)[1][0]+(x)[2][0]-(x)[3][0];
			b = (x)[0][1]-(x)[1][1]+(x)[2][1]-(x)[3][1];
			c = (x)[0][2]-(x)[1][2]+(x)[2][2]-(x)[3][2];
			if ((local)[0] > (local)[1]) {
			JT(0,0) = (x)[1][0]-(x)[0][0]+(local)[1]*a;
			JT(0,1) = (x)[1][1]-(x)[0][1]+(local)[1]*b;
			JT(0,2) = (x)[1][2]-(x)[0][2]+(local)[1]*c;
			JT(1,0) = (x)[3][0]-(x)[0][0]+((local)[0]+(local)[2])*a;
			JT(1,1) = (x)[3][1]-(x)[0][1]+((local)[0]+(local)[2])*b;
			JT(1,2) = (x)[3][2]-(x)[0][2]+((local)[0]+(local)[2])*c;
			JT(2,0) = (x)[4][0]-(x)[0][0]+(local)[1]*a;
			JT(2,1) = (x)[4][1]-(x)[0][1]+(local)[1]*b;
			JT(2,2) = (x)[4][2]-(x)[0][2]+(local)[1]*c;}
			else {
			JT(0,0) = (x)[1][0]-(x)[0][0]+((local)[1]+(local)[2])*a;
			JT(0,1) = (x)[1][1]-(x)[0][1]+((local)[1]+(local)[2])*b;
			JT(0,2) = (x)[1][2]-(x)[0][2]+((local)[1]+(local)[2])*c;
			JT(1,0) = (x)[3][0]-(x)[0][0]+(local)[0]*a;
			JT(1,1) = (x)[3][1]-(x)[0][1]+(local)[0]*b;
			JT(1,2) = (x)[3][2]-(x)[0][2]+(local)[0]*c;
			JT(2,0) = (x)[4][0]-(x)[0][0]+(local)[0]*a;
			JT(2,1) = (x)[4][1]-(x)[0][1]+(local)[0]*b;
			JT(2,2) = (x)[4][2]-(x)[0][2]+(local)[0]*c;}
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


template <>
class reference_element_traits<Pyramid>
{
	public:
		typedef ReferencePyramid reference_element_type;
};


}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_PYRAMID__ */
