/*
 * reference_tetrahedron.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef m__H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__
#define m__H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__

namespace ug{


class ReferenceTetrahedron{
	public:
		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_TETRAHEDRON;
		static const int dim = 3;
		static const int num_corners = 4;
		static const int num_edges = 6;
		static const int num_faces = 4;
		static const int num_volumes = 1;

		ReferenceTetrahedron(){initializeArrays();}

		/* Dimension where reference element lives */
		int dimension() const{return dim;}

		/* size of reference triangle */
		number size() const	{return 1./6.;}

		/* coordinates of reference corner Nr i (i=0..numberOfCorners) */
		const MathVector<dim>& corner(size_t i) const {return m_corner[i];}

		/* number of reference elements with id of dimension 'dim' of this reference element */
		size_t num_ref_elem(ReferenceElementType type) const {return m_ref_elem[type];}

		/* reference element type of subObject nr i of dimension dim_i */
		ReferenceElementType ref_elem_type(int dim_i, int i) const{	return m_ref_elem_type[dim_i][i];}

		unsigned int num_obj(int dim) const	{return m_num_obj[dim];}

		unsigned int num_obj_of_obj(int dim_i, int i, int dim_j) const
			{return m_num_obj_of_obj[dim_i][i][dim_j];}

		int id(int dim_i, int i, int dim_j, int j) const
			{return m_id[dim_i][i][dim_j][j];}

		~ReferenceTetrahedron()
		{}

	private:
		// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME= 3};
		enum{MAXOBJECTS = 6};

		/* number of Geometric Objects of Reference Element
		 * (m_num_obj[dim] = number of GeomObjects of dimension dim) */
		unsigned int m_num_obj[dim+1];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int m_num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
		/* coordinates of Reference Corner */
		MathVector<dim> m_corner[num_corners];
		// indices of GeomObjects
		int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

		unsigned int m_ref_elem[NUM_REFERENCE_ELEMENTS];
		ReferenceElementType m_ref_elem_type[dim+1][MAXOBJECTS];

		void initializeArrays()
		{
			//number of Geometric Objects
		 	m_num_obj[POINT] = 4;
		 	m_num_obj[EDGE] = 6;
		 	m_num_obj[FACE] = 4;
		 	m_num_obj[VOLUME] = 1;

			// number of Geometric Objects
		 	m_num_obj_of_obj[VOLUME][0][POINT] = 4;
		 	m_num_obj_of_obj[VOLUME][0][EDGE] = 6;
		 	m_num_obj_of_obj[VOLUME][0][FACE] = 4;
		 	m_num_obj_of_obj[VOLUME][0][VOLUME] = 1;

		 	for(unsigned int i = 0; i < m_num_obj[FACE]; ++i)
		 	{
				m_num_obj_of_obj[FACE][i][POINT] = 3;
				m_num_obj_of_obj[FACE][i][EDGE] = 3;
				m_num_obj_of_obj[FACE][i][FACE] = 1;
				m_num_obj_of_obj[FACE][i][VOLUME] = 1;

				m_ref_elem_type[FACE][i] = RET_TRIANGLE;
		 	}

		 	for(unsigned int i = 0; i < m_num_obj[EDGE]; ++i)
		 	{
			 	m_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	m_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	m_num_obj_of_obj[EDGE][i][FACE] = 2;
			 	m_num_obj_of_obj[EDGE][i][VOLUME] = 1;

			 	m_ref_elem_type[EDGE][i] = RET_EDGE;
		 	}

		 	for(unsigned int i = 0; i < m_num_obj[POINT]; ++i)
		 	{
			 	m_num_obj_of_obj[POINT][i][POINT] = 1;
			 	m_num_obj_of_obj[POINT][i][EDGE] = 3;
			 	m_num_obj_of_obj[POINT][i][FACE] = 3;
			 	m_num_obj_of_obj[POINT][i][VOLUME] = 1;

			 	m_ref_elem_type[POINT][i] = RET_POINT;
		 	}

			//reset m_id to -1
			for(int i=0; i<=dim; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(int k=0; k<=dim; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
						 	m_id[i][j][k][l] = -1;
						}

			//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
			for(int i=0; i<=dim; ++i)
				for(unsigned int j=0; j<m_num_obj[i]; ++j)
				{
				 	m_id[i][j][i][0] = j;
				}

			// Face <-> Volume
			for(unsigned int i=0; i<m_num_obj[VOLUME]; ++i)
			{
			 	m_id[VOLUME][0][FACE][i] = i;
			 	m_id[FACE][i][VOLUME][0] = 0;
			}

			// Edge <-> Volume
			for(unsigned int i=0; i<m_num_obj[EDGE]; ++i)
			{
			 	m_id[VOLUME][0][EDGE][i] = i;
			 	m_id[EDGE][i][VOLUME][0] = 0;
			}

			// Point <-> Volume
			for(unsigned int i=0; i<m_num_obj[POINT]; ++i)
			{
			 	m_id[VOLUME][0][POINT][i] = i;
			 	m_id[POINT][i][VOLUME][0] = 0;
			}

			// Points <-> Faces
		 	m_id[FACE][0][POINT][0] = 0;
		 	m_id[FACE][0][POINT][1] = 2;
		 	m_id[FACE][0][POINT][2] = 1;

		 	m_id[FACE][1][POINT][0] = 1;
		 	m_id[FACE][1][POINT][1] = 2;
		 	m_id[FACE][1][POINT][2] = 3;

		 	m_id[FACE][2][POINT][0] = 0;
		 	m_id[FACE][2][POINT][1] = 3;
		 	m_id[FACE][2][POINT][2] = 2;

		 	m_id[FACE][3][POINT][0] = 0;
		 	m_id[FACE][3][POINT][1] = 1;
		 	m_id[FACE][3][POINT][2] = 3;


		 	m_id[POINT][0][FACE][0] = 0;
		 	m_id[POINT][0][FACE][1] = 2;
		 	m_id[POINT][0][FACE][2] = 3;

		 	m_id[POINT][1][FACE][0] = 0;
		 	m_id[POINT][1][FACE][1] = 1;
		 	m_id[POINT][1][FACE][2] = 3;

		 	m_id[POINT][2][FACE][0] = 0;
		 	m_id[POINT][2][FACE][1] = 1;
		 	m_id[POINT][2][FACE][2] = 2;

		 	m_id[POINT][3][FACE][0] = 1;
		 	m_id[POINT][3][FACE][1] = 2;
		 	m_id[POINT][3][FACE][2] = 3;

		 	// Edges <-> Faces
		 	m_id[FACE][0][EDGE][0] = 0;
		 	m_id[FACE][0][EDGE][1] = 2;
		 	m_id[FACE][0][EDGE][2] = 1;

		 	m_id[FACE][1][EDGE][0] = 1;
		 	m_id[FACE][1][EDGE][1] = 5;
		 	m_id[FACE][1][EDGE][2] = 4;

		 	m_id[FACE][2][EDGE][0] = 2;
		 	m_id[FACE][2][EDGE][1] = 3;
		 	m_id[FACE][2][EDGE][2] = 5;

		 	m_id[FACE][3][EDGE][0] = 0;
		 	m_id[FACE][3][EDGE][1] = 4;
		 	m_id[FACE][3][EDGE][2] = 3;

		 	m_id[EDGE][0][FACE][0] = 0;
		 	m_id[EDGE][0][FACE][1] = 3;

		 	m_id[EDGE][1][FACE][0] = 0;
		 	m_id[EDGE][1][FACE][1] = 1;

		 	m_id[EDGE][2][FACE][0] = 0;
		 	m_id[EDGE][2][FACE][1] = 2;

		 	m_id[EDGE][3][FACE][0] = 3;
		 	m_id[EDGE][3][FACE][1] = 2;

		 	m_id[EDGE][4][FACE][0] = 1;
		 	m_id[EDGE][4][FACE][1] = 3;

		 	m_id[EDGE][5][FACE][0] = 2;
		 	m_id[EDGE][5][FACE][1] = 1;


			// Points of Edges
			// edge 0 = (0,1)
		 	m_id[EDGE][0][POINT][0] = 0;
		 	m_id[EDGE][0][POINT][1] = 1;
		 	// edge 1 = (1,2)
		 	m_id[EDGE][1][POINT][0] = 1;
		 	m_id[EDGE][1][POINT][1] = 2;
			// edge 2 = (2,0)
		 	m_id[EDGE][2][POINT][0] = 2;
		 	m_id[EDGE][2][POINT][1] = 0;
			// edge 3 = (0,3)
		 	m_id[EDGE][3][POINT][0] = 0;
		 	m_id[EDGE][3][POINT][1] = 3;
			// edge 4 = (1,3)
		 	m_id[EDGE][4][POINT][0] = 1;
		 	m_id[EDGE][4][POINT][1] = 3;
			// edge 5 = (2,3)
		 	m_id[EDGE][5][POINT][0] = 2;
		 	m_id[EDGE][5][POINT][1] = 3;

		 	// Edges of Point
		 	m_id[POINT][0][EDGE][0] = 2;
		 	m_id[POINT][0][EDGE][1] = 0;
		 	m_id[POINT][0][EDGE][2] = 3;

		 	m_id[POINT][1][EDGE][0] = 0;
		 	m_id[POINT][1][EDGE][1] = 1;
		 	m_id[POINT][1][EDGE][2] = 4;

		 	m_id[POINT][2][EDGE][0] = 1;
		 	m_id[POINT][2][EDGE][1] = 2;
		 	m_id[POINT][2][EDGE][2] = 5;

		 	m_id[POINT][3][EDGE][0] = 3;
		 	m_id[POINT][3][EDGE][1] = 4;
		 	m_id[POINT][3][EDGE][2] = 5;

			// Reference Corners
		 	m_corner[0] = MathVector<dim>(0.0, 0.0, 0.0);
		 	m_corner[1] = MathVector<dim>(1.0, 0.0, 0.0);
		 	m_corner[2] = MathVector<dim>(0.0, 1.0, 0.0);
		 	m_corner[4] = MathVector<dim>(0.0, 0.0, 1.0);

		 	// Reference Element Types
		 	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
		 	{
				m_ref_elem[i] = 0;
		 	}
		 	m_ref_elem[RET_POINT] = 4;
		 	m_ref_elem[RET_EDGE] = 6;
		 	m_ref_elem[RET_TRIANGLE] = 4;
		 	m_ref_elem[RET_TETRAHEDRON] = 1;
		}


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

			if( (world_dim == 3) && (dim==3) )
			{
				typename MathMatrix<3,3>::value_type det;
				det = JT(0,0)*JT(1,1)*JT(2,2)
				+ JT(0,1)*JT(1,2)*JT(2,0)
				+ JT(0,2)*JT(1,0)*JT(2,1)
				- JT(0,0)*JT(1,2)*JT(2,1)
				- JT(0,1)*JT(1,0)*JT(2,2)
				- JT(0,2)*JT(1,1)*JT(2,0);

				assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
				typename MathMatrix<3,3>::value_type invdet = 1./det;

				JTInv(0,0) = ( JT(1,1)*JT(2,2) - JT(1,2)*JT(2,1)) * invdet;
				JTInv(0,1) = (-JT(0,1)*JT(2,2) + JT(0,2)*JT(2,1)) * invdet;
				JTInv(0,2) = ( JT(0,1)*JT(1,2) - JT(0,2)*JT(1,1)) * invdet;
				JTInv(1,0) = (-JT(1,0)*JT(2,2) + JT(1,2)*JT(2,0)) * invdet;
				JTInv(1,1) = ( JT(0,0)*JT(2,2) - JT(0,2)*JT(2,0)) * invdet;
				JTInv(1,2) = (-JT(0,0)*JT(1,2) + JT(0,2)*JT(1,0)) * invdet;
				JTInv(2,0) = ( JT(1,0)*JT(2,1) - JT(1,1)*JT(2,0)) * invdet;
				JTInv(2,1) = (-JT(0,0)*JT(2,1) + JT(0,1)*JT(2,0)) * invdet;
				JTInv(2,2) = ( JT(0,0)*JT(1,1) - JT(0,1)*JT(1,0)) * invdet;
				return true;
			}

			//TODO: Implement pseudo inverse
			return false;
		}

		bool jacobian_det(const MathVector<dim>& loc_pos, number& det) const
		{
			MathMatrix<dim, world_dim> JT;
			if(!jacobian_transposed(loc_pos, JT)) return false;
			det = JT(0,0)*JT(1,1)*JT(2,2)
			+ JT(0,1)*JT(1,2)*JT(2,0)
			+ JT(0,2)*JT(1,0)*JT(2,1)
			- JT(0,0)*JT(1,2)*JT(2,1)
			- JT(0,1)*JT(1,0)*JT(2,2)
			- JT(0,2)*JT(1,1)*JT(2,0);
			return true;
		}

	private:
		const MathVector<world_dim>* m_corners;

		MathVector<world_dim> a10, a20, a30;
};


template <>
class reference_element_traits<Tetrahedron>
{
	public:
		typedef ReferenceTetrahedron reference_element_type;
};


}

#endif /* m__H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__ */
