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
		static const int dim = 3;

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
		enum{POINT = 0, EDGE = 1, FACE = 2};
		enum{MAXOBJECTS = 3};

		/* number of Geometric Objects of Reference Element
		 * (m_num_obj[dim] = number of GeomObjects of dimension dim) */
		unsigned int m_num_obj[dim+1];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int m_num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
		/* coordinates of Reference Corner */
		MathVector<dim> m_corner[3];
		// indices of GeomObjects
		int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

		unsigned int m_ref_elem[NUM_REFERENCE_ELEMENTS];
		ReferenceElementType m_ref_elem_type[dim+1][MAXOBJECTS];

		void initializeArrays()
		{
			//number of Geometric Objects
		 	m_num_obj[POINT] = 3;
		 	m_num_obj[EDGE] = 3;
		 	m_num_obj[FACE] = 1;

			// number of Geometric Objects
		 	m_num_obj_of_obj[FACE][0][POINT] = 3;
		 	m_num_obj_of_obj[FACE][0][EDGE] = 3;
		 	m_num_obj_of_obj[FACE][0][FACE] = 1;

		 	m_ref_elem_type[FACE][0] = RET_TRIANGLE;

		 	for(unsigned int i = 0; i < m_num_obj[EDGE]; ++i)
		 	{
			 	m_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	m_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	m_num_obj_of_obj[EDGE][i][FACE] = 1;

			 	m_ref_elem_type[EDGE][i] = RET_EDGE;
		 	}

		 	for(unsigned int i = 0; i < m_num_obj[POINT]; ++i)
		 	{
			 	m_num_obj_of_obj[POINT][i][POINT] = 1;
			 	m_num_obj_of_obj[POINT][i][EDGE] = 2;
			 	m_num_obj_of_obj[POINT][i][FACE] = 1;

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

			//Edges <-> Face
			for(unsigned int i=0; i<m_num_obj[EDGE]; ++i)
			{
			 	m_id[FACE][0][EDGE][i] = i;
			 	m_id[EDGE][i][FACE][0] = 0;
			}

			// Points <-> Face
			for(unsigned int i=0; i<m_num_obj[POINT]; ++i)
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
			// edge 2 = (2,0)
		 	m_id[EDGE][2][POINT][0] = 2;
		 	m_id[EDGE][2][POINT][1] = 0;

		 	// Edges of Point
		 	m_id[POINT][0][EDGE][0] = 2;
		 	m_id[POINT][0][EDGE][1] = 0;

		 	m_id[POINT][1][EDGE][0] = 0;
		 	m_id[POINT][1][EDGE][1] = 1;

		 	m_id[POINT][2][EDGE][0] = 1;
		 	m_id[POINT][2][EDGE][1] = 2;


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
		 	m_ref_elem[RET_POINT] = 3;
		 	m_ref_elem[RET_EDGE] = 3;
		 	m_ref_elem[RET_TRIANGLE] = 1;
		}


};

template <>
class reference_element_traits<ReferenceTetrahedron>
{
	public:
		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_TRIANGLE;

		static const int dim = 2;
		static const int num_corners = 3;
		static const int num_edges = 3;
		static const int num_faces = 1;
		static const int num_volumes = 0;

		typedef MathVector<dim> position_type;
};

template <>
class reference_element_traits<Tetrahedron>
{
	public:
		typedef ReferenceTriangle reference_element_type;
		typedef reference_element_traits<ReferenceTriangle> reference_element_type_traits;

		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = reference_element_type_traits::REFERENCE_ELEMENT_TYPE;
		static const int dim = reference_element_type_traits::dim;
		static const int num_corners = reference_element_type_traits::num_corners;
		static const int num_edges = reference_element_type_traits::num_edges;
		static const int num_faces = reference_element_type_traits::num_faces;
		static const int num_volumes = reference_element_type_traits::num_volumes;
		typedef reference_element_type_traits::position_type position_type;
};


}

#endif /* m__H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_TETRAHEDRON__ */
