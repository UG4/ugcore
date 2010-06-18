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
	static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_EDGE;
	static const int dim = 1;
	static const int num_corners = 2;
	static const int num_edges = 1;
	static const int num_faces = 0;
	static const int num_volumes = 0;
public:
	ReferenceEdge(){initializeArrays();}

	/* Dimension where reference element lives */
	int dimension() const {return dim;}

	/* size of reference triangle */
	number size() const	{return 1.0;}

	/* coordinates of reference corner Nr i (i=0..numberOfCorners) */
	const MathVector<dim>& corner(int i) const {return m_corner[i];}

	/* number of reference elements with id of dimension 'dim' of this reference element */
	unsigned int num_ref_elem(ReferenceElementType type) const {return m_ref_elem[type];}

	/* reference element type of subObject nr i of dimension dim_i */
	ReferenceElementType ref_elem_type(int dim_i, int i) const{	return m_ref_elem_type[dim_i][i];}

	unsigned int num_obj(int dim) const	{return m_num_obj[dim];}

	unsigned int num_obj_of_obj(int dim_i, int i, int dim_j) const
		{return m_num_obj_of_obj[dim_i][i][dim_j];}

	int id(int dim_i, int i, int dim_j, int j) const
		{return m_id[dim_i][i][dim_j][j];}

	~ReferenceEdge()
	{}

private:
	// to make it more readable
	enum{POINT = 0, EDGE = 1};
	enum{MAXOBJECTS = 2};

	/* number of Geometric Objects of Reference Element
	 * (m_num_obj[dim] = number of GeomObjects of dimension dim) */
	unsigned int m_num_obj[dim+1];
	/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
	unsigned int m_num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
	/* coordinates of Reference Corner */
	MathVector<dim> m_corner[2];
	// indices of GeomObjects
	int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

	unsigned int m_ref_elem[NUM_REFERENCE_ELEMENTS];
	ReferenceElementType m_ref_elem_type[dim+1][MAXOBJECTS];

	void initializeArrays()
	{
		//number of Geometric Objects
	 	m_num_obj[POINT] = 2;
	 	m_num_obj[EDGE] = 1;

		// number of Geometric Objects
	 	m_num_obj_of_obj[EDGE][0][POINT] = 2;
	 	m_num_obj_of_obj[EDGE][0][EDGE] = 1;

	 	m_ref_elem_type[EDGE][0] = RET_EDGE;

	 	for(unsigned int i = 0; i < m_num_obj[POINT]; ++i)
	 	{
		 	m_num_obj_of_obj[POINT][i][POINT] = 1;
		 	m_num_obj_of_obj[POINT][i][EDGE] = 1;

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

		// Points <-> Face
		for(unsigned int i=0; i<m_num_obj[POINT]; ++i)
		{
		 	m_id[EDGE][0][POINT][i] = i;
		 	m_id[POINT][i][EDGE][0] = 0;
		}

		// Reference Corners
	 	m_corner[0][0] = 0.0;
	 	m_corner[1][0] = 1.0;

	 	// Reference Element Types
	 	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	 	{
			m_ref_elem[i] = 0;
	 	}
	 	m_ref_elem[RET_POINT] = 2;
	 	m_ref_elem[RET_EDGE] = 1;
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
				JT[0][i] = a10[i];
			}
			return true;
		}

		bool jacobian_transposed_inverse(	const MathVector<dim>& loc_pos,
											MathMatrix<world_dim, dim>& JTInv) const
		{
			MathMatrix<dim, world_dim> JT;

			if(!jacobian_transposed(loc_pos, JT)) return false;

			if( (world_dim == 1) && (dim==1) )
			{
				JTInv[0][0] = 1./ JT[0][0];
				return true;
			}

			//TODO: Implement pseudo inverse
			return false;
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
