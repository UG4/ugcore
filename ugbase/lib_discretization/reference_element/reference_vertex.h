/*
 * reference_vertex.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_VERTEX__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_VERTEX__

namespace ug{

class ReferenceVertex
{
public:
	static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
	static const int dim = 0;
	static const int num_corners = 1;
	static const int num_edges = 0;
	static const int num_faces = 0;
	static const int num_volumes = 0;

public:
	ReferenceVertex() {initializeArrays();}

	/// reference object id
	ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// Dimension where reference element lives
	int dimension() const {return dim;}

	/// size of reference triangle
	number size() const	{return 1.0;}

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

private:
	// to make it more readable
	enum{POINT = 0, EDGE = 1};
	enum{MAXOBJECTS = 1};

	/* number of Geometric Objects of Reference Element
	 * (m_num_obj[dim] = number of GeomObjects of dimension dim) */
	size_t m_num_obj[dim+1];
	/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
	size_t m_num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
	// indices of GeomObjects
	int m_id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

	size_t m_ref_elem[NUM_REFERENCE_OBJECTS];
	ReferenceObjectID m_ref_elem_type[dim+1][MAXOBJECTS];

	void initializeArrays()
	{
		//number of Geometric Objects
		m_num_obj[POINT] = 1;

		// number of Geometric Objects
		m_num_obj_of_obj[POINT][0][POINT] = 1;

		m_ref_elem_type[POINT][0] = ROID_VERTEX;

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

		// Reference Element Types
		for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		{
			m_ref_elem[i] = 0;
		}
		m_ref_elem[ROID_VERTEX] = 1;
	}
};

template <>
class reference_element_traits<VertexBase>
{
	public:
		typedef ReferenceVertex reference_element_type;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_VERTEX__ */
