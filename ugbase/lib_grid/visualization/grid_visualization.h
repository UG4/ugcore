// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.10.2011 (m,d,y)

#ifndef __H__UG__grid_visualization__
#define __H__UG__grid_visualization__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

///	This class generates arrays of primitives, which can be used for the visualization of a grid
/**	It is suggested that you keep one GridVisualization object for each grid.
 * If you change contents of the associated grid (positions, new elements, ...),
 * you should call the method update_visuals().
 */
template <class TNumber, class TInt, class TAVrtPos>
class GridVisualization
{
	public:
	///	Element constants are used to enable or disable rendering of elements.
		enum Elements{
			E_VERTEX = 1,
			E_EDGE = 1<<1,
			E_FACE = 1<<2,
			E_VOLUME = 1<<3
		};

	///	The visual type tells, what a visual represents.
		enum VisualType{
			E_GRID_SUBSET,
			E_SELECTION,
			E_MARKS
		};

	public:
		GridVisualization();
		~GridVisualization();

		void set_grid(Grid& grid, TAVrtPos& aVrtPos);
		//void set_subset_handler(Grid& sh);

	///	Set visibility for a subset
		//void hide_subset(int subsetIndex);

	///	returns whether a subset is hidden
		//bool subset_is_hidden(int subsetIndex);

	///	defines which elements will be considered for geometry updates.
	/**	Pass or-combinations of constants enumerated in GridVisualization::Elements.*/
		//void set_considered_elements(uint elems);

	///	returns which elements are considered for geometry updates.
	/**	Returns or-combinations of constants enumerated in GridVisualization::Elements.*/
		//uint considered_elements();

	///	updates and creates visuals.
	/**	Call this method if the underlying grid has changed. All visuals will
	 * be recomputed. Note that all arrays returned by methods of this class
	 * are invalidated by a call to update_visuals. You thus should make sure
	 * to re-retrieve those arrays by calls to the appropriate methods.*/
		void update_visuals();


	///	returns the number of vertices
		int num_vertices();

	///	returns the array of vertex positions.
	/**	The returned array contains num_vertices() * 3 entries. Three consecutive
	 * entries represent the x, y and z coordinate of one vertex (xyzxyzxyz...)*/
		const TNumber* vertex_positions();


	///	returns the number of visuals
		int num_visuals();

	///	returns the visual type. One of the constants enumerated in GridVisualization::VisualType
		int visual_type(int visInd);

	///	returns an array of 4 numbers, describing the color of the specified visual (rgba).
		const TNumber* visual_color(int visInd);

	///	returns an array of face normals for the specified visual
	/**	The returned array contains num_faces() * 3 entries. Three consecutive
	 * entries represent the x, y and z component of the normal of a face.
	 * The first part of the array (num_triangles(visInd) * 3 entries) contain
	 * the normals for the triangles, the second part
	 * (num_quadrilaterals(visInd) * 3 entries) contains the normals for the
	 * quadrilaterals.*/
		const TNumber* face_normals(int visInd);

	///	returns the number of points that shall be rendered
		//int num_points(int visInd);

	///	returns a list describing which vertices shall be rendered.
	/**	Returns a list of vertex indices.
	 * The returned array has the size num_points(visInd).*/
		//const TInt* point_list(int visInd);

	///	returns the number of edges that shall be rendered
		//int num_edges(int visInd);

	///	returns a list describing which edges shall be rendered.
	/**	Returns a list of vertex indices.
	 * Two consecutive indices represent an edge. The returned array has the size
	 * num_edges(visInd) * 2.*/
		//const TInt* edge_list(int visInd);

	///	returns the number of faces that shall be rendered.
	/**	num_faces(visInd) == num_triangles(visInd) + num_quadrilaterals(visInd);*/
		//int num_faces(int visInd);

	///	returns the number of triangles that shall be rendered
		int num_triangles(int visInd);

	///	returns a list describing which triangles shall be rendered.
	/**	Returns a list of vertex indices.
	 * Three consecutive indices represent a triangle in counter-clockwise order.
	 * The returned array has the size num_triangles(visInd) * 3.*/
		const TInt* triangle_list(int visInd);

	///	returns the number of quadrilaterals that shall be rendered
		//int num_quadrilaterals(int visInd);

	///	returns a list describing which quadrilaterals shall be rendered.
	/**	Returns a list of vertex indices.
	 * Four consecutive indices represent a quadrilateral in counter-clockwise order.
	 * The returned array has the size num_quadrilaterals(visInd) * 4.*/
		//const TInt* quadrilateral_list(int visInd);
	
	protected:
		void release_grid();
		
	protected:		
		struct Visual{
			std::vector<TNumber>	m_faceNormals;
			std::vector<TInt>		m_triIndices;
			TNumber					m_color[4];
			int 					m_visType;
		};
		typedef SmartPtr<Visual>	SPVisual;
		
		Grid*	m_pGrid;
		std::vector<TNumber>	m_vrtPositions;
		std::vector<SPVisual>	m_visuals;


		AInt		m_aVrtInd;
		TAVrtPos	m_aVrtPos;
		AVector3	m_aNormal;

		Grid::VertexAttachmentAccessor<AInt>		m_aaIndVRT; ///< indices the first coordinate in m_vrtPositions
		Grid::VertexAttachmentAccessor<TAVrtPos>	m_aaPosVRT;
		Grid::FaceAttachmentAccessor<AVector3>		m_aaNormFACE;
};

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "grid_visualization_impl.hpp"

#endif
