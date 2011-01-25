// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)

#ifndef __H__UG__HANGIN_NODE_REFINER_MULTI_GRID__
#define __H__UG__HANGIN_NODE_REFINER_MULTI_GRID__

#include "hanging_node_refiner_base.h"

namespace ug
{

class HangingNodeRefiner_MultiGrid : public HangingNodeRefinerBase
{
	public:
		HangingNodeRefiner_MultiGrid(IRefinementCallback* refCallback = NULL);
		HangingNodeRefiner_MultiGrid(MultiGrid& mg,
									IRefinementCallback* refCallback = NULL);

		virtual ~HangingNodeRefiner_MultiGrid();

		void assign_grid(MultiGrid& mg);
		virtual Grid* get_associated_grid()		{return m_pMG;}

	protected:
	///	performs registration and deregistration at a grid.
	/**	Initializes all grid related variables.
	 *  call set_grid(NULL) to unregister the observer from a grid.
	 *
	 * 	Please note that though the base grid features a set_grid method,
	 *  it is not declared virtual. This is because we want to call it
	 *  during construction and destruction.*/
		void set_grid(MultiGrid* mg);

	///	prepares selection and calls the base implementation
	/**	Makes sure that no elements with children are selected,
	 *  Additionally vertices are marked, which have no children but
	 *  associated elements which are marked for refinement.*/
		virtual void collect_objects_for_refine();

	///	creates required vertices in higher levels.
		virtual void pre_refine();

	/**	Calls the base implementation and passes the mg-child vertices as
	 *  newCornerVrts. If newCornerVrts were passed to this method, they
	 *  are ignored.
	 *  \{*/
		virtual void refine_edge_with_normal_vertex(EdgeBase* e,
											VertexBase** newCornerVrts = NULL);
		virtual void refine_edge_with_hanging_vertex(EdgeBase* e,
											VertexBase** newCornerVrts = NULL);

		virtual void refine_face_with_normal_vertex(Face* f,
											VertexBase** newCornerVrts = NULL);
		virtual void refine_face_with_hanging_vertex(Face* f,
											VertexBase** newCornerVrts = NULL);

		virtual void refine_volume_with_normal_vertex(Volume* v,
											VertexBase** newVolumeVrts = NULL);
	/*	\} */

	///	Returns the vertex associated with the edge
		virtual VertexBase* get_center_vertex(EdgeBase* e);

	///	Associates a vertex with the edge.
		virtual void set_center_vertex(EdgeBase* e, VertexBase* v);

	///	Returns the vertex associated with the face
		virtual VertexBase* get_center_vertex(Face* f);

	///	Associates a vertex with the face.
		virtual void set_center_vertex(Face* f, VertexBase* v);

	///	calls base implementation and replaces cge with a normal edge.
		virtual void refine_constraining_edge(ConstrainingEdge* cge);


	///	collects corner vertices and fills them into the associated vector
	/**	The size of cornersOut is automatically adjusted.
	 *  The i-th element of corners out will contain the child vertex of the
	 *  i-th vertex of elem.
	 */
		template <class TElem>
		void collect_child_corners(std::vector<VertexBase*>& cornersOut, TElem* elem)
		{
			cornersOut.resize(elem->num_vertices());
			for(size_t i = 0; i < elem->num_vertices(); ++i){
				//UG_ASSERT(m_pMG->get_child_vertex(elem->vertex(i)), "A child vertex has to exists!");
				cornersOut[i] = m_pMG->get_child_vertex(elem->vertex(i));
			}
		}

	private:
		MultiGrid*	m_pMG;
};

}//	end of namespace

#endif
