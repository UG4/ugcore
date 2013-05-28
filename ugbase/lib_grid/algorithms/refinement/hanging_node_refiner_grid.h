// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)

#ifndef __H__UG__HANGING_NODE_REFINER_GRID__
#define __H__UG__HANGING_NODE_REFINER_GRID__

#include "hanging_node_refiner_base.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	Specialization of ug::HangingNodeRefiner for ug::Grid
/**	This class should be used, if hanging node refinement shall be
 * applied on a flat grid (ug::Grid).
 *
 * Marked elements will be replaced by their newly created children.
 *
 * Take a look at ug::HangingNodeRefinerBase for a more in-depth documentation.
 *
 * \sa ug::HangingNodeRefinerBase, ug::HangingNodeRefiner_Grid
 */
class HangingNodeRefiner_Grid : public HangingNodeRefinerBase
{
	public:
		using HangingNodeRefinerBase::mark;

	public:
		HangingNodeRefiner_Grid(IRefinementCallback* refCallback = NULL);
		HangingNodeRefiner_Grid(Grid& grid,
								IRefinementCallback* refCallback = NULL);

		virtual ~HangingNodeRefiner_Grid();

		virtual void grid_to_be_destroyed(Grid* grid);

		virtual void assign_grid(Grid& grid);
		virtual Grid* get_associated_grid()		{return m_pGrid;}
		virtual Grid* grid()					{return m_pGrid;}

		virtual bool adaptivity_supported() const	{return true;}
		virtual bool coarsening_supported() const	{return false;}

	///	Marks a vertex for refinement (ignores RM_COARSEN).
		virtual bool mark(VertexBase* v, RefinementMark refMark = RM_REFINE);

	///	Marks an edge for refinement (ignores RM_COARSEN).
		virtual bool mark(EdgeBase* e, RefinementMark refMark = RM_REFINE);

	///	Marks a face for refinement (ignores RM_COARSEN).
		virtual bool mark(Face* f, RefinementMark refMark = RM_REFINE);

	///	Marks a volume for refinement (ignores RM_COARSEN).
		virtual bool mark(Volume* v, RefinementMark refMark = RM_REFINE);

	protected:
	///	performs registration and deregistration at a grid.
	/**	Initializes all grid related variables.
	 *  call set_grid(NULL) to unregister the observer from a grid.
	 *
	 * 	Please note that though the base grid features a set_grid method,
	 *  it is not declared virtual. This is because we want to call it
	 *  during construction and destruction.*/
		void set_grid(Grid* grid);

	///	erases unused refined elements
		virtual void post_refine();

		virtual void process_constraining_edge(ConstrainingEdge* cge);
		virtual void refine_edge_with_normal_vertex(EdgeBase* e,
											VertexBase** newCornerVrts = NULL);

		virtual void refine_face_with_normal_vertex(Face* f,
											VertexBase** newCornerVrts = NULL);
		virtual void process_constraining_face(ConstrainingFace* cgf);

		virtual void refine_volume_with_normal_vertex(Volume* v,
											VertexBase** newVolumeVrts = NULL);

	///	Returns the vertex associated with the edge
		virtual VertexBase* get_center_vertex(EdgeBase* e);

	///	Associates a vertex with the edge.
		virtual void set_center_vertex(EdgeBase* e, VertexBase* v);

	///	Returns the vertex associated with the face
		virtual VertexBase* get_center_vertex(Face* f);

	///	Associates a vertex with the face.
		virtual void set_center_vertex(Face* f, VertexBase* v);

	private:
		Grid* 			m_pGrid;
		AVertexBase		m_aVertex;
		Grid::EdgeAttachmentAccessor<AVertexBase>		m_aaVertexEDGE;
		Grid::FaceAttachmentAccessor<AVertexBase>		m_aaVertexFACE;
};

/// @}	// end of add_to_group command

}//	end of namespace

#endif
