// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y09 m03 d16

#ifndef __H__LIB_GRID__HANGING_NODE_REFINER_2D_IRN__
#define __H__LIB_GRID__HANGING_NODE_REFINER_2D_IRN__

#include <queue>
#include "lib_grid/lg_base.h"
#include "refinement_callbacks.h"
#include "refiner_interface.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
//	HangingNodeRefiner2D_IRN
///	Performs adaptive grid refinement with hanging nodes (DEPRECIATED AND UNSUPPORTED).
/** DEPRECIATED AND UNSUPPORTED
 *
 * IRN stands for "irregularity rule n", which intends to say that
 * an arbitrary number of hanging nodes may exist on each edge.
 *
 * Refines flat 2d grids (no multigrids) using hanging node refinement
 * with an arbitrary number of hanging nodes per edge.
 *
 * Most commonly you won't use this refiner, but the ug::HangingNodeRefiner_Grid
 * or ug::HanginNodeRefiner_MultiGrid, which only support one hanging node per
 * edge, but which are more robust and work in 3d as well.
 *
 * Initialize the Refiner with a grid or manually assign one later on.
 * by marking elements of the grid using mark_for_refinement,
 * you can define where the grid shall be refined.
 * A call to refine() will finally create the new geometry.
 *
 * Make sure, that you won't delete any marked elements from the grid.
 * Please note, that a grid, which was refined by this class,
 * will most commonly contain HangingNodes, Constrained- and
 * Constraining- Edges/Faces. Such a grid has to be treated with care,
 * since those elements are tightly connected and reference each other.
 *
 * You may specify a refinement callback which can be used to calculate
 * the new positions of newly generated vertices. You can specify it in
 * the constructor or through an explicit call to set_refinement_callback.
 * If no refinement-callback is specified, HangingNodeRefiner will
 * attempt to create a standard linear refinement callback for an attached
 * standard position attachment (ug::aPosition or ug::aPosition2).
 */

class HangingNodeRefiner2D_IRN : public IRefiner, public GridObserver
{
	public:
		using IRefiner::mark_for_refinement;

	public:
		HangingNodeRefiner2D_IRN(IRefinementCallback* refCallback = NULL);
		HangingNodeRefiner2D_IRN(Grid& grid, IRefinementCallback* refCallback = NULL);
	//todo: make copy-constructor public
		virtual ~HangingNodeRefiner2D_IRN();

		virtual void grid_to_be_destroyed(Grid* grid);

		void assign_grid(Grid& grid);
		virtual Grid* get_associated_grid()		{return m_pGrid;}
		
		virtual void clear_marks();
		virtual void mark_for_refinement(Edge* e);
		virtual void mark_for_refinement(Face* f);
		virtual void mark_for_refinement(Volume* v);

		bool set_irregularity_rule(uint irregularityRule);
		uint get_irregularity_rule();

	///	performs refinement on the marked elements.
	/**
	 * irregularityRule specified the number of hanging nodes that may coexist
	 * on one edge. If there are more, then the edge will be refined.
	 */
		virtual void refine();

	protected:
		typedef Attachment<int> ARefinementMark;
		//typedef Attachment<RefinementInfo> ARefinementInfo;

		typedef std::vector<Edge*>	EdgeVec;
		typedef std::vector<Face*>		FaceVec;
		typedef std::vector<Volume*>	VolumeVec;

		typedef std::queue<Edge*> EdgeQueue;
		typedef std::queue<Face*> FaceQueue;
		typedef std::queue<Volume*> VolumeQueue;

	protected:
	///	performs registration and deregistration at a grid.
	/**	call set_grid(NULL) to unregister the observer from a grid.*/
		void set_grid(Grid* grid);
		
	///	marks unmarked elements that have to be refined due to marked neighbors.
	/**
	 * all elements that have to be refined will be written to the passed queues.
	 * Note that this will most likely be more elements than just the marked ones.
	 *
	 * This method is virtual to allow derivates to mark additional elements as required.
	 * Normally a a derived class will first call the method of its this class and
	 * the perform its own operations.
	 */
		virtual void collect_objects_for_refine();

	///	this method helps derived classes to perform operations directly before actual element refinment is performed.
	/**	Called from the refine() method in each refinement-iteration after
	 *	collect_objects_for_refine().
	 *	Default implementation is empty.*/
		virtual void refinement_step_begins()	{};

	///	this method helps derived classes to perform operations directly after actual element refinment took place.
	/**	Called from the refine() method in each refinement-iteration after
	 *	all scheduled elements had been refined.
	 *	The refine process will either terminate after this method or will
	 *	start a new iteration, if new elements had been marked during refine.
	 *	Default implementation is empty.*/
		virtual void refinement_step_ends()		{};

	///	two constrained edges will be created. The old edge remains and has to be deleted later on.
		void refine_constrained_edge(ConstrainedEdge* constrainedEdge);

	///	two edges will be created. The old edge remains and has to be deleted later on.
	/**
	 * Both returned edges can possibly be constraining-edges that have to be refined
	 * once more (due to the irregularity-rule). Be sure to check this.
	 */
		void refine_constraining_edge(ConstrainingEdge* constrainingEdge,
									Edge** ppEdge1Out, Edge** ppEdge2Out,
									bool& scheduledEdgeReplaced1Out, bool& scheduledEdgeReplaced2Out);

	///	two normal edges will be created. The old edge remains and has to be deleted later on.
		void refine_edge_with_normal_vertex(Edge* e);

	///	e will be replaced by a constraining edge. Two new constrained edges will be created.
		void refine_edge_with_hanging_vertex(Edge* e);

	///	f will be refined by Face::refine. If f has more than 3 corners, a new vertex will be inserted on f.
		void refine_face_with_normal_vertex(Face* f);

	///	f will be replaced by a constraining face. new constrained faces will be created.
		void refine_face_with_hanging_vertex(Face* f);

	///	new faces will be created.
	/**
	 * returned faces can possibly be constraining-faces that have to be refined
	 * once more (due to the irregularity-rule). Be sure to check this.
	 */
		void refine_constraining_face(std::vector<Face*>& vNewFacesOut,
									ConstrainingFace* constrainingFace);

	///	v will be refined by Volume::refine. Eventually a new vertex will be inserted on v.
		void refine_volume_with_normal_vertex(Volume* v);

	////////////////////////////////////////////////////////////////////////
	//	helpers. Make sure that everything is initialized properly
	//	before calling these methods.
	//	you should use this methods instead of directly marking elements.
		inline bool is_marked(Edge* e)							{return m_selMarkedElements.is_selected(e);}
		inline void mark(Edge* e)								{m_selMarkedElements.select(e);}
		inline Vertex* get_center_vertex(Edge* e)			{return m_aaVertexEDGE[e];}
		inline void set_center_vertex(Edge* e, Vertex* v)	{m_aaVertexEDGE[e] = v;}

		inline bool is_marked(Face* f)								{return m_selMarkedElements.is_selected(f);}
		inline void mark(Face* f)									{m_selMarkedElements.select(f);}
		inline Vertex* get_center_vertex(Face* f)				{return m_aaVertexFACE[f];}
		inline void set_center_vertex(Face* f, Vertex* v)		{m_aaVertexFACE[f] = v;}

		inline bool is_marked(Volume* v)							{return m_selMarkedElements.is_selected(v);}
		inline void mark(Volume* v)									{m_selMarkedElements.select(v);}

	private:
		HangingNodeRefiner2D_IRN(const HangingNodeRefiner2D_IRN& hnr);
		
	protected:
	///	RefinementMarks are used throughout the refinement-process to indicate how an element shall be processed.
	/*
		enum RefinementMark
		{
			RM_NONE = 0,
			RM_REFINE,
			RM_LOCKED
		};

	///	The RefinementInfo will be attached to edges and faces and is used to store temporary refinement information.
		class RefinementInfo
		{
			public:
				RefinementInfo()	{reset();}

				inline void reset()			{m_mark = RM_NONE; m_vertex = NULL;}

				inline void set_vertex(Vertex* v)		{m_vertex = v;}
				inline void set_mark(RefinementMark mark)	{m_mark = mark;}

				inline Vertex* get_vertex()	{return m_vertex;}
				inline int get_mark()			{return m_mark;}

			protected:
				Vertex* m_vertex;
				int m_mark;
		};
	*/


	protected:
		Grid*		m_pGrid;

		AVertex		m_aVertex;

		uint 		m_irregularityRule;
		Selector	m_selMarkedElements;
		Selector	m_selScheduledElements;

		//ARefinementInfo m_aRefinementInfo;

		//Grid::VertexAttachmentAccessor<APosition>		m_aaPos;
		//Grid::EdgeAttachmentAccessor<ARefinementInfo>	m_aaRefinementInfoEDGE;
		//Grid::FaceAttachmentAccessor<ARefinementInfo>	m_aaRefinementInfoFACE;
		Grid::EdgeAttachmentAccessor<AVertex>		m_aaVertexEDGE;
		Grid::FaceAttachmentAccessor<AVertex>		m_aaVertexFACE;
		//Grid::VolumeAttachmentAccessor<ARefinementMark>	m_aaRefinementMarkVOL;
};


/// @}

}// end of namespace

#endif
