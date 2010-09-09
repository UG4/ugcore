// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y09 m03 d16

#ifndef __H__LIB_GRID__HANGING_NODE_REFINER_IR1__
#define __H__LIB_GRID__HANGING_NODE_REFINER_IR1__

#include <queue>
#include "lib_grid/lg_base.h"
#include "refinement_callbacks.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
//	HangingNodeRefiner_IR1
///	Performs adaptive grid refinement with hanging nodes.
/**
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

class HangingNodeRefiner_IR1 : public GridObserver
{
	public:
		HangingNodeRefiner_IR1(IRefinementCallback* refCallback = NULL);
		HangingNodeRefiner_IR1(Grid& grid, IRefinementCallback* refCallback = NULL);
	//todo: make copy-constructor public
		virtual ~HangingNodeRefiner_IR1();

		virtual void grid_to_be_destroyed(Grid* grid);

		void assign_grid(Grid& grid);
		void set_refinement_callback(IRefinementCallback* refCallback);
		
		void clear_marks();
		void mark_for_refinement(EdgeBase* e);
		void mark_for_refinement(Face* f);
		void mark_for_refinement(Volume* v);

	///	the value-type of TIterator has to be a pointer to a type derived from either EdgeBase, Face or Volume.
		template <class TIterator>
		void mark_for_refinement(const TIterator& iterBegin, const TIterator& iterEnd);

	///	performs refinement on the marked elements.
	/**
	 * automatically extends the refinement to avoid multiple hanging nodes
	 * on a single edge or face.
	 */
		void refine();

	protected:
		typedef Attachment<int> ARefinementMark;
		//typedef Attachment<RefinementInfo> ARefinementInfo;

		typedef std::vector<EdgeBase*>	EdgeVec;
		typedef std::vector<Face*>		FaceVec;
		typedef std::vector<Volume*>	VolumeVec;

		typedef std::queue<EdgeBase*> EdgeQueue;
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

	////////////////////////////////////////////////////////////////////////
	//	refine methods
		void refine_constraining_edge(ConstrainingEdge* cge);
		void refine_edge_with_normal_vertex(EdgeBase* e);
		void refine_edge_with_hanging_vertex(EdgeBase* e);

		void refine_face_with_normal_vertex(Face* f);
		void refine_face_with_hanging_vertex(Face* f);
		void refine_constraining_face(ConstrainingFace* cgf);
		
		void refine_volume_with_normal_vertex(Volume* v);
		
	////////////////////////////////////////////////////////////////////////
	//	helpers. Make sure that everything is initialized properly
	//	before calling these methods.
	//	you should use this methods instead of directly marking elements.
		inline bool is_marked(EdgeBase* e)							{return m_selMarkedElements.is_selected(e);}
		inline void mark(EdgeBase* e)								{m_selMarkedElements.select(e);}
		inline VertexBase* get_center_vertex(EdgeBase* e)			{return m_aaVertexEDGE[e];}
		inline void set_center_vertex(EdgeBase* e, VertexBase* v)	{m_aaVertexEDGE[e] = v;}

		inline bool is_marked(Face* f)								{return m_selMarkedElements.is_selected(f);}
		inline void mark(Face* f)									{m_selMarkedElements.select(f);}
		inline VertexBase* get_center_vertex(Face* f)				{return m_aaVertexFACE[f];}
		inline void set_center_vertex(Face* f, VertexBase* v)		{m_aaVertexFACE[f] = v;}

		inline bool is_marked(Volume* v)							{return m_selMarkedElements.is_selected(v);}
		inline void mark(Volume* v)									{m_selMarkedElements.select(v);}

	/**	used during collect_objects_for_refine.
	 *	unmarked associated elements of the elements between elemsBegin and
	 *	elemsEnd are pushed to the queue.
	 *	Please note that a single object may appear multiple times in the queue.
	 *	\{ */
		template <class TIterator>
		void collect_associated_unmarked_edges(
							std::queue<EdgeBase*>& qEdgesOut, Grid& grid,
							TIterator elemsBegin, TIterator elemsEnd);
		
		template <class TIterator>
		void collect_associated_unmarked_faces(
							std::queue<Face*>& qFacesOut, Grid& grid,
							TIterator elemsBegin, TIterator elemsEnd);
							
		template <class TIterator>
		void collect_associated_unmarked_volumes(
							std::queue<Volume*>& qVolsOut, Grid& grid,
							TIterator elemsBegin, TIterator elemsEnd);
	/** \} */
	
	private:
		HangingNodeRefiner_IR1(const HangingNodeRefiner_IR1& hnr);

	protected:
		Grid*		m_pGrid;

		AVertexBase		m_aVertex;

		Selector	m_selMarkedElements;

		IRefinementCallback*	m_refCallback;

		Grid::EdgeAttachmentAccessor<AVertexBase>		m_aaVertexEDGE;
		Grid::FaceAttachmentAccessor<AVertexBase>		m_aaVertexFACE;
};


/// @}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of template methods.
template <class TIterator>
void HangingNodeRefiner_IR1::mark_for_refinement(const TIterator& iterBegin,
												const TIterator& iterEnd)
{
	TIterator iter = iterBegin;
	while(iter != iterEnd)
	{
		mark_for_refinement(*iter);
		++iter;
	}
}

}// end of namespace

#endif
