// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011

#ifndef __H__UG__HANGING_NODE_REFINER_BASE__
#define __H__UG__HANGING_NODE_REFINER_BASE__

#include <queue>
#include "lib_grid/lg_base.h"
#include "refinement_callbacks.h"

namespace ug
{

class HangingNodeRefinerBase : public GridObserver
{
	public:
		HangingNodeRefinerBase(IRefinementCallback* refCallback = NULL);
		virtual ~HangingNodeRefinerBase();

		virtual void grid_to_be_destroyed(Grid* grid);

		void set_refinement_callback(IRefinementCallback* refCallback);
		IRefinementCallback* get_refinement_callback()	{return m_refCallback;}

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
	 *
	 * refine calls several virtual methods, which allow to influence the
	 * refinement process. Most notably the methods
	 *
	 * 	- collect_objects_for_refine
	 * 	- pre_refine
	 *  - post_refine
	 *
	 * are called in the given order. During element refinement further
	 * virtual methods are called, which perform the actual element refinement.
	 */
		void refine();

	protected:
	///	performs registration and deregistration at a grid.
	/**	Sets a grid and performs registration at the given grid.
	 * 	The associated selector is also initialised with the given grid.
	 * 	It is cruical to call this method or everything will fail.
	 *
	 *  call set_grid(NULL) to unregister the observer from a grid.
	 *
	 *  Please note that this method is not declared virtual, since it
	 *  is called during construction and destruction.*/
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

	/**	called by refine after collect_objects_for_refine and before
	 *	actual refinement begins.*/
		virtual void pre_refine()	{}

	/**	called by refine after refinement is done.*/
		virtual void post_refine()	{}

	////////////////////////////////////////////////////////////////////////
	//	refine methods
	///	called to refine the specified element.
	/**	Refines the element. Corner vertices of the newly created element
	 *  can be specified through newCornerVrts. newCornerVrts = NULL (default)
	 *  means, that the corner vertices of the original element shall be taken.
	 *  \{
	 */
		virtual void refine_constraining_edge(ConstrainingEdge* cge);
		virtual void refine_edge_with_normal_vertex(EdgeBase* e,
											VertexBase** newCornerVrts = NULL);
		virtual void refine_edge_with_hanging_vertex(EdgeBase* e,
											VertexBase** newCornerVrts = NULL);

		virtual void refine_face_with_normal_vertex(Face* f,
											VertexBase** newCornerVrts = NULL);
		virtual void refine_face_with_hanging_vertex(Face* f,
											VertexBase** newCornerVrts = NULL);
		virtual void refine_constraining_face(ConstrainingFace* cgf);

		virtual void refine_volume_with_normal_vertex(Volume* v,
											VertexBase** newVolumeVrts = NULL);
	/**	\} */

	////////////////////////////////////////////////////////////////////////
	//	helpers. Make sure that everything is initialized properly
	//	before calling these methods.
	//	you should use this methods instead of directly marking elements.
		inline bool is_marked(EdgeBase* e)					{return m_selMarkedElements.is_selected(e);}
		virtual void mark(EdgeBase* e)						{m_selMarkedElements.select(e);}

	///	Returns the vertex associated with the edge
	/**	pure virtual method.
	 *	Has to return the center vertex which was set to the edge via
	 *	set_center_vertex. If no vertex was set, NULL has to be returned.*/
		virtual VertexBase* get_center_vertex(EdgeBase* e) = 0;

	///	Associates a vertex with the edge (pure virtual).
		virtual void set_center_vertex(EdgeBase* e, VertexBase* v) = 0;

		inline bool is_marked(Face* f)						{return m_selMarkedElements.is_selected(f);}
		virtual void mark(Face* f)							{m_selMarkedElements.select(f);}

	///	Returns the vertex associated with the face
	/**	pure virtual method.
	 *	Has to return the center vertex which was set to the face via
	 *	set_center_vertex. If no vertex was set, NULL has to be returned.*/
		virtual VertexBase* get_center_vertex(Face* f) = 0;

	///	Associates a vertex with the face (pure virtual).
		virtual void set_center_vertex(Face* f, VertexBase* v) = 0;

		inline bool is_marked(Volume* v)							{return m_selMarkedElements.is_selected(v);}
		virtual void mark(Volume* v)								{m_selMarkedElements.select(v);}

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
	///	private copy constructor to avoid copy construction
		HangingNodeRefinerBase(const HangingNodeRefinerBase&);

	protected:
		Selector	m_selMarkedElements;
		IRefinementCallback*	m_refCallback;

	private:
		Grid*		m_pGrid;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of template methods.
template <class TIterator>
void HangingNodeRefinerBase::mark_for_refinement(const TIterator& iterBegin,
												const TIterator& iterEnd)
{
	TIterator iter = iterBegin;
	while(iter != iterEnd)
	{
		mark_for_refinement(*iter);
		++iter;
	}
}

}//	end of namespace

#endif
