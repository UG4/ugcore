// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011

#ifndef __H__UG__HANGING_NODE_REFINER_BASE__
#define __H__UG__HANGING_NODE_REFINER_BASE__

#include <queue>
#include <vector>
#include "lib_grid/lg_base.h"
#include "refinement_callbacks.h"
#include "refiner_interface.h"
#include "ref_mark_adjuster_interface.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	Base class for a hanging-node refiner.
/**	A hanging node refiner allows to adaptively refine grid elements
 * by inserting hanging nodes at points where T-junctions would occur
 * (ug::HangingVertex).
 * If further refinement is performed, the refiner automatically
 * takes care to adjust the refined area, so that no more than
 * 1 hanging vertex resides on any edge.
 *
 * Each edge on which a hanging node lies is replaced by a
 * ug::ConstrainingEdge. Child-edges of a constraining edge have the
 * type ug::ConstrainedEdge.
 *
 * If volumes exist, the refiner may also create elements of type
 * ug::ConstrainingTriangle or ug::ConstrainingQuadrilateral, as well
 * as ug::ConstrainedTriangle and ug::ConstrainedQuadrilateral.
 *
 * Use the mark methods to mark elements which
 * shall be refined. A call to refine will then perform the refinement
 * and clear all marks.
 *
 * Please note: If you're using a hanging node refiner, you have to
 * be careful to not destroy the connections between hanging-vertices,
 * and constrained / constraining objects.
 *
 * Specializations of this class exist to support hanging node refinement
 * on flat and on hierarchical grids.
 *
 * Note that you may set a refinement callback, which will be
 * responsible to calculate new positions of newly created vertices.
 * By default a linear refinement callback is created for one of
 * the standard position attachments (ug::aPosition, ug::aPosition2,
 * ug::aPosition1 - whichever is present).
 *
 * This class can't be instantiated directly. Use one of its
 * specializations instead.
 *
 * For Developers: Note that HangingNodeRefinerBase stores flags together with
 * the refinement marks. Precisely the value 128 is used to flag whether an
 * element has to be refined using hanging-node-rules.
 *
 * \sa ug::HangingNodeRefiner_Grid, ug::HangingNodeRefiner_MultiGrid
 */
template <class TSelector>
class HangingNodeRefinerBase : public IRefiner, public GridObserver
{
	public:
		using IRefiner::mark;

		typedef TSelector	selector_t;

	/**	additional mark to RefinementMarks. Used to flag whether an element
	 * will be refined with constraining.*/
		enum HNodeRefMarks{
			HNRM_TO_NORMAL = 		1 << 5,
			HNRM_TO_CONSTRAINED =	1 << 6,
			HNRM_TO_CONSTRAINING =	1 << 7,
			HNRM_MAX
		};

	public:
		HangingNodeRefinerBase(IRefinementCallback* refCallback = NULL);
		virtual ~HangingNodeRefinerBase();

		virtual void grid_to_be_destroyed(Grid* grid);

	///	enables or disables node-dependency-order-1.
	/**	\{
	 * If enabled, hanging nodes may only depend on non-hanging nodes.
	 * An edge containing a hanging node thus will not have a hanging-node
	 * as a corner vertex.
	 *
	 * Enabled by default.*/
		void enable_node_dependency_order_1(bool bEnable);
		bool node_dependency_order_1_enabled()				{return m_nodeDependencyOrder1;}
	/**	\} */

		virtual void clear_marks();

	///	Marks a element for refinement.
	/**	\{ */
		virtual bool mark(Vertex* v, RefinementMark refMark = RM_REFINE);
		virtual bool mark(Edge* e, RefinementMark refMark = RM_REFINE);
		virtual bool mark(Face* f, RefinementMark refMark = RM_REFINE);
		virtual bool mark(Volume* v, RefinementMark refMark = RM_REFINE);
	/**	\} */

	///	Marks the neighborhood of the current selection.
	/**	\sa ISelector::mark_neighborhood*/
		virtual void mark_neighborhood(size_t numIterations, RefinementMark refMark);

	///	Returns the mark of a given element.
	/**	\{ */
		virtual RefinementMark get_mark(Vertex* v);
		virtual RefinementMark get_mark(Edge* e);
		virtual RefinementMark get_mark(Face* f);
		virtual RefinementMark get_mark(Volume* v);
	/**	\} */

		virtual bool save_marks_to_file(const char* filename);

	///	Add a refmark adjuster, which will be called while marks are adjusted during refinement / coarsening
		void add_ref_mark_adjuster(SPIRefMarkAdjuster adjuster)		{m_refMarkAdjusters.push_back(adjuster);}

	protected:
		typedef typename TSelector::template traits<Vertex>::iterator	sel_vrt_iter;
		typedef typename TSelector::template traits<Edge>::iterator		sel_edge_iter;
		typedef typename TSelector::template traits<Face>::iterator			sel_face_iter;
		typedef typename TSelector::template traits<Volume>::iterator		sel_vol_iter;

	///	performs refinement on the marked elements.
	/**
	 * The grid's message hub is informed using a "GridAdaption" message,
	 * passing an instance of GridMessage_Adapation, with values
	 * GMAT_HNODE_REFINEMENT_BEGINS and GMAT_HNODE_REFINEMENT_ENDS.
	 * See lib_grid/lib_grid_messages.h for more details.
	 *
	 * as defined in lib_grid/
	 *
	 * automatically extends the refinement to avoid multiple hanging nodes
	 * on a single edge or face.
	 *
	 * refine calls several virtual methods, which allow to influence the
	 * refinement process. Most notably the methods
	 *
	 *		- collect_objects_for_refine
	 *		- pre_refine
	 * 		- post_refine
	 *
	 * are called in the given order. During element refinement further
	 * virtual methods are called, which perform the actual element refinement.
	 */
		void perform_refinement();


	///	a callback that allows to deny refinement of special vertices
		virtual bool refinement_is_allowed(Vertex* elem)	{return true;}
	///	a callback that allows to deny refinement of special edges
		virtual bool refinement_is_allowed(Edge* elem)		{return true;}
	///	a callback that allows to deny refinement of special faces
		virtual bool refinement_is_allowed(Face* elem)			{return true;}
	///	a callback that allows to deny refinement of special volumes
		virtual bool refinement_is_allowed(Volume* elem)		{return true;}

	///	performs registration and deregistration at a grid.
	/**	Sets a grid and performs registration at the given grid.
	 * 	The associated selector is also initialised with the given grid.
	 * 	It is cruical to call this method or everything will fail.
	 *
	 *  call set_grid(NULL) to unregister the observer from a grid.
	 *
	 *  Please note that this method is not declared virtual, since it
	 *  is called during construction and destruction.*/
		void set_grid(typename TSelector::grid_type* grid);

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

	/**	after each iteration in collet_objects_for_refine, this method determines
	 * whether the iteration shall be continued. Important for parallel refiners.
	 * The default implementation simply returns the specified value. This is fine
	 * for serial environments.*/
		virtual bool continue_collect_objects_for_refine(bool continueRequired)
		{return continueRequired;}

	/**	This callback is called during execution of the refine() method after
	 * collect_objects_for_refine has returned. It is responsible to mark
	 * elements for hnode refinement. That means all elements on which a hanging
	 * node or constrained children shall be created have to be marked using
	 * mark_for_hnode_refinement during this method. The default implementation
	 * performs this marking for all local elements.
	 * Make sure to not mark any new elements for refinement during this method.*/
		virtual void assign_hnode_marks();

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
	 *  \{ */
		virtual void process_constrained_vertex(ConstrainedVertex* cdv);
		virtual void process_constrained_edge(ConstrainedEdge* cde);
		virtual void process_constraining_edge(ConstrainingEdge* cge);
		virtual void refine_edge_with_normal_vertex(Edge* e,
											Vertex** newCornerVrts = NULL);
		virtual void refine_edge_with_hanging_vertex(Edge* e,
											Vertex** newCornerVrts = NULL);

		virtual void process_constrained_face(ConstrainedFace* cdf);
		virtual void process_constraining_face(ConstrainingFace* cgf);
		virtual void refine_face_with_normal_vertex(Face* f,
											Vertex** newCornerVrts = NULL);
		virtual void refine_face_with_hanging_vertex(Face* f,
											Vertex** newCornerVrts = NULL);

		virtual void refine_volume_with_normal_vertex(Volume* v,
											Vertex** newVolumeVrts = NULL);
	/**	\} */

	////////////////////////////////////////////////////////////////////////
	//	helpers. Make sure that everything is initialized properly
	//	before calling these methods.
	//	you should use this methods instead of directly marking elements.
		inline bool is_marked(Vertex* v)				{return m_selMarkedElements.is_selected(v);}
		//inline void mark(Vertex* v)						{mark(v);}

		inline bool is_marked(Edge* e)					{return m_selMarkedElements.is_selected(e);}
		//inline void mark(Edge* e)						{mark(e);}

	///	Returns the vertex associated with the edge
	/**	pure virtual method.
	 *	Has to return the center vertex which was set to the edge via
	 *	set_center_vertex. If no vertex was set, NULL has to be returned.*/
		virtual Vertex* get_center_vertex(Edge* e) = 0;

	///	Associates a vertex with the edge (pure virtual).
		virtual void set_center_vertex(Edge* e, Vertex* v) = 0;

		inline bool is_marked(Face* f)						{return m_selMarkedElements.is_selected(f);}
		//inline void mark(Face* f)							{mark(f);}

	///	Returns the vertex associated with the face
	/**	pure virtual method.
	 *	Has to return the center vertex which was set to the face via
	 *	set_center_vertex. If no vertex was set, NULL has to be returned.*/
		virtual Vertex* get_center_vertex(Face* f) = 0;

	///	Associates a vertex with the face (pure virtual).
		virtual void set_center_vertex(Face* f, Vertex* v) = 0;

		inline bool is_marked(Volume* v)					{return m_selMarkedElements.is_selected(v);}
		//inline void mark(Volume* v)						{mark(v);}


		template <class TElem>
		inline bool marked_copy(TElem* elem)				{return (m_selMarkedElements.get_selection_status(elem) & RM_COPY) == RM_COPY;}

		template <class TElem>
		inline bool marked_refine(TElem* elem)				{return (m_selMarkedElements.get_selection_status(elem) & (RM_REFINE | RM_ANISOTROPIC)) != 0;}

		template <class TElem>
		inline bool marked_anisotropic(TElem* elem)			{return (m_selMarkedElements.get_selection_status(elem) & RM_ANISOTROPIC) == RM_ANISOTROPIC;}

		template <class TElem>
		inline bool marked_regular(TElem* elem)				{return marked_refine(elem) && (!(marked_anisotropic(elem) || marked_copy(elem)));}

		template <class TElem>
		inline bool marked_coarsen(TElem* elem)				{return (m_selMarkedElements.get_selection_status(elem) & RM_COARSEN) == RM_COARSEN;}

		template <class TElem>
		inline bool marked_to_normal(TElem* elem)			{return (m_selMarkedElements.get_selection_status(elem) & HNRM_TO_NORMAL) == HNRM_TO_NORMAL;}

		template <class TElem>
		inline bool marked_to_constrained(TElem* elem)		{return (m_selMarkedElements.get_selection_status(elem) & HNRM_TO_CONSTRAINED) == HNRM_TO_CONSTRAINED;}

		template <class TElem>
		inline bool marked_to_constraining(TElem* elem)		{return (m_selMarkedElements.get_selection_status(elem) & HNRM_TO_CONSTRAINING) == HNRM_TO_CONSTRAINING;}

		template <class TElem>
		void add_hmark(TElem* elem, HNodeRefMarks mark);

		template <class TElem>
		void remove_hmark(TElem* elem, uint mark);

	///	removes coarsen marks from the selection
	/**	Note that derived classes are not informed about those deselections!*/
		template <class TElem>
		bool remove_coarsen_marks();

	///	returns the selector which is internally used to mark elements.
	/**	Be sure to use it carefully!*/
		TSelector& get_refmark_selector()	{return m_selMarkedElements;}

	private:
	///	private copy constructor to avoid copy construction
		HangingNodeRefinerBase(const HangingNodeRefinerBase&);

	protected:
		TSelector							m_selMarkedElements;
		std::vector<SPIRefMarkAdjuster>		m_refMarkAdjusters;

	private:
		Grid*		m_pGrid;
		std::vector<Vertex*>	m_newlyMarkedRefVrts;
		std::vector<Edge*>		m_newlyMarkedRefEdges;
		std::vector<Face*>			m_newlyMarkedRefFaces;
		std::vector<Volume*>		m_newlyMarkedRefVols;
		//todo:	Use the following vectors during coarsening...
		/*
		std::vector<Vertex*>	m_newlyMarkedCoarseVrts;
		std::vector<Edge*>		m_newlyMarkedCoarseEdges;
		std::vector<Face*>			m_newlyMarkedCoarseFaces;
		std::vector<Volume*>		m_newlyMarkedCoarseVols;
		*/
		bool		m_nodeDependencyOrder1;
		//bool		m_automarkHigherDimensionalObjects; <-- unused
		bool		m_adjustingRefMarks;///<	true during collect_objects_for_refine
};

/// @}	// end of add_to_group command

}//	end of namespace

#endif
