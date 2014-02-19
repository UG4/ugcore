// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.01.2012 (m,d,y)

#ifndef __H__UG__bool_marker__
#define __H__UG__bool_marker__

#include "lib_grid/grid/grid.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

/** \ingroup lib_grid_tools
 *  \{ */

///	Allows to mark elements.
/** This class allows to mark elements of a grid.
 * The BoolMarker associates a bool with each element.
 * Note that clearing the marks has a runtime complexity of O(n). If you need
 * marks for repeatedly called local algorithms you may want to use Grid::mark
 * instead, which has a clear_marks method with runtime complexity of O(1).
 *
 * Note that methods like mark, unmark, is_marked, clear, ... may only be invoked,
 * if a grid was assigned through either assign_grid or through the constructor.
 *
 * Marks can be passed to children, when an element is refined. You have to enable
 * mark inheritance in order to activate this behavior (use enable_mark_inheritance).
 * Mark inheritance is enabled by default.
 *
 * \todo	Allow to restrict marking to vertices, edges, faces or volumes
 * \todo	Add is_marked, mark, unmark for GridObject
 * \todo	Refactor to template \<class T\> Marker.
 */
class BoolMarker : public GridObserver
{
	public:
		BoolMarker();
		BoolMarker(Grid& g);

		virtual ~BoolMarker();

	///	Assign the grid on which the marker shall operate.
	/**	NULL is a valid argument and sets the marker into an unassigned state.
	 * The marker may only be used, if it is associated with a grid instance.*/
		void assign_grid(Grid* g);
	///	Assign the grid on which the marker shall operate.
		void assign_grid(Grid& g)					{assign_grid(&g);}

		Grid* grid()								{return m_pGrid;}

	///	set the mark which is applied when a new element is created
	/**	By default the default-mark is set to false. Note that the default mark
	 * only has an effect on child elements, if mark-inheritance is disabled.
	 * @param mark	this mark is set to new elements on creation.*/
		void set_default_mark(bool mark)			{m_defaultMark = mark;}
	///	returns the default mark.
		bool default_mark()							{return m_defaultMark;}

	///	if enabled, marks are passed from parents on to their children
	/**	\{ */
		void enable_mark_inheritance(bool enable)	{m_markInheritanceEnabled = enable;}
		bool mark_inheritance_enabeld()				{return m_markInheritanceEnabled;}
	/**	\} */

	/**	restricts mark inheritance so that new elements derive their selection
	 * status only from parents with the same base-type. Disabled by default.
	 * 	NOTE: strict inheritance only has an effect if selection inheritance is enabled.
	 * 	\{ */
		void enable_strict_inheritance(bool enable)	{m_strictInheritanceEnabled = enable;}
		bool strict_inheritance_enabled()			{return m_strictInheritanceEnabled;}
	/**	\} */


		bool is_marked(GridObject* e) const;
		bool is_marked(Vertex* e) const			{assert(m_pGrid); return m_aaMarkVRT[e];}
		bool is_marked(EdgeBase* e)	const			{assert(m_pGrid); return m_aaMarkEDGE[e];}
		bool is_marked(Face* e)	const				{assert(m_pGrid); return m_aaMarkFACE[e];}
		bool is_marked(Volume* e) const				{assert(m_pGrid); return m_aaMarkVOL[e];}

		void mark(Vertex* e, bool mark = true)	{assert(m_pGrid); m_aaMarkVRT[e] = mark;}
		void mark(EdgeBase* e, bool mark = true)	{assert(m_pGrid); m_aaMarkEDGE[e] = mark;}
		void mark(Face* e, bool mark = true)		{assert(m_pGrid); m_aaMarkFACE[e] = mark;}
		void mark(Volume* e, bool mark = true)		{assert(m_pGrid); m_aaMarkVOL[e] = mark;}

		template <class TIter>
		void mark(TIter begin, TIter end, bool mark = true)
		{
			for(TIter iter = begin; iter != end; ++iter) BoolMarker::mark(*iter, mark);
		}

		void unmark(Vertex* e)	{mark(e, false);}
		void unmark(EdgeBase* e)	{mark(e, false);}
		void unmark(Face* e)		{mark(e, false);}
		void unmark(Volume* e)		{mark(e, false);}

		template <class TIter>
		void unmark(TIter begin, TIter end)			{mark(begin, end, false);}

	///	Sets all marks to false. O(n).
		void clear();

	///	derived from GridObserver
		virtual void grid_to_be_destroyed(Grid* grid);

	//	element callbacks
		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void edge_created(Grid* grid, EdgeBase* e,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void volume_created(Grid* grid, Volume* vol,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void vertices_to_be_merged(Grid* grid, Vertex* target,
										 Vertex* elem1, Vertex* elem2);

		virtual void edges_to_be_merged(Grid* grid, EdgeBase* target,
										 EdgeBase* elem1, EdgeBase* elem2);

		virtual void faces_to_be_merged(Grid* grid, Face* target,
										 Face* elem1, Face* elem2);

		virtual void volumes_to_be_merged(Grid* grid, Volume* target,
										 Volume* elem1, Volume* elem2);

	protected:
		Grid*	m_pGrid;
		ABool	m_aBool;
		bool	m_defaultMark;
		bool	m_markInheritanceEnabled;
		bool	m_strictInheritanceEnabled;
		Grid::AttachmentAccessor<Vertex, ABool>	m_aaMarkVRT;
		Grid::AttachmentAccessor<EdgeBase, ABool>	m_aaMarkEDGE;
		Grid::AttachmentAccessor<Face, ABool>		m_aaMarkFACE;
		Grid::AttachmentAccessor<Volume, ABool>		m_aaMarkVOL;
};

/** \} */

}//	end of namespace

#endif
