// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.01.2012 (m,d,y)

#ifndef __H__UG__bool_marker__
#define __H__UG__bool_marker__

#include "lib_grid/grid/grid.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

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
 * \todo	Allow to restrict marking to vertices, edges, faces or volumes
 * \todo	Support mark-inheritance
 * \todo	Support default-marks
 * \todo	Add is_marked, mark, unmark for GeometricObject
 * \todo	Refactor to template <class T> Marker.
 */
class BoolMarker : public GridObserver
{
	public:
		BoolMarker();
		BoolMarker(Grid& g);

		virtual ~BoolMarker();

		void assign_grid(Grid* g);
		void assign_grid(Grid& g)					{assign_grid(&g);}

		Grid* grid()								{return m_pGrid;}

		bool is_marked(VertexBase* e) const			{assert(m_pGrid); return m_aaMarkVRT[e];}
		bool is_marked(EdgeBase* e)	const			{assert(m_pGrid); return m_aaMarkEDGE[e];}
		bool is_marked(Face* e)	const				{assert(m_pGrid); return m_aaMarkFACE[e];}
		bool is_marked(Volume* e) const				{assert(m_pGrid); return m_aaMarkVOL[e];}

		void mark(VertexBase* e, bool mark = true)	{assert(m_pGrid); m_aaMarkVRT[e] = mark;}
		void mark(EdgeBase* e, bool mark = true)	{assert(m_pGrid); m_aaMarkEDGE[e] = mark;}
		void mark(Face* e, bool mark = true)		{assert(m_pGrid); m_aaMarkFACE[e] = mark;}
		void mark(Volume* e, bool mark = true)		{assert(m_pGrid); m_aaMarkVOL[e] = mark;}

		template <class TIter>
		void mark(TIter begin, TIter end, bool mark = true)
		{
			for(TIter iter = begin; iter != end; ++iter) BoolMarker::mark(*iter, mark);
		}

		void unmark(VertexBase* e)	{mark(e, false);}
		void unmark(EdgeBase* e)	{mark(e, false);}
		void unmark(Face* e)		{mark(e, false);}
		void unmark(Volume* e)		{mark(e, false);}

		template <class TIter>
		void unmark(TIter begin, TIter end)			{mark(begin, end, false);}

	///	Sets all marks to false. O(n).
		void clear();

	///	derived from GridObserver
		virtual void grid_to_be_destroyed(Grid* grid);

	protected:
		Grid*	m_pGrid;
		ABool	m_aBool;
		Grid::AttachmentAccessor<VertexBase, ABool>	m_aaMarkVRT;
		Grid::AttachmentAccessor<EdgeBase, ABool>	m_aaMarkEDGE;
		Grid::AttachmentAccessor<Face, ABool>		m_aaMarkFACE;
		Grid::AttachmentAccessor<Volume, ABool>		m_aaMarkVOL;
};


}//	end of namespace

#endif
