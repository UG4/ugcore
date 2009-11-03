//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d09

#ifndef __H__LIB_GRID__MULTI_GRID_REFINER__
#define __H__LIB_GRID__MULTI_GRID_REFINER__

#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"

namespace ug
{

class MultiGridRefiner : public GridObserver
{
	public:
		MultiGridRefiner();
		MultiGridRefiner(MultiGrid& mg);

		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);

		void assign_grid(MultiGrid& mg);

		void clear_marks();

		template <class TElem>
		inline void mark_for_refinement(TElem* elem)	{m_selMarks.select(elem);}

	///	the value-type of TIterator has to be a pointer to a type derived from either EdgeBase, Face or Volume.
		template <class TIterator>
		inline void mark_for_refinement(const TIterator& iterBegin,
										const TIterator& iterEnd)
					{m_selMarks.select(iterBegin, iterEnd);}

		template <class TElem>
		inline bool is_marked(TElem* elem)	{return m_selMarks.is_selected(elem);}

	///	performs refinement on the marked elements.
		void refine();

	protected:
		virtual void collect_objects_for_refine();

	protected:
		MultiGrid*	m_pMG;
		Selector	m_selMarks;
};

}//	end of namespace

#endif
