// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SELECTOR_CALLBACKS__
#define __H__LIB_GRID__SELECTOR_CALLBACKS__

#include "lib_grid/tools/selector_interface.h"
#include "callback_definitions.h"

namespace ug
{

class IsSelected
{
	public:
		IsSelected(const Selector& sel) :
			m_sel(sel)	{}

		bool operator() (VertexBase* v)	{return m_sel.is_selected(v);}
		bool operator() (EdgeBase* e)	{return m_sel.is_selected(e);}
		bool operator() (Face* f)		{return m_sel.is_selected(f);}
		bool operator() (Volume* v)		{return m_sel.is_selected(v);}

	private:
		const Selector&		m_sel;
};

}//	end of namespace

#endif
