// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_selection_callbacks
#define __H__UG_selection_callbacks

#include "lib_grid/tools/selector_interface.h"

namespace ug{
/** \ingroup lib_grid_element_callbacks
 * \{ */

///	Element callback that returns true, if an element is selected
class IsSelected
{
	public:
		IsSelected(const ISelector& sel) :
			m_sel(sel)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}
		bool operator() (Volume* v)		{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return m_sel.is_selected(e);}

	private:
		const ISelector&	m_sel;
};

///	Element callback that returns true, if an element is not selected
class IsNotSelected
{
	public:
		IsNotSelected(const ISelector& sel) :
			m_sel(sel)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}
		bool operator() (Volume* v)		{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return !m_sel.is_selected(e);}

	private:
		const ISelector&	m_sel;
};

/** \} */ //lib_grid_element_callbacks
}//	end of namespace

#endif	//__H__UG_selection_callbacks
