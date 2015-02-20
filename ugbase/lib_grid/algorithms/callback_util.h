// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.08.2012 (m,d,y)

#ifndef __H__UG__callback_util__
#define __H__UG__callback_util__

#include "lib_grid/tools/selector_interface.h"
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug
{

/**
 * \brief contains grid element callback definitions for common usage
 *
 * \defgroup lib_grid_algorithms_element_callbacks element_callbacks
 * \ingroup lib_grid_algorithms
 * \{
 */

class ConsiderAll{
	public:
		bool operator() (Vertex* v) const	{return true;}
		bool operator() (Edge* e) const		{return true;}
		bool operator() (Face* f) const		{return true;}
		bool operator() (Volume* v) const	{return true;}
};


class ConsiderNone{
	public:
		bool operator() (Vertex* v) const	{return false;}
		bool operator() (Edge* e) const		{return false;}
		bool operator() (Face* f) const		{return false;}
		bool operator() (Volume* v) const	{return false;}
};

////////////////////////////////////////////////////////////////////////////////
///	Element callback that returns true, if an element is marked
class IsMarked
{
	public:
		IsMarked(const Grid& grid) :
			m_grid(grid)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)	{return callback(f);}
		bool operator() (Volume* v)	{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return m_grid.is_marked(e);}

	private:
		const Grid&	m_grid;
};

///	Element callback that returns true, if an element is not marked
class IsNotMarked
{
	public:
		IsNotMarked(const Grid& grid) :
			m_grid(grid)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)	{return callback(f);}
		bool operator() (Volume* v)	{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return !m_grid.is_marked(e);}

	private:
		const Grid&	m_grid;
};


////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
///	Element callback that returns true, if an element is contained in a subset
class IsInSubset
{
	public:
		IsInSubset(const ISubsetHandler& sh, int subsetIndex) :
			m_sh(sh),
			m_si(subsetIndex)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}
		bool operator() (Volume* v)		{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return m_sh.get_subset_index(e) == m_si;}

	private:
		const ISubsetHandler& m_sh;
		int m_si;
};

///	Element callback that returns true, if an element is not contained in a subset
class IsNotInSubset
{
	public:
		IsNotInSubset(const ISubsetHandler& sh, int subsetIndex) :
			m_sh(sh),
			m_si(subsetIndex)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}
		bool operator() (Volume* v)		{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return m_sh.get_subset_index(e) != m_si;}

	private:
		const ISubsetHandler& m_sh;
		int m_si;
};


////////////////////////////////////////////////////////////////////////////////
///	Element callback that returns true, if an element lies on the grids boundary
class IsOnBoundary
{
	public:
		IsOnBoundary(Grid& g) :
			m_grid(g)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return LiesOnBoundary(m_grid, e);}

	private:
		Grid&	m_grid;
};

///	Element callback that returns true, if an element does not lie on the grids boundary
class IsNotOnBoundary
{
	public:
		IsNotOnBoundary(Grid& g) :
			m_grid(g)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return !LiesOnBoundary(m_grid, e);}

	private:
		Grid&	m_grid;
};

/** \} */ // lib_grid_algorithms_callbacks
}//	end of namespace

#endif
