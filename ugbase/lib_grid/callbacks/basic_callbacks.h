// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_basic_callbacks
#define __H__UG_basic_callbacks

#include "lib_grid/grid/grid.h"

namespace ug{

/**
 * \brief contains grid element callback definitions for common usage
 *
 * \defgroup lib_grid_element_callbacks element_callbacks
 * \ingroup lib_grid
 * \{
 */

///	callback that always returns true
class ConsiderAll{
	public:
		bool operator() (Vertex* v) const	{return true;}
		bool operator() (Edge* e) const		{return true;}
		bool operator() (Face* f) const		{return true;}
		bool operator() (Volume* v) const	{return true;}
};

///	callback that always returns false
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

/** \} */ // lib_grid_element_callbacks

}//	end of namespace

#endif	//__H__UG_basic_callbacks
