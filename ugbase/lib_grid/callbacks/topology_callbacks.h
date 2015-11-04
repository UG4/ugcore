// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_topology_callbacks
#define __H__UG_topology_callbacks

#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

namespace ug{

/** \ingroup lib_grid_element_callbacks
 * \{ */

///	Element callback that returns true, if an element lies on the grids boundary
class IsOnBoundary
{
	public:
		IsOnBoundary(Grid& g) :
			m_grid(g)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)	{return callback(f);}

	private:
		template <class TElem>
		bool callback(TElem* e)		{return LiesOnBoundary(m_grid, e);}

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
		bool operator() (Face* f)	{return callback(f);}

	private:
		template <class TElem>
		bool callback(TElem* e)		{return !LiesOnBoundary(m_grid, e);}

	private:
		Grid&	m_grid;
};

/** \} */ //lib_grid_element_callbacks

}//	end of namespace

#endif	//__H__UG_topology_callbacks
