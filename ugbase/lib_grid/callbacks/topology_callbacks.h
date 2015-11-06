/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
