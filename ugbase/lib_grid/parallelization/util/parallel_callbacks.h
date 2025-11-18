/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__parallel_callbacks__
#define __H__UG__parallel_callbacks__

#include "../distributed_grid.h"
#include "common/assert.h"

namespace ug
{

///	Returns true if an element is a regular surface element.
/**	Regular surface elements are elements which lie on the
 * surface and are not a ghost element.
 * (Ghost are vertical masters, which do not lie in any other interface).
 */
class IsRegularSurfaceElem
{
	public:
		IsRegularSurfaceElem(const DistributedGridManager& dgm) :
			m_dgm(dgm), m_mg(dgm.get_assigned_grid())
		{
			UG_ASSERT(m_mg, "A grid has to be assigned to the distributed grid manager.");
		}

		bool operator () (Vertex* v) {return is_ok(v);}
		bool operator () (Edge* e) {return is_ok(e);}
		bool operator () (Face* f) {return is_ok(f);}
		bool operator () (Volume* v) {return is_ok(v);}

	private:
		template <class TElem>
		inline bool is_ok(TElem* e)
		{
			return !(m_mg->has_children(e) || m_dgm.is_ghost(e));
		}

		const DistributedGridManager& 	m_dgm;
		const MultiGrid*				m_mg;
};

}//	end of namespace

#endif
