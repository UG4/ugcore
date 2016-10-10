/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

#include "domain.h"
#include "common/util/table.h"


using namespace std;

namespace ug{

std::string DomainInfo::
to_string() const
{
	if((m_numElems.size() != m_numLocalElems.size())
	   || (m_numElems.size() != m_numLocalGhosts.size()))
	{
		UG_THROW("elem-arrays have to have the same number of entries!");
	}

	StringStreamTable t;

	t(0, 0) << "lvl";
	t(0, 1) << "#total-elems";
	t(0, 2) << "#local-elems";
	t(0, 3) << "(% of total)";
	t(0, 4) << "#local-ghosts";
	t(0, 5) << "#min-local-elems";
	t(0, 6) << "#max-local-elems";

	for(size_t i = 0; i < m_numElems.size(); ++i){
		int r = i+1;
		t(r, 0) << i;
		t(r, 1) << m_numElems[i];
		t(r, 2) << m_numLocalElems[i];
		if(m_numElems[i] > 0)
			t(r, 3) << (float)m_numLocalElems[i] / (float)m_numElems[i];
		else
			t(r, 3) << "-";
		t(r, 4) << m_numLocalGhosts[i];
		t(r, 5) << m_minNumLocalElems[i];
		t(r, 6) << m_maxNumLocalElems[i];
	}

	return t.to_string();
}

}//	end of namespace
