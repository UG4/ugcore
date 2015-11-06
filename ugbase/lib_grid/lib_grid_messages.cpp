/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "lib_grid_messages.h"

namespace ug{


bool GridMessage_Adaption::adaption_begins() const
{
	return	(m_adaptionType == GMAT_GLOBAL_ADAPTION_BEGINS)
			|| (m_adaptionType == GMAT_HNODE_ADAPTION_BEGINS);
}

bool GridMessage_Adaption::adaption_ends() const
{
	return	(m_adaptionType == GMAT_GLOBAL_ADAPTION_ENDS)
			|| (m_adaptionType == GMAT_HNODE_ADAPTION_ENDS);
}

bool GridMessage_Adaption::step_begins() const
{
	switch(m_adaptionType){
		case GMAT_GLOBAL_REFINEMENT_BEGINS:
		case GMAT_HNODE_REFINEMENT_BEGINS:
		case GMAT_GLOBAL_COARSENING_BEGINS:
		case GMAT_HNODE_COARSENING_BEGINS:
			return true;
		default:
			return false;
	}
}

bool GridMessage_Adaption::step_ends() const
{
	switch(m_adaptionType){
		case GMAT_GLOBAL_REFINEMENT_ENDS:
		case GMAT_HNODE_REFINEMENT_ENDS:
		case GMAT_GLOBAL_COARSENING_ENDS:
		case GMAT_HNODE_COARSENING_ENDS:
			return true;
		default:
			return false;
	}
}

bool GridMessage_Adaption::adaptive() const
{
	switch(m_adaptionType){
		case GMAT_HNODE_ADAPTION_BEGINS:
		case GMAT_HNODE_ADAPTION_ENDS:
		case GMAT_HNODE_REFINEMENT_BEGINS:
		case GMAT_HNODE_REFINEMENT_ENDS:
		case GMAT_HNODE_COARSENING_BEGINS:
		case GMAT_HNODE_COARSENING_ENDS:
			return true;
		default:
			return false;
	}
}

bool GridMessage_Adaption::global() const
{
	switch(m_adaptionType){
		case GMAT_GLOBAL_ADAPTION_BEGINS:
		case GMAT_GLOBAL_ADAPTION_ENDS:
		case GMAT_GLOBAL_REFINEMENT_BEGINS:
		case GMAT_GLOBAL_REFINEMENT_ENDS:
		case GMAT_GLOBAL_COARSENING_BEGINS:
		case GMAT_GLOBAL_COARSENING_ENDS:
			return true;
		default:
			return false;
	}
}

bool GridMessage_Adaption::refinement() const
{
	switch(m_adaptionType){
		case GMAT_GLOBAL_REFINEMENT_BEGINS:
		case GMAT_HNODE_REFINEMENT_BEGINS:
		case GMAT_GLOBAL_REFINEMENT_ENDS:
		case GMAT_HNODE_REFINEMENT_ENDS:
			return true;
		default:
			return false;
	}
}

bool GridMessage_Adaption::coarsening() const
{
	switch(m_adaptionType){
		case GMAT_GLOBAL_COARSENING_BEGINS:
		case GMAT_HNODE_COARSENING_BEGINS:
		case GMAT_GLOBAL_COARSENING_ENDS:
		case GMAT_HNODE_COARSENING_ENDS:
			return true;
		default:
			return false;
	}
}

}// end of namespace
