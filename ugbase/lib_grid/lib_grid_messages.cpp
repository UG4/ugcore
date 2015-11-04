// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 02.03.2012 (m,d,y)
 
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
