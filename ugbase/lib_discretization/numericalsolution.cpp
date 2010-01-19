/*
 * numericalsolution.cpp
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#include "numericalsolution.h"

namespace ug{

/////////////////////////////////////
/// pattern
/////////////////////////////////////

DoFPattern::DoFPattern(std::string name, ISubsetHandler& sh)
{
	GridSubsetHandler* grid_sh = dynamic_cast<GridSubsetHandler*>(&sh);
	MultiGridSubsetHandler* mg_sh = dynamic_cast<MultiGridSubsetHandler*>(&sh);

	_name = name;
	_sh = &sh;
	if(grid_sh != NULL)
	{
		_num_sh = grid_sh->num_subsets();
	}
	else if(mg_sh != NULL)
	{
		_num_sh = mg_sh->num_subsets();
	}
	else { assert(0); }

	_lock = false;
};

bool DoFPattern::add_solution(std::string name, TrialSpaceType TrialSpace, std::vector<int>& SubsetIndices)
{
	if(_lock == true)
	{
		std::cout << "Pattern already locked. Can not add solution" << std::endl;
		return false;
	}

	// Create SingleSolutionInfo
	SingleSolutionInfo* info = new SingleSolutionInfo();
	info->name = name;
	info->TrialSpace = TrialSpace;
	(info->SubsetIndex).resize(SubsetIndices.size());
	(info->group_comp).resize(SubsetIndices.size());
	for(uint i = 0; i < SubsetIndices.size(); ++i)
	{
		(info->SubsetIndex)[i] = SubsetIndices[i];
		(info->group_comp)[i] = 0;
	}
	_SingleSolutionInfoVec.push_back(info);

	// Create Grouping with only on component: the SingleSolution
	GroupSolutionInfo* groupInfo = new GroupSolutionInfo();
	groupInfo->name = name;
	(groupInfo->SingleSolutions).resize(1);
	for(uint i = 0; i < 1; ++i)
	{
		(groupInfo->SingleSolutions)[i] = _SingleSolutionInfoVec.size() - 1;
	}
	(groupInfo->SubsetIndex).resize(SubsetIndices.size());
	(groupInfo->num_comp).resize(SubsetIndices.size());
	for(uint i = 0; i < SubsetIndices.size(); ++i)
	{
		(groupInfo->SubsetIndex)[i] = SubsetIndices[i];
		(groupInfo->num_comp)[i] = 1;
	}
	_GroupSolutionInfoVec.push_back(groupInfo);

	return true;
}

bool DoFPattern::group_solution(std::string name, std::vector<uint>& GroupSolutions)
{
	if(_lock == true)
	{
		std::cout << "Pattern already locked. Can not add solution" << std::endl;
		return false;
	}

	// Sort GroupSolutions to be grouped
	sort(GroupSolutions.begin(), GroupSolutions.end());

	// collect subsets and num_Comp of all solutions that are to be grouped
	std::vector<int> SubsetIndex;
	std::vector<int> num_comp;
	std::vector<int> SingleSolutions;

	// loop all groups, that will be merged to a new group
	for(uint i = 0; i < GroupSolutions.size(); ++i)
	{
		// check that group exist
		if(GroupSolutions[i] >= _GroupSolutionInfoVec.size() || GroupSolutions[i] < 0)
		{
			std::cout << "ERROR in group_solution: Trying to read Group " << GroupSolutions[i] << ", but grouping does not exist. Aborting" << std::endl;
			return false;
		}

		// remember current old group
		uint nr = GroupSolutions[i];

		// loop all subsets of current old group
		for(uint j = 0; j < (_GroupSolutionInfoVec[nr]->SubsetIndex).size(); ++j)
		{
			uint k;
			for(k = 0; k < SubsetIndex.size(); ++k)
			{
				if(SubsetIndex[k] == _GroupSolutionInfoVec[nr]->SubsetIndex[j])
				{
					num_comp[k] += _GroupSolutionInfoVec[nr]->num_comp[j];
					break;
				}
			}
			if(k == SubsetIndex.size())
			{
				SubsetIndex.push_back(_GroupSolutionInfoVec[nr]->SubsetIndex[j]);
				num_comp.push_back(_GroupSolutionInfoVec[nr]->num_comp[j]);
			}
		}
		for(uint j = 0; j < (_GroupSolutionInfoVec[nr]->SingleSolutions).size(); ++j)
		{
			for(uint k = 0; k < SingleSolutions.size(); ++k)
			{
				if(SingleSolutions[k] == _GroupSolutionInfoVec[nr]->SingleSolutions[j])
				{
					// This should never happen
					std::cout << "ERROR in group_solutions: Multiple occurance of a single solution in grouping. Internal error. Aborting" << std::endl;
					return false;
				}
			}
			SingleSolutions.push_back(_GroupSolutionInfoVec[nr]->SingleSolutions[j]);
		}
	}

	sort(SingleSolutions.begin(), SingleSolutions.end());

	// create new grouping
	GroupSolutionInfo* info = new GroupSolutionInfo();
	info->name = name;
	(info->SubsetIndex).resize(SubsetIndex.size());
	(info->num_comp).resize(SubsetIndex.size());
	for(uint i = 0; i < SubsetIndex.size(); ++i)
	{
		info->SubsetIndex[i] = SubsetIndex[i];
		info->num_comp[i] = num_comp[i];
	}

	(info->SingleSolutions).resize(SingleSolutions.size());
	for(uint i = 0; i < SingleSolutions.size(); ++i)
	{
		info->SingleSolutions[i] = SingleSolutions[i];
	}

	for(uint i = 0; i < SubsetIndex.size(); ++i)
	{
		num_comp[i] = 0;
		for(uint j = 0; j < SingleSolutions.size(); ++j)
		{
			uint nr = SingleSolutions[j];
			uint k;
			for(uint k = 0; k < (_SingleSolutionInfoVec[nr]->SubsetIndex).size(); ++k)
			{
				if(_SingleSolutionInfoVec[nr]->SubsetIndex[k] == SubsetIndex[i])
				{
					_SingleSolutionInfoVec[nr]->group_comp[k] = (num_comp[i])++;
				}
			}
		}
	}
	for(uint i = 0; i < SubsetIndex.size(); ++i)
	{
		if(info->num_comp[i] != num_comp[i])
		{
			std::cout << "ERROR in group_solutions: Something wrong with num_comp. Aborting" << std::endl;
			return false;
		}
	}

	_GroupSolutionInfoVec.push_back(info);

	for(uint i = 0; i < GroupSolutions.size(); ++i)
	{
		uint nr = GroupSolutions[i];
		delete _GroupSolutionInfoVec[nr];
	}

	for(uint i = GroupSolutions.size()-1; i >= 0; --i)
	{
		uint nr = GroupSolutions[i];
		_GroupSolutionInfoVec.erase(_GroupSolutionInfoVec.begin() + nr);
	}

	return true;
}

bool DoFPattern::clear_grouping()
{
	if(_lock == true)
	{
		std::cout << "Pattern already locked. Can not clear grouping" << std::endl;
		return false;
	}

	// Clear old grouping
	for(int i = 0; i< _GroupSolutionInfoVec.size();++i)
	{
		delete _GroupSolutionInfoVec[i];
	}
	_GroupSolutionInfoVec.clear();


	// restore single solution grouping
	for(uint i = 0; i < _SingleSolutionInfoVec.size(); ++i)
	{
		GroupSolutionInfo* groupInfo = new GroupSolutionInfo();
		groupInfo->name = _SingleSolutionInfoVec[i]->name;
		(groupInfo->SingleSolutions).resize(1);
		for(uint j = 0; j < 1; ++j)
		{
			groupInfo->SingleSolutions[j] = j;
		}
		(groupInfo->SubsetIndex).resize((_SingleSolutionInfoVec[i]->SubsetIndex).size());
		(groupInfo->num_comp).resize((_SingleSolutionInfoVec[i]->SubsetIndex).size());
		for(int j = 0; j < (_SingleSolutionInfoVec[i]->SubsetIndex).size(); ++j)
		{
			groupInfo->SubsetIndex[j] = _SingleSolutionInfoVec[i]->SubsetIndex[j];
			groupInfo->num_comp[j] = 1;
		}
		_GroupSolutionInfoVec.push_back(groupInfo);
	}

	return true;
}

bool DoFPattern::print_info()
{
	std::cout << "Pattern contains the following single Solutions:" << std::endl;
	for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
	{
		std::cout << i << ": '"<< _SingleSolutionInfoVec[i]->name << "' on Subset(s) ";
		for(uint j = 0; j < (_SingleSolutionInfoVec[i]->SubsetIndex).size(); ++j)
		{
			std::cout << _SingleSolutionInfoVec[i]->SubsetIndex[j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "Grouped in the following components:" << std::endl;
	for(int i = 0; i < _GroupSolutionInfoVec.size(); ++i)
	{
		std::cout << i << ": '"<< _GroupSolutionInfoVec[i]->name << "' grouping single Solution(s) ";
		for(int j = 0; j < (_GroupSolutionInfoVec[i]->SingleSolutions).size(); ++j)
		{
			std::cout << _GroupSolutionInfoVec[i]->SingleSolutions[j] << " ";
		}
		//std::cout << std::endl;
		std::cout << "with DoF (Subset id, Nr. of Comps): ";
		for(int j = 0; j < (_GroupSolutionInfoVec[i]->SubsetIndex).size(); ++j)
		{
			std::cout << "(" << _GroupSolutionInfoVec[i]->SubsetIndex[j] << ",";
			std::cout << _GroupSolutionInfoVec[i]->num_comp[j] << ") ";
		}
		std::cout << std::endl;
	}

	return true;
}

bool DoFPattern::finalize()
{
	// Check if only one Grouping present
	if(_GroupSolutionInfoVec.size() > 1)
	{
		std::cout << "ERROR in finalize: Currently only one Group Solution (with node values) allowed. Can not finalize." << std::endl;
		return false;
	}

	_lock = true;
	Grid* grid = _sh->get_assigned_grid();
	grid->register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER);

	grid->attach_to_vertices(_aDoF);
	_aaDoFVRT.access(*grid, _aDoF);

	// assign DoFs
	VertexIterator iterBegin, iterEnd, iter;

	for(uint i=0; i<_GroupSolutionInfoVec.size(); ++i)
	{
		for(uint j = 0; j < _num_sh; ++j)
		{
			uint ncomp = _GroupSolutionInfoVec[i]->num_comp[j];
			std::cout << "Assigning "<< ncomp << " Dof(s) per Vertex in Subset " << j << std::endl;
		}

		MultiGrid* mg = dynamic_cast<MultiGrid*>(grid);
		if(mg != NULL)
		{
			_num_levels = mg->num_levels();
			std::cout << "Assigning DoFs to Multigrid with " <<_num_levels << " Level(s)." << std::endl;

			_num_dofs.resize(_num_levels);
			for(uint l = 0; l < _num_levels; ++l)
			{
				iterBegin = mg->begin<Vertex>(l);
				iterEnd = mg->end<Vertex>(l);

				index_type n = 0;
				for(iter = iterBegin; iter != iterEnd; iter++)
				{
					Vertex* vrt = *iter;
					int SubsetIndex = _sh->get_subset_index(vrt);
					comp_type ncomp = _GroupSolutionInfoVec[i]->num_comp[SubsetIndex];
					_aaDoFVRT[vrt].index = n;
					_aaDoFVRT[vrt].ncomp = ncomp;
					n += ncomp;
				}
				std::cout << n << " DoFs assigned to level " << l << std::endl;
				_num_dofs[l] = n;
			}
		}
		else
		{
			std::cout << "Assigning DoFs to Grid." << std::endl;

			_num_levels = 1;
			(_num_dofs).resize(1);

			iterBegin = grid->begin<Vertex>();
			iterEnd = grid->end<Vertex>();

			index_type n = 0;
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				Vertex* vrt = *iter;
				int SubsetIndex = _sh->get_subset_index(vrt);
				comp_type ncomp = _GroupSolutionInfoVec[i]->num_comp[SubsetIndex];
				_aaDoFVRT[vrt].index = n;
				_aaDoFVRT[vrt].ncomp = ncomp;
				n += ncomp;
			}
			std::cout << n << " DoFs assigned to level " << 0 << std::endl;
			_num_dofs[0] = n;
		}
	}

	return true;
}

uint DoFPattern::num_levels()
{
	return _num_levels;
}

uint DoFPattern::num_dofs(uint level)
{
	return _num_dofs[level];
}

/*uint DoFPattern::num_dofs(int s, uint level)
{
	// TODO: include subset index in evaluation
	return _num_dofs[level];
}*/

TrialSpaceType DoFPattern::get_TrialSpaceType(uint nr_solution)
{
	return _SingleSolutionInfoVec[nr_solution]->TrialSpace;
}


DoFPattern::index_type DoFPattern::get_index(VertexBase* vrt, uint nr_solution)
{
	assert(nr_solution < _SingleSolutionInfoVec.size());
	int subsetIndex = _sh->get_subset_index(vrt);
	uint i;
	for(i = 0; i < (_SingleSolutionInfoVec[nr_solution]->SubsetIndex).size(); ++i)
	{
		if(subsetIndex == _SingleSolutionInfoVec[nr_solution]->SubsetIndex[i]) break;
	}
	if(i ==  (_SingleSolutionInfoVec[nr_solution]->SubsetIndex).size())
	{
		std::cout << "Pattern does not recognize subsetIndex " << subsetIndex << ". Aborting. \n";
		assert(0 && "ERROR in get_index. Solution not defined for SubsetIndex.");
	}

	comp_type ncomp = _SingleSolutionInfoVec[nr_solution]->group_comp[i];

	return _aaDoFVRT[vrt].index + ncomp;
}

DoFPattern::~DoFPattern()
{
	for(uint i = 0; i < _SingleSolutionInfoVec.size(); ++i)
	{
		delete _SingleSolutionInfoVec[i];
	}

	for(uint i = 0; i < _GroupSolutionInfoVec.size(); ++i)
	{
		delete _GroupSolutionInfoVec[i];
	}

	Grid* grid = _sh->get_assigned_grid();
	grid->detach_from_vertices(_aDoF);
	_aaDoFVRT.invalidate();
	grid->unregister_observer(this);

}

std::string DoFPattern::get_name(uint nr_fct)
{
	assert(nr_fct < _SingleSolutionInfoVec.size());
	return _SingleSolutionInfoVec[nr_fct]->name;
}

std::string DoFPattern::get_name()
{
	return _name;
}

int DoFPattern::num_comp()
{
	return _SingleSolutionInfoVec.size();
}

bool DoFPattern::comp_def_in_subset(uint nr_fct, int subsetIndex)
{
	for(uint i = 0; i < (_SingleSolutionInfoVec[nr_fct]->SubsetIndex).size(); ++i)
	{
		if(subsetIndex == _SingleSolutionInfoVec[nr_fct]->SubsetIndex[i])
			return true;
	}
	return false;
}

ISubsetHandler* DoFPattern::get_assigned_subset()
{
	return _sh;
}




}
