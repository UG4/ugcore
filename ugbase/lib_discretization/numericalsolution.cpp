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

NumericalSolutionPattern::NumericalSolutionPattern(std::string name, SubsetHandler& sh)
{
	_name = name;
	_sh = &sh;
	_num_sh = _sh->num_subsets();
	_lock = false;
};

bool NumericalSolutionPattern::add_solution(std::string name, TrialSpaceType Ansatzspace, int* SubsetIndex, int n)
{
	if(_lock == true)
	{
		std::cout << "Pattern already locked. Can not add solution" << std::endl;
		return false;
	}

	// Create SingleSolutionInfo
	SingleSolutionInfo* info = new SingleSolutionInfo();
	info->name = name;
	info->Ansatzspace = Ansatzspace;
	info->num_SubsetIndex = n;
	info->SubsetIndex = new int[n];
	info->group_Comp = new int[n];
	for(int i = 0; i < n; ++i)
	{
		info->SubsetIndex[i] = SubsetIndex[i];
		info->group_Comp[i] = 0;
	}
	_SingleSolutionInfoVec.push_back(info);

	// Create Grouping with only on component: the SingleSolution
	GroupSolutionInfo* groupInfo = new GroupSolutionInfo();
	groupInfo->name = name;
	groupInfo->num_SingleSolutions = 1;
	groupInfo->SingleSolutions = new int[1];
	for(int i = 0; i < 1; ++i)
	{
		groupInfo->SingleSolutions[i] = _SingleSolutionInfoVec.size() - 1;
	}
	groupInfo->num_SubsetIndex = n;
	groupInfo->SubsetIndex = new int[n];
	groupInfo->num_Comp = new int[n];
	for(int i = 0; i < n; ++i)
	{
		groupInfo->SubsetIndex[i] = SubsetIndex[i];
		groupInfo->num_Comp[i] = 1;
	}
	_GroupSolutionInfoVec.push_back(groupInfo);

	return true;
}

bool NumericalSolutionPattern::group_solution(std::string name, int* GroupSolutions, int n)
{
	if(_lock == true)
	{
		std::cout << "Pattern already locked. Can not add solution" << std::endl;
		return false;
	}

	// Sort GroupSolutions to be grouped
	std::sort(GroupSolutions, GroupSolutions + n);

	// collect subsets and num_Comp of all solutions that are to be grouped
	std::vector<int> SubsetIndex;
	std::vector<int> num_Comp;
	std::vector<int> SingleSolutions;
	for(int i = 0; i < n; ++i)
	{
		if(GroupSolutions[i] >= _GroupSolutionInfoVec.size() || GroupSolutions[i] < 0)
		{
			std::cout << "ERROR in group_solution: Trying to read Group " << GroupSolutions[i] << ", but grouping does not exist. Aborting" << std::endl;
			return false;
		}

		int nr = GroupSolutions[i];

		for(int j = 0; j < _GroupSolutionInfoVec[nr]->num_SubsetIndex; ++j)
		{
			int k;
			for(k = 0; k < SubsetIndex.size(); ++k)
			{
				if(SubsetIndex[k] == _GroupSolutionInfoVec[nr]->SubsetIndex[j])
				{
					num_Comp[k] += _GroupSolutionInfoVec[nr]->num_Comp[j];
					break;
				}
			}
			if(k == SubsetIndex.size())
			{
				SubsetIndex.push_back(_GroupSolutionInfoVec[nr]->SubsetIndex[j]);
				num_Comp.push_back(_GroupSolutionInfoVec[nr]->num_Comp[j]);
			}
		}
		for(int j = 0; j < _GroupSolutionInfoVec[nr]->num_SingleSolutions; ++j)
		{
			for(int k = 0; k < SingleSolutions.size(); ++k)
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

	std::sort(SingleSolutions.begin(), SingleSolutions.end());

	// create new grouping
	GroupSolutionInfo* info = new GroupSolutionInfo();
	info->name = name;
	info->num_SubsetIndex = SubsetIndex.size();
	info->SubsetIndex = new int[SubsetIndex.size()];
	info->num_Comp = new int[SubsetIndex.size()];
	for(int i = 0; i < SubsetIndex.size(); ++i)
	{
		info->SubsetIndex[i] = SubsetIndex[i];
		info->num_Comp[i] = num_Comp[i];
	}

	info->num_SingleSolutions = SingleSolutions.size();
	info->SingleSolutions = new int[SingleSolutions.size()];
	for(int i = 0; i < SingleSolutions.size(); ++i)
	{
		int nr = SingleSolutions[i];
		info->SingleSolutions[i] = nr;
	}

	for(int i = 0; i < SubsetIndex.size(); ++i)
	{
		num_Comp[i] = 0;
		for(int j = 0; j < SingleSolutions.size(); ++j)
		{
			int nr = SingleSolutions[j];
			int k;
			for(int k = 0; k < _SingleSolutionInfoVec[nr]->num_SubsetIndex; ++k)
			{
				//_SingleSolutionInfoVec[nr]->group_Comp[k] = -100;
				if(_SingleSolutionInfoVec[nr]->SubsetIndex[k] == SubsetIndex[i])
				{
					_SingleSolutionInfoVec[nr]->group_Comp[k] = (num_Comp[i])++;
				}
			}
		}
	}
	for(int i = 0; i < SubsetIndex.size(); ++i)
	{
		if(info->num_Comp[i] != num_Comp[i])
		{
			std::cout << "ERROR in group_solutions: Something wrong with num_Comp. Aborting" << std::endl;
			return false;
		}
	}

	_GroupSolutionInfoVec.push_back(info);

	for(int i = 0; i < n; ++i)
	{
		int nr = GroupSolutions[i];
		delete _GroupSolutionInfoVec[nr]->SubsetIndex;
		delete _GroupSolutionInfoVec[nr]->num_Comp;
		delete _GroupSolutionInfoVec[nr]->SingleSolutions;
		delete _GroupSolutionInfoVec[nr];
	}

	for(int i = n-1; i >= 0; --i)
	{
		int nr = GroupSolutions[i];
		_GroupSolutionInfoVec.erase(_GroupSolutionInfoVec.begin() + nr);
	}

	return true;
}

bool NumericalSolutionPattern::clear_grouping()
{
	if(_lock == true)
	{
		std::cout << "Pattern already locked. Can not clear grouping" << std::endl;
		return false;
	}

	// Clear old grouping
	for(int i = 0; i< _GroupSolutionInfoVec.size();++i)
	{
		delete _GroupSolutionInfoVec[i]->SubsetIndex;
		delete _GroupSolutionInfoVec[i]->num_Comp;
		delete _GroupSolutionInfoVec[i]->SingleSolutions;
		delete _GroupSolutionInfoVec[i];
	}
	_GroupSolutionInfoVec.clear();


	// restore single solution grouping
	for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
	{
		GroupSolutionInfo* groupInfo = new GroupSolutionInfo();
		groupInfo->name = _SingleSolutionInfoVec[i]->name;
		groupInfo->num_SingleSolutions = 1;
		groupInfo->SingleSolutions = new int[1];
		for(int j = 0; j < 1; ++j)
		{
			groupInfo->SingleSolutions[j] = j;
		}
		groupInfo->num_SubsetIndex = _SingleSolutionInfoVec[i]->num_SubsetIndex;
		groupInfo->SubsetIndex = new int[groupInfo->num_SubsetIndex];
		groupInfo->num_Comp = new int[groupInfo->num_SubsetIndex];
		for(int j = 0; j < groupInfo->num_SubsetIndex; ++j)
		{
			groupInfo->SubsetIndex[j] = _SingleSolutionInfoVec[i]->SubsetIndex[j];
			groupInfo->num_Comp[j] = 1;
		}
		_GroupSolutionInfoVec.push_back(groupInfo);
	}

	return true;
}

bool NumericalSolutionPattern::print_info()
{
	std::cout << "Pattern contains the following single Solutions:" << std::endl;
	for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
	{
		std::cout << i << ": '"<< _SingleSolutionInfoVec[i]->name << "' on Subset(s) ";
		for(int j = 0; j < _SingleSolutionInfoVec[i]->num_SubsetIndex; ++j)
		{
			std::cout << _SingleSolutionInfoVec[i]->SubsetIndex[j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "Grouped in the following components:" << std::endl;
	for(int i = 0; i < _GroupSolutionInfoVec.size(); ++i)
	{
		std::cout << i << ": '"<< _GroupSolutionInfoVec[i]->name << "' grouping single Solution(s) ";
		for(int j = 0; j < _GroupSolutionInfoVec[i]->num_SingleSolutions; ++j)
		{
			std::cout << _GroupSolutionInfoVec[i]->SingleSolutions[j] << " ";
		}
		//std::cout << std::endl;
		std::cout << "with DoF (Subset id, Nr. of Comps): ";
		for(int j = 0; j < _GroupSolutionInfoVec[i]->num_SubsetIndex; ++j)
		{
			std::cout << "(" << _GroupSolutionInfoVec[i]->SubsetIndex[j] << ",";
			std::cout << _GroupSolutionInfoVec[i]->num_Comp[j] << ") ";
		}
		std::cout << std::endl;
	}

	return true;
}

bool NumericalSolutionPattern::finalize()
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

	for(int i=0; i<_GroupSolutionInfoVec.size(); ++i)
	{
		for(int j = 0; j < _sh->num_subsets(); ++j)
		{
			int ncomp = _GroupSolutionInfoVec[i]->num_Comp[j];
			std::cout << "Assigning "<< ncomp << " Dof(s) per Vertex in Subset " << j << std::endl;
		}

		MultiGrid* mg = dynamic_cast<MultiGrid*>(grid);
		if(mg != NULL)
		{
			_num_levels = mg->num_levels();
			std::cout << "Assigning DoFs to Multigrid with " <<_num_levels << " Level(s)." << std::endl;

			_num_doubles = new int[_num_levels];
			for(int l = 0; l < _num_levels; ++l)
			{
				iterBegin = mg->begin<Vertex>(l);
				iterEnd = mg->end<Vertex>(l);

				int n = 0;
				for(iter = iterBegin; iter != iterEnd; iter++)
				{
					Vertex* vrt = *iter;
					int SubsetIndex = _sh->get_subset_index(vrt);
					int ncomp = _GroupSolutionInfoVec[i]->num_Comp[SubsetIndex];
					_aaDoFVRT[vrt].nr = n;
					_aaDoFVRT[vrt].ncomp = ncomp;
					n += ncomp;
				}
				std::cout << n << " DoFs assigned to level" << l << std::endl;
				_num_doubles[l] = n;
			}
		}
		else
		{
			std::cout << "Assigning DoFs to Grid." << std::endl;

			_num_levels = 1;
			_num_doubles = new int[1];

			iterBegin = grid->begin<Vertex>();
			iterEnd = grid->end<Vertex>();

			int n = 0;
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				Vertex* vrt = *iter;
				int SubsetIndex = _sh->get_subset_index(vrt);
				int ncomp = _GroupSolutionInfoVec[i]->num_Comp[SubsetIndex];
				_aaDoFVRT[vrt].nr = n;
				_aaDoFVRT[vrt].ncomp = ncomp;
				n += ncomp;
			}
			std::cout << n << " DoFs assigned to level" << 0 << std::endl;
			_num_doubles[0] = n;
		}
	}

	return true;
}

int NumericalSolutionPattern::num_levels()
{
	return _num_levels;
}

int NumericalSolutionPattern::num_doubles(int i)
{
	assert(i >= 0 && "ERROR in num_doubles: negative level index");
	return _num_doubles[i];
}

int NumericalSolutionPattern::num_doubles()
{
	return _num_doubles[0];
}

TrialSpaceType NumericalSolutionPattern::get_trial_space_type(int nr_solution)
{
	return _SingleSolutionInfoVec[nr_solution]->Ansatzspace;
}


int NumericalSolutionPattern::get_index(VertexBase* vrt, int nr_solution)
{
	assert(nr_solution < _SingleSolutionInfoVec.size());
	int subsetIndex = _sh->get_subset_index(vrt);
	int i;
	for(i = 0; i < _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex; ++i)
	{
		if(subsetIndex == _SingleSolutionInfoVec[nr_solution]->SubsetIndex[i]) break;
	}
	if(i ==  _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex)
	{
		std::cout << "Pattern does not recognize subsetIndex " << subsetIndex << ". Aborting. \n";
		assert(0 && "ERROR in get_index. Solution not defined for SubsetIndex.");
	}

	int ncomp = _SingleSolutionInfoVec[nr_solution]->group_Comp[i];

	return _aaDoFVRT[vrt].nr + ncomp;
}

NumericalSolutionPattern::~NumericalSolutionPattern()
{
	for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
	{
		delete _SingleSolutionInfoVec[i]->SubsetIndex;
		delete _SingleSolutionInfoVec[i]->group_Comp;
		delete _SingleSolutionInfoVec[i];
	}

	for(int i = 0; i < _GroupSolutionInfoVec.size(); ++i)
	{
		delete _GroupSolutionInfoVec[i]->SubsetIndex;
		delete _GroupSolutionInfoVec[i]->num_Comp;
		delete _GroupSolutionInfoVec[i]->SingleSolutions;
		delete _GroupSolutionInfoVec[i];
	}

	Grid* grid = _sh->get_assigned_grid();
	grid->detach_from_vertices(_aDoF);
	_aaDoFVRT.invalidate();
	grid->unregister_observer(this);

}

std::string NumericalSolutionPattern::get_name(int comp)
{
	assert(comp < _SingleSolutionInfoVec.size());
	return _SingleSolutionInfoVec[comp]->name;
}

std::string NumericalSolutionPattern::get_name()
{
	return _name;
}

int NumericalSolutionPattern::num_comp()
{
	return _SingleSolutionInfoVec.size();
}

bool NumericalSolutionPattern::comp_def_in_subset(int nrComp, int subsetIndex)
{
	for(int i = 0; i < _SingleSolutionInfoVec[nrComp]->num_SubsetIndex; ++i)
	{
		if(subsetIndex == _SingleSolutionInfoVec[nrComp]->SubsetIndex[i])
			return true;
	}
	return false;
}

SubsetHandler* NumericalSolutionPattern::get_assigned_subset()
{
	return _sh;
}



////////////////////////////
// numerical solution
////////////////////////////




NumericalSolution::NumericalSolution(std::string shortname, std::string longname, std::string description, NumericalSolutionPattern& pattern, Grid& grid)
{
	m_shortname = shortname;
	m_longname = longname;
	m_description = description;
	m_grid = &grid;

	_pattern = &pattern;

	_GridVector.resize(_pattern->num_levels());
	for(int i = 0; i < _pattern->num_levels(); ++i)
	{
		_GridVector[i] = new Vector();
		_GridVector[i]->create_vector(_pattern->num_doubles(i));
	}

	return;
}

bool NumericalSolution::set_values(ValueFunction fct, int nr_func, SubsetHandler& sh, uint subsetIndex)
{
	return set_values(fct, nr_func, sh, subsetIndex, 0);
}

bool NumericalSolution::set_values(ValueFunction fct, int nr_func, SubsetHandler& sh, uint subsetIndex, int level)
{
	VertexBaseIterator iterBegin, iterEnd, iter;

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<APosition> aaPos(*grid, aPosition);

	MultiGrid* mg = dynamic_cast<MultiGrid*>(grid);

	if(mg == NULL)
	{
		if(level != 0)
		{
			std::cout << "Numerical Solution is only defined on a Grid, no MultiGrid assignment possible." << std::endl;
			return false;
		}
		iterBegin = grid->begin<VertexBase>();
		iterEnd = grid->end<VertexBase>();
	}
	else
	{
		iterBegin = mg->begin<VertexBase>(level);
		iterEnd = mg->end<VertexBase>(level);
	}

	int index;
	number val;
	MathVector<3> corner;
	std::vector<int> indices;
	std::vector<number> values;

	/* loop over all Vertices */
	int nvertices = 0;
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		if(sh.get_subset_index(vert) != subsetIndex) continue;

		corner = aaPos[vert];
		index = _pattern->get_index(vert, nr_func);
		fct(corner, val);
		values.push_back(val);
		indices.push_back(index);

		nvertices++;
	}

	double* valueArray = new double[indices.size()];
	int* indexArray = new int[indices.size()];

	for(int i=0; i< indices.size(); i++)
	{
		valueArray[i] = values[i];
		indexArray[i] = indices[i];
	}

	if(_GridVector[level]->set_values((int) indices.size(), indexArray, valueArray) != true)
		return false;
	delete valueArray;
	delete indexArray;

	return true;
}

Vector* NumericalSolution::GridVector(int level)
{
	return _GridVector[level];
}

Vector* NumericalSolution::GridVector()
{
	return _GridVector[0];
}

NumericalSolutionPattern* NumericalSolution::get_pattern()
{
	return _pattern;
}


bool NumericalSolution::print()
{
	std::cout << "shortname: \"" << m_shortname << "\", ";
	std::cout << "longname: \"" << m_longname <<"\", ";
	std::cout << "description: \"" << m_description << "\"" << std::endl;
	return true;
}

NumericalSolution::~NumericalSolution()
{
	for(int i = 0; i < _pattern->num_levels(); ++i)
	{
		delete _GridVector[i];
	}

}

std::string NumericalSolution::get_name(int comp)
{
	return _pattern->get_name(comp);
}



}
