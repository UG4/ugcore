/*
 * dofpattern.cpp
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#include "p1conform_dof_manager.h"

namespace ug{

/////////////////////////////////////
/// P1ConformDoFManager
/////////////////////////////////////


P1ConformDoFManager::
P1ConformDoFManager(grid_type& grid) : m_grid(grid)
{
	m_num_levels = grid.num_levels();
	m_bLocked = false;
}

bool
P1ConformDoFManager::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
{
	if(id != LSFS_LAGRANGEP1)
	{
		std::cout << "P1ConformDoFManager: Only LSFS_LAGRANGEP1 functions are supported. Can not add function." << std::endl;
		return false;
	}

	m_SingleDiscreteFunctionNames.push_back(name);
	return true;
}


bool
P1ConformDoFManager::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices, int dim)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager currently only works on a grid. No Subset support." << std::endl;
	return false;
}

bool
P1ConformDoFManager::
group_discrete_functions(std::vector<uint>& selected_solutions)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager groups solutions automatically.\n";

	return true;
}

bool
P1ConformDoFManager::
reset_grouping()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager groups solutions automatically.\n";

	return true;
}

bool
P1ConformDoFManager::
print_info()
{
	using namespace std;

	cout << "\n================\n";
	cout << "P1 Conform DoF-Pattern INFO: " << endl;
	cout << "================\n";

	// print discrete functions
	cout << "Number of discrete functions: " << m_SingleDiscreteFunctionNames.size() << endl;
	for(uint nr_fct = 0; nr_fct < m_SingleDiscreteFunctionNames.size(); ++nr_fct)
	{
		cout << "# " << nr_fct << ":  '" << m_SingleDiscreteFunctionNames[nr_fct] << endl;
	}

	cout << "================\n";
	return true;
}

bool
P1ConformDoFManager::
finalize()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is already locked.\n";
		return true;
	}

	// Attach indices
	m_grid.attach_to<VertexBase>(m_aIndex);
	m_aaIndex.access(m_grid, m_aIndex);

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// get number of single grid functions
	m_num_single_discrete_functions = m_SingleDiscreteFunctionNames.size();

	// loop levels
	m_num_dof_index.resize(m_grid.num_levels());
	for(uint level = 0; level < m_grid.num_levels(); ++level)
	{
		iterBegin =  m_grid.begin<VertexBase>(level);
		iterEnd = m_grid.end<VertexBase>(level);

		single_index_type i = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			VertexBase* vrt = *iter;
			m_aaIndex[vrt] = i;
			i += m_num_single_discrete_functions;
		}
		m_num_dof_index[level] = i;

		std::cout << i << " DoF indices distributed on grid " <<
					" on level " << level << " [each carrying " << 1 << " component(s)]" << std::endl;
	}

	///////////////////
	// update IndexInfo

	// loop all levels
	m_IndexInfo.resize(m_num_levels);
	for(uint level = 0; level < m_num_levels; ++level)
	{
		m_IndexInfo[level].set_leaf(true);
		m_IndexInfo[level].set_num_index(m_num_dof_index[level]);
		m_IndexInfo[level].set_num_comp(1);
	}
	return true;
}

uint
P1ConformDoFManager::
num_dofs(uint level)
{
	return m_num_dof_index[level];
}

uint
P1ConformDoFManager::
num_levels()
{
	return m_num_levels;
}

LocalShapeFunctionSetID
P1ConformDoFManager::
get_local_shape_function_set_id(uint nr_fct)
{
	return LSFS_LAGRANGEP1;
}

std::string
P1ConformDoFManager::
get_name(uint fct)
{
	return m_SingleDiscreteFunctionNames[fct];
}

uint
P1ConformDoFManager::
num_fct()
{
	return m_num_single_discrete_functions;
}

bool
P1ConformDoFManager::
fct_def_in_subset(uint nr_fct, int s)
{
	return false;
}

P1ConformDoFManager::element_container_type&
P1ConformDoFManager::
get_assigned_element_container()
{
	return m_grid;
}

P1ConformDoFManager::
~P1ConformDoFManager()
{}

}
