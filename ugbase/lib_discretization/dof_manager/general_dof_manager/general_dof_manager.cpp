/*
 * dofpattern.cpp
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#include "general_dof_manager.h"

namespace ug{

/////////////////////////////////////
/// GeneralDoFManager
/////////////////////////////////////


GeneralDoFManager::
GeneralDoFManager(subset_handler_type& sh) : m_sh(sh)
{
	m_num_subsets = sh.num_subsets();
	m_num_levels = sh.num_levels();
	m_bLocked = false;
}

bool
GeneralDoFManager::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
{
	std::vector<int> SubsetIndices;

	for(int i = 0; i < (int) m_sh.num_subsets(); ++i)
	{
		SubsetIndices.push_back(i);
	}

	return add_discrete_function(name, id, SubsetIndices, dim);
}


bool
GeneralDoFManager::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices, int dim)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	// check if all subset indices are correct
	for(uint i = 0; i < SubsetIndices.size(); ++i)
	{
		if(SubsetIndices[i] >= m_num_subsets) return false;
	}

	// check, that dimension is correct
	if(dim > 3) return false;

	// remember nr of added function
	uint nr_fct = m_SingleDiscreteFunctionVec.size();

	// create new single solution
	m_SingleDiscreteFunctionVec.push_back(SingleDiscreteFunction());

	// set member functions of newly created Discrete Function
	m_SingleDiscreteFunctionVec[nr_fct].set_name(name);
	m_SingleDiscreteFunctionVec[nr_fct].set_subset_handler(m_sh);
	m_SingleDiscreteFunctionVec[nr_fct].set_local_shape_function_set_id(id);
	m_SingleDiscreteFunctionVec[nr_fct].set_dim(dim);
	m_SingleDiscreteFunctionVec[nr_fct].set_subset_indices(SubsetIndices, m_num_subsets);

	// create need ContinuousDoFManager
	ContinuousDoFPattern contDoFPattern;

	if(dim == 2 )
	{
		const LocalShapeFunctionSet<ReferenceTriangle>& p1 =
			LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ReferenceTriangle>(id);
		contDoFPattern.add_local_dof_pattern(p1.local_dof_pattern());
		const LocalShapeFunctionSet<ReferenceQuadrilateral>& p2 =
			LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ReferenceQuadrilateral>(id);
		contDoFPattern.add_local_dof_pattern(p2.local_dof_pattern());
	}
	if(dim == 3 || dim == 1)
	{
		std::cout << "Not implemented" << std::endl;
		return false;
	}

	m_SingleDiscreteFunctionVec[nr_fct].set_continuous_dof_pattern(contDoFPattern);

	// Now we have to add a new DoFGroup for this Discrete Function.
	for(uint i = 0; i < SubsetIndices.size(); ++i)
	{
		const int subsetIndex = SubsetIndices[i];

		uint num_dofs = contDoFPattern.num_dofs(ROID_VERTEX);
		add_dof_group<VertexBase>(nr_fct, subsetIndex, num_dofs);

		// same for other objects
		// [ ... ]
	}
	return true;
}

bool
GeneralDoFManager::
group_discrete_functions(std::vector<uint>& selected_solutions)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	if(group_discrete_functions<VertexBase>(selected_solutions) != true) return false;

	// same for other objects
	// [ ... ]
	return true;
}

bool
GeneralDoFManager::
reset_grouping()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	if(reset_grouping<VertexBase>() != true) return false;

	// same for other objects
	// [ ... ]
	return true;
}

bool
GeneralDoFManager::
print_info()
{
	using namespace std;

	cout << "\n================\n";
	cout << "DoF-Pattern INFO: " << endl;
	cout << "================\n";

	// print discrete functions
	cout << "Number of discrete functions: " << m_SingleDiscreteFunctionVec.size() << endl;
	for(uint nr_fct = 0; nr_fct < m_SingleDiscreteFunctionVec.size(); ++nr_fct)
	{
		cout << "# " << nr_fct << ":  '" << m_SingleDiscreteFunctionVec[nr_fct].get_name() << "' on Subset(s): ";
		for(int j = 0; j < m_num_subsets; ++j)
		{
			if(m_SingleDiscreteFunctionVec[nr_fct].is_def_in_subset(j))
				cout << j << " ";
		}
		cout << endl;
	}

	// print dof groups
	cout << endl << "Vertex DoFGroup distribution on the subset is as follows:" << endl;
	for(int subsetIndex = 0; subsetIndex < m_num_subsets; ++subsetIndex)
	{
		cout << "Subset # " << subsetIndex << ": " << endl;
		for(uint nr_group = 0; nr_group < m_VertexBaseDoFGroups.size(); ++nr_group)
		{
			const DoFGroup<VertexBase>& group = m_VertexBaseDoFGroups[nr_group];
			if(group.get_subset_index() != subsetIndex)	continue;

			cout << "  DoFGroup # " << nr_group << " is a Vertex Vector of this components: [ ";
			for(uint dof = 0; dof < group.num_dofs(); ++dof)
			{
				const int nr = group.get_nr_fct(dof);
				cout << nr << "(" << m_SingleDiscreteFunctionVec[nr].get_name() << ") ";
				if(dof < group.num_dofs() - 1) cout << "; ";
			}
			cout << "]" << endl;
		}
	}

	cout << "================\n";
	return true;
}

// returns true if group_id found
bool
GeneralDoFManager::
get_dof_group_index_info(int subsetIndex, uint group_id, uint level, uint& num_indices, uint& num_comp)
{
	if(get_dof_group_index_info<VertexBase>(subsetIndex, group_id, level, num_indices, num_comp)) return true;
	// same for other objects
	// [ ... ]

	return false;
}

bool
GeneralDoFManager::
finalize()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is already locked.\n";
		return true;
	}

	if(finalize<VertexBase>() != true) return false;

	// same for other objects
	// [ ... ]


	///////////////////
	// update IndexInfo

	// loop all levels
	m_IndexInfo.resize(m_num_levels);

	for(uint level = 0; level < m_num_levels; ++level)
	{
		m_IndexInfo[level].set_leaf(false);
		m_IndexInfo[level].set_num_index(m_num_subsets);

		// loop all subsets
		for(int subsetIndex = 0; subsetIndex < m_num_subsets; ++subsetIndex)
		{
			IndexInfo& subsetIndexInfo = m_IndexInfo[level].get_index_info(subsetIndex);

			// give ids to dof groups and count dof groups
			uint num_dof_groups = 0;
			num_dof_groups += set_dof_group_ids<VertexBase>(subsetIndex, num_dof_groups);

			// same for other objects
			// [ ... ]

			subsetIndexInfo.set_leaf(false);
			subsetIndexInfo.set_num_index(num_dof_groups);

			// set group infos
			for(uint group = 0; group < num_dof_groups; ++group)
			{
				uint num_indices;
				uint num_comp;

				if(get_dof_group_index_info(subsetIndex, group, level, num_indices, num_comp) != true) return false;

				IndexInfo& groupIndexInfo = subsetIndexInfo.get_index_info(group);
				groupIndexInfo.set_leaf(true);
				groupIndexInfo.set_num_index(num_indices);
				groupIndexInfo.set_num_comp(num_comp);
			}
		}
	}
	return true;
}

uint
GeneralDoFManager::
num_dofs(uint level)
{
	return 0;
}

uint
GeneralDoFManager::
num_levels()
{
	return m_num_levels;
}

LocalShapeFunctionSetID
GeneralDoFManager::
get_local_shape_function_set_id(uint nr_fct)
{
	return m_SingleDiscreteFunctionVec[nr_fct].get_local_shape_function_set_id();
}

std::string
GeneralDoFManager::
get_name(uint nr_fct)
{
	return m_SingleDiscreteFunctionVec[nr_fct].get_name();
}

uint
GeneralDoFManager::
num_fct()
{
	return m_SingleDiscreteFunctionVec.size();
}

bool
GeneralDoFManager::
fct_def_in_subset(uint nr_fct, int s)
{
	return false;
}

GeneralDoFManager::subset_handler_type&
GeneralDoFManager::
get_assigned_subset()
{
	return m_sh;
}

GeneralDoFManager::
~GeneralDoFManager()
{}



/////////////////////////
// DoF Group Map
/////////////////////////
GeneralDoFManager::DoFGroupMap::
DoFGroupMap() : m_num_dofs(0)
{}

GeneralDoFManager::DoFGroupMap::
DoFGroupMap(const DoFGroupMap& map)
{
	m_num_dofs = map.m_num_dofs;
	m_nr_dofgroup = new uint[m_num_dofs];
	//if(m_nr_dofgroup == 0) return false;
	m_comp_in_dofgroup = new uint[m_num_dofs];
	//if(m_comp_in_dofgroup == 0) return false;
	for(uint i = 0; i < m_num_dofs; ++i)
	{
		m_nr_dofgroup[i] = map.m_nr_dofgroup[i];
		m_comp_in_dofgroup[i] = map.m_comp_in_dofgroup[i];
	}
}

GeneralDoFManager::DoFGroupMap&
GeneralDoFManager::DoFGroupMap::
operator= (const DoFGroupMap& map)
{
	m_num_dofs = map.m_num_dofs;
	m_nr_dofgroup = new uint[m_num_dofs];
	//if(m_nr_dofgroup == 0) return false;
	m_comp_in_dofgroup = new uint[m_num_dofs];
	//if(m_comp_in_dofgroup == 0) return false;
	for(uint i = 0; i < m_num_dofs; ++i)
	{
		m_nr_dofgroup[i] = map.m_nr_dofgroup[i];
		m_comp_in_dofgroup[i] = map.m_comp_in_dofgroup[i];
	}
	return *this;
}

bool
GeneralDoFManager::DoFGroupMap::
set_num_dofs(uint num_dofs)
{
	if(m_num_dofs != 0)
	{
		delete[] m_nr_dofgroup;
		delete[] m_comp_in_dofgroup;
	}
	m_num_dofs = num_dofs;
	m_nr_dofgroup = new uint[m_num_dofs];
	if(m_nr_dofgroup == 0) return false;
	m_comp_in_dofgroup = new uint[m_num_dofs];
	if(m_comp_in_dofgroup == 0) return false;
	return true;
}

bool
GeneralDoFManager::DoFGroupMap::
set_nr_dofgroup(uint dof, uint group)
{
	if(dof >= num_dofs()) return false;
	m_nr_dofgroup[dof] = group;
	return true;
}

bool
GeneralDoFManager::DoFGroupMap::
set_comp_in_dofgroup(uint dof, uint comp)
{
	if(dof >= num_dofs()) return false;
	m_comp_in_dofgroup[dof] = comp;
	return true;
}

GeneralDoFManager::DoFGroupMap::
~DoFGroupMap()
{
	if(m_num_dofs != 0)
	{
		delete[] m_nr_dofgroup;
		delete[] m_comp_in_dofgroup;
	}
}


///////////////////////////////
// Single Discrete Function
///////////////////////////////


void
GeneralDoFManager::SingleDiscreteFunction::
set_name(std::string name)
{
	m_name = name;
}

std::string
GeneralDoFManager::SingleDiscreteFunction::
get_name() const
{
	return m_name;
}

void
GeneralDoFManager::SingleDiscreteFunction::
set_subset_handler(subset_handler_type& sh)
{
	m_sh = &sh;
}

GeneralDoFManager::subset_handler_type*
GeneralDoFManager::SingleDiscreteFunction::
get_subset_handler() const
{
	return m_sh;
}

void
GeneralDoFManager::SingleDiscreteFunction::
set_continuous_dof_pattern(ContinuousDoFPattern& p)
{
	m_continuousDoFManager = p;
}

const ContinuousDoFPattern&
GeneralDoFManager::SingleDiscreteFunction::
get_continuous_dof_pattern() const
{
	return m_continuousDoFManager;
}

void
GeneralDoFManager::SingleDiscreteFunction::
set_local_shape_function_set_id(LocalShapeFunctionSetID id)
{
	m_fct_set_id = id;
}

LocalShapeFunctionSetID
GeneralDoFManager::SingleDiscreteFunction::
get_local_shape_function_set_id() const
{
	return m_fct_set_id;
}

void
GeneralDoFManager::SingleDiscreteFunction::
set_dim(int dim)
{
	m_dim = dim;
}

int
GeneralDoFManager::SingleDiscreteFunction::
get_dim() const
{
	return m_dim;
}

void
GeneralDoFManager::SingleDiscreteFunction::
set_subset_indices(std::vector<int> subsetIndices, int num_subsets)
{
	m_subsetIndices = subsetIndices;
	for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
	{
		m_dofGroupMap[i].resize(num_subsets);
	}
}

bool
GeneralDoFManager::SingleDiscreteFunction::
is_def_in_subset(int s) const
{
	std::vector<int>::const_iterator iter;
	iter = find (m_subsetIndices.begin(), m_subsetIndices.end(), s);
	if(iter == m_subsetIndices.end()) return false;
	return true;
}



}
