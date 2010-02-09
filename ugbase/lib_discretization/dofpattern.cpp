/*
 * dofpattern.cpp
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#include "dofpattern.h"

namespace ug{

/////////////////////////////////////
/// DoFPattern
/////////////////////////////////////


DoFPattern::
DoFPattern(std::string name, subset_handler_type& sh) : m_sh(sh), m_name(name)
{
	m_num_subsets = sh.num_subsets();
	m_num_levels = sh.num_levels();
	m_bLocked = false;
}

bool
DoFPattern::
add_discrete_function(std::string name, TrialSpaceType TrialSpace, std::vector<int>& SubsetIndices, uint dim)
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
	m_SingleDiscreteFunctionVec[nr_fct].set_trial_space(TrialSpace);
	m_SingleDiscreteFunctionVec[nr_fct].set_trial_space_dim(dim);
	m_SingleDiscreteFunctionVec[nr_fct].set_subset_indices(SubsetIndices, m_num_subsets);

	// get Element DoF Pattern
	const ElementDoFPattern& elemDoFPattern = TrialSpaces::get_element_dof_pattern(TrialSpace, dim);

	// Now we have to add a new DoFGroup for this Discrete Function.
	for(uint i = 0; i < SubsetIndices.size(); ++i)
	{
		const int subsetIndex = SubsetIndices[i];

		uint num_dofs = elemDoFPattern.num_dofs<Vertex>();
		add_dof_group<VertexBase>(nr_fct, subsetIndex, num_dofs);

		// same for other objects
		// [ ... ]
	}
	return true;
}

bool
DoFPattern::
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
DoFPattern::
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
DoFPattern::
print_info()
{
	using namespace std;

	cout << "\n================\n";
	cout << "DoF-Pattern INFO: " << m_name << endl;
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

bool
DoFPattern::
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

	return true;
}

uint
DoFPattern::
num_dofs(uint level)
{
	return 0;
}

uint
DoFPattern::
num_levels()
{
	return m_num_levels;
}

TrialSpaceType
DoFPattern::
get_TrialSpaceType(uint nr_fct)
{
	return m_SingleDiscreteFunctionVec[nr_fct].get_trial_space();
}

std::string
DoFPattern::
get_name(uint nr_fct)
{
	return m_SingleDiscreteFunctionVec[nr_fct].get_name();
}

std::string
DoFPattern::
get_name()
{
	return m_name;
}

uint
DoFPattern::
num_fct()
{
	return m_SingleDiscreteFunctionVec.size();
}

bool
DoFPattern::
fct_def_in_subset(uint nr_fct, int s)
{
	return false;
}

DoFPattern::subset_handler_type&
DoFPattern::
get_assigned_subset()
{
	return m_sh;
}

DoFPattern::
~DoFPattern()
{}



/////////////////////////
// DoF Group Map
/////////////////////////
DoFPattern::DoFGroupMap::
DoFGroupMap() : m_num_dofs(0)
{}

DoFPattern::DoFGroupMap::
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

DoFPattern::DoFGroupMap&
DoFPattern::DoFGroupMap::
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
}

bool
DoFPattern::DoFGroupMap::
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
DoFPattern::DoFGroupMap::
set_nr_dofgroup(uint dof, uint group)
{
	if(dof >= num_dofs()) return false;
	m_nr_dofgroup[dof] = group;
	return true;
}

bool
DoFPattern::DoFGroupMap::
set_comp_in_dofgroup(uint dof, uint comp)
{
	if(dof >= num_dofs()) return false;
	m_comp_in_dofgroup[dof] = comp;
	return true;
}

DoFPattern::DoFGroupMap::
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
DoFPattern::SingleDiscreteFunction::
set_name(std::string name)
{
	m_name = name;
}

std::string
DoFPattern::SingleDiscreteFunction::
get_name() const
{
	return m_name;
}

void
DoFPattern::SingleDiscreteFunction::
set_subset_handler(subset_handler_type& sh)
{
	m_sh = &sh;
}

DoFPattern::subset_handler_type*
DoFPattern::SingleDiscreteFunction::
get_subset_handler() const
{
	return m_sh;
}

void
DoFPattern::SingleDiscreteFunction::
set_element_dof_pattern(ElementDoFPattern& elemPattern)
{
	m_elementDoFPattern = elemPattern;
}

const ElementDoFPattern&
DoFPattern::SingleDiscreteFunction::
get_element_dof_pattern() const
{
	return m_elementDoFPattern;
}

void
DoFPattern::SingleDiscreteFunction::
set_trial_space(TrialSpaceType tsp)
{
	m_trialSpace = tsp;
}

TrialSpaceType
DoFPattern::SingleDiscreteFunction::
get_trial_space() const
{
	return m_trialSpace;
}

void
DoFPattern::SingleDiscreteFunction::
set_trial_space_dim(uint dim)
{
	m_trialSpaceDim = dim;
}

uint
DoFPattern::SingleDiscreteFunction::
get_trial_space_dim() const
{
	return m_trialSpaceDim;
}

void
DoFPattern::SingleDiscreteFunction::
set_subset_indices(std::vector<int> subsetIndices, uint num_subsets)
{
	m_subsetIndices = subsetIndices;
	for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
	{
		m_dofGroupMap[i].resize(num_subsets);
	}
}

bool
DoFPattern::SingleDiscreteFunction::
is_def_in_subset(int s) const
{
	std::vector<int>::const_iterator iter;
	iter = find (m_subsetIndices.begin(), m_subsetIndices.end(), s);
	if(iter == m_subsetIndices.end()) return false;
	return true;
}



}
