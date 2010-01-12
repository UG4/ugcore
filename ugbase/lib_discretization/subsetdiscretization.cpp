/*
 * problem.cpp
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#include "subsetdiscretization.h"

namespace ug {

template <int d>
SubsetDiscretization<d>::SubsetDiscretization()
{
	m_name = "No Name";
}
template <int d>
SubsetDiscretization<d>::SubsetDiscretization(std::string name)
{
	m_name = name;
}

template <int d>
void SubsetDiscretization<d>::set_name(std::string name)
{
	m_name = name;
}

template <int d>
std::string SubsetDiscretization<d>::name()
{
	return m_name;
}

template <int d>
void SubsetDiscretization<d>::add_Subset(SubsetHandler& sh, uint subsetIndex)
{
	m_sh = &sh;
	m_subsetIndex = subsetIndex;
}

template <int d>
void SubsetDiscretization<d>::clear_Subset()
{
	m_sh = NULL;
	m_subsetIndex = 0;
}

template <int d>
void SubsetDiscretization<d>::add_System(SystemDiscretization<d>& sys)
{
	m_SystemDiscretization.push_back(&sys);
}

template <int d>
void SubsetDiscretization<d>::delete_System(int nr)
{
	m_SystemDiscretization.erase(m_SystemDiscretization.begin() + nr);
}

template <int d>
void SubsetDiscretization<d>::clear_Systems()
{
	m_SystemDiscretization.clear();
}

template <int d>
int SubsetDiscretization<d>::numberOfSystems()
{
	return m_SystemDiscretization.size();
}

template <int d>
void SubsetDiscretization<d>::print_info()
{
	std::cout << "SubsetDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
}


// force code creation for dim d=1,2,3
	template class SubsetDiscretization<1>;
	template class SubsetDiscretization<2>;
	template class SubsetDiscretization<3>;

}
