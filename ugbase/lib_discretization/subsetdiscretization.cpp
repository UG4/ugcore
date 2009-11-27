/*
 * problem.cpp
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#include "subsetdiscretization.h"

namespace ug {

SubsetDiscretization::SubsetDiscretization()
{
	m_name = "No Name";
}
SubsetDiscretization::SubsetDiscretization(std::string name)
{
	m_name = name;
}

void SubsetDiscretization::set_name(std::string name)
{
	m_name = name;
}

std::string SubsetDiscretization::name()
{
	return m_name;
}

void SubsetDiscretization::add_Subset(SubsetHandler& sh, uint subsetIndex)
{
	m_sh = &sh;
	m_subsetIndex = subsetIndex;
}

void SubsetDiscretization::clear_Subset()
{
	m_sh = NULL;
	m_subsetIndex = 0;
}

void SubsetDiscretization::add_System(SystemDiscretization& sys)
{
	m_SystemDiscretization.push_back(&sys);
}

void SubsetDiscretization::delete_System(int nr)
{
	m_SystemDiscretization.erase(m_SystemDiscretization.begin() + nr);
}

void SubsetDiscretization::clear_Systems()
{
	m_SystemDiscretization.clear();
}

int SubsetDiscretization::numberOfSystems()
{
	return m_SystemDiscretization.size();
}

void SubsetDiscretization::print_info()
{
	std::cout << "SubsetDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
}


}
