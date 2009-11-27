/*
 * spacialdiscretization.cpp
 *
 *  Created on: 22.10.2009
 *      Author: andreasvogel
 */

#include "spacialdiscretization.h"

namespace ug {

SpacialDiscretization::SpacialDiscretization()
{
	m_name = "No Name";
}
SpacialDiscretization::SpacialDiscretization(std::string name)
{
	m_name = name;
}

void SpacialDiscretization::set_name(std::string name)
{
	m_name = name;
}

std::string SpacialDiscretization::name()
{
	return m_name;
}

void SpacialDiscretization::add_SubsetDiscretization(SubsetDiscretization& psd)
{
	m_SubsetDiscretization.push_back(&psd);
}

void SpacialDiscretization::delete_SubsetDiscretization(int nr)
{
	m_SubsetDiscretization.erase(m_SubsetDiscretization.begin() + nr);
}

void SpacialDiscretization::clear_SubsetDiscretizations()
{
	m_SubsetDiscretization.clear();
}

int SpacialDiscretization::numberOfSubsetDiscretizations()
{
	return m_SubsetDiscretization.size();
}

void SpacialDiscretization::print_info()
{
	std::cout << "SpacialDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
	std::cout << this->numberOfSubsetDiscretizations() <<" SubsetDiscretization(s)" << std::endl;
}

}
