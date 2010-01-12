/*
 * spacialdiscretization.cpp
 *
 *  Created on: 22.10.2009
 *      Author: andreasvogel
 */

#include "spacialdiscretization.h"

namespace ug {

template <int d>
SpacialDiscretization<d>::SpacialDiscretization()
{
	m_name = "No Name";
}
template <int d>
SpacialDiscretization<d>::SpacialDiscretization(std::string name)
{
	m_name = name;
}

template <int d>
void SpacialDiscretization<d>::set_name(std::string name)
{
	m_name = name;
}

template <int d>
std::string SpacialDiscretization<d>::name()
{
	return m_name;
}

template <int d>
void SpacialDiscretization<d>::add_SubsetDiscretization(SubsetDiscretization<d>& psd)
{
	m_SubsetDiscretization.push_back(&psd);
}

template <int d>
void SpacialDiscretization<d>::delete_SubsetDiscretization(int nr)
{
	m_SubsetDiscretization.erase(m_SubsetDiscretization.begin() + nr);
}

template <int d>
void SpacialDiscretization<d>::clear_SubsetDiscretizations()
{
	m_SubsetDiscretization.clear();
}

template <int d>
int SpacialDiscretization<d>::numberOfSubsetDiscretizations()
{
	return m_SubsetDiscretization.size();
}

template <int d>
void SpacialDiscretization<d>::print_info()
{
	std::cout << "SpacialDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
	std::cout << this->numberOfSubsetDiscretizations() <<" SubsetDiscretization(s)" << std::endl;
}

// force code creation for dim d=1,2,3
template class SpacialDiscretization<1>;
template class SpacialDiscretization<2>;
template class SpacialDiscretization<3>;


}
