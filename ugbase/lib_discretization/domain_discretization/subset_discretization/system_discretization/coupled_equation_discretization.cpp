/*
 * systemdiscretization.cpp
 *
 *  Created on: 16.11.2009
 *      Author: andreasvogel
 */

#include "systemdiscretization.h"

namespace ug{

template <int d>
EquationSystem<d>::EquationSystem()
{
	m_name = "No Name";
}
template <int d>
EquationSystem<d>::EquationSystem(std::string name)
{
	m_name = name;
}

template <int d>
void EquationSystem<d>::set_name(std::string name)
{
	m_name = name;
}

template <int d>
std::string EquationSystem<d>::name()
{
	return m_name;
}


template <int d>
void EquationSystem<d>::add_Equation(Equation<d>& eq)
{
	m_Equation.push_back(&eq);
}

template <int d>
void EquationSystem<d>::delete_Equation(int nr)
{
	m_Equation.erase(m_Equation.begin() + nr);
}

template <int d>
void EquationSystem<d>::clear_Equations()
{
	m_Equation.clear();
}

template <int d>
int EquationSystem<d>::numberOfEquations()
{
	return m_Equation.size();
}

template <int d>
void EquationSystem<d>::print_info()
{
	std::cout << "SubsetDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
	std::cout << this->numberOfEquations() <<" Equation(s)" << std::endl;
}


// force code creation for dim d= 1,2,3
template class EquationSystem<1>;
template class EquationSystem<2>;
template class EquationSystem<3>;


}
