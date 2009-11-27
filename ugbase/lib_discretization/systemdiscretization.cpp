/*
 * systemdiscretization.cpp
 *
 *  Created on: 16.11.2009
 *      Author: andreasvogel
 */

#include "systemdiscretization.h"

namespace ug{

EquationSystem::EquationSystem()
{
	m_name = "No Name";
}
EquationSystem::EquationSystem(std::string name)
{
	m_name = name;
}

void EquationSystem::set_name(std::string name)
{
	m_name = name;
}

std::string EquationSystem::name()
{
	return m_name;
}


void EquationSystem::add_Equation(Equation& eq)
{
	m_Equation.push_back(&eq);
}

void EquationSystem::delete_Equation(int nr)
{
	m_Equation.erase(m_Equation.begin() + nr);
}

void EquationSystem::clear_Equations()
{
	m_Equation.clear();
}

int EquationSystem::numberOfEquations()
{
	return m_Equation.size();
}

void EquationSystem::print_info()
{
	std::cout << "SubsetDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
	std::cout << this->numberOfEquations() <<" Equation(s)" << std::endl;
}

}
