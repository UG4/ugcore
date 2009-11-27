/*
 * rhs.cpp
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#include "rhs.h"

namespace ug {


RHS::RHS(std::string name)
{
	m_name = name;
}

void RHS::set_name(std::string name)
{
	m_name = name;
}

std::string RHS::name()
{
	return m_name;
}


}
