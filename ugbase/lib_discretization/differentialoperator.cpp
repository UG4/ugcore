/*
 * differentialoperator.cpp
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#include "differentialoperator.h"

namespace ug {


DifferentialOperator::DifferentialOperator(std::string name)
{
	m_name = name;
}

void DifferentialOperator::set_name(std::string name)
{
	m_name = name;
}

std::string DifferentialOperator::name()
{
	return m_name;
}


}
