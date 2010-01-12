/*
 * rhs.cpp
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

namespace ug {

template <int d>
RHS<d>::RHS(std::string name)
{
	m_name = name;
}

template <int d>
void RHS<d>::set_name(std::string name)
{
	m_name = name;
}

template <int d>
std::string RHS<d>::name()
{
	return m_name;
}


}
