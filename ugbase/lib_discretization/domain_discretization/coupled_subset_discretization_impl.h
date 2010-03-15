/*
 * spacialdiscretization_impl.h
 *
 *  Created on: 22.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__SPACIALDISCRETIZATION__IMPL__
#define __H__LIBDISCRETIZATION__SPACIALDISCRETIZATION__IMPL__

namespace ug {

template <typename TDomain, typename TAlgebra>
SpacialDiscretization<TDomain, TAlgebra>::
SpacialDiscretization()
{
	m_name = "No Name";
}
template <typename TDomain, typename TAlgebra>
SpacialDiscretization<TDomain, TAlgebra>::
SpacialDiscretization(std::string name)
{
	m_name = name;
}

template <typename TDomain, typename TAlgebra>
void
SpacialDiscretization<TDomain, TAlgebra>::
set_name(std::string name)
{
	m_name = name;
}

template <typename TDomain, typename TAlgebra>
std::string
SpacialDiscretization<TDomain, TAlgebra>::
name()
{
	return m_name;
}

template <typename TDomain, typename TAlgebra>
void
SpacialDiscretization<TDomain, TAlgebra>::
add_SubsetDiscretization(SubsetDiscretization<d>& psd)
{
	m_SubsetDiscretization.push_back(&psd);
}

template <typename TDomain, typename TAlgebra>
void
SpacialDiscretization<TDomain, TAlgebra>::
delete_SubsetDiscretization(int nr)
{
	m_SubsetDiscretization.erase(m_SubsetDiscretization.begin() + nr);
}

template <typename TDomain, typename TAlgebra>
void
SpacialDiscretization<TDomain, TAlgebra>::
clear_SubsetDiscretizations()
{
	m_SubsetDiscretization.clear();
}

template <typename TDomain, typename TAlgebra>
int
SpacialDiscretization<TDomain, TAlgebra>::
numberOfSubsetDiscretizations()
{
	return m_SubsetDiscretization.size();
}

template <typename TDomain, typename TAlgebra>
void
SpacialDiscretization<TDomain, TAlgebra>::
print_info()
{
	std::cout << "SpacialDiscretization \"" << this->name() << "\" has the following member:" << std::endl;
	std::cout << this->numberOfSubsetDiscretizations() <<" SubsetDiscretization(s)" << std::endl;
}

}


#endif /* #ifndef __H__LIBDISCRETIZATION__SPACIALDISCRETIZATION__IMPL__ */
