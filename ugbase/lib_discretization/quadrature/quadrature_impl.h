/*
 * quadrature_impl.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__QUADRATURE_IMPL__
#define __H__LIBDISCRETIZATION__QUADRATURE_IMPL__

namespace ug{


template <typename TRefElem>
GaussQuadrature<TRefElem>::~GaussQuadrature()
{
	if(this->m_points != 0) delete[] this->m_points;
	if(this->m_weights != 0) delete[] this->m_weights;
}

template <typename TRefElem>
inline bool GaussQuadrature<TRefElem>::allocate_memory(std::size_t n)
{
	this->m_points = new typename QuadratureRule<TRefElem>::position_type[this->m_num_points];
	this->m_weights = new typename QuadratureRule<TRefElem>::weight_type[this->m_num_points];
	if(this->m_points == 0 || this->m_weights != 0)
		return false;
	return true;
}

}


#endif /* __H__LIBDISCRETIZATION__QUADRATURE_IMPL__ */
