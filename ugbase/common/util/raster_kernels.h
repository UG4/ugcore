/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_raster_kernels
#define __H__UG_raster_kernels

#include "raster.h"

namespace ug {
namespace raster_kernels {

///	Kernel which counts the number of times it was run on valid data values
/** For use with ug::Raster.
 * 
 * This class defines a default constructor, the type 'result_t',
 * and the method 'result_t result() const' and can thus be used like this:
 *
 * \code
 * Raster<T,TDIM> raster;
 * //...
 * size_t countAll = raster.run_on_all<Count<T,TDIM> > ();
 * size_t countNbrs = raster.run_on_nbrs<Count<T,TDIM> > (someMultiIndex);
 * \endcode
 *
 * \note	if the kernel is run on a 'no_data_value', the counter is not increased.
 */
template <typename T, int TDIM>
class Count {
public:
	using result_t = size_t;

	Count () :
		m_count (0)
	{}

	inline void operator () (Raster<T, TDIM>& raster,
							 const typename Raster<T, TDIM>::MultiIndex& cur)
	{
		if(raster.node_value(cur) != raster.no_data_value())
			++m_count;
	}

	inline size_t result() const
	{return m_count;}

private:
	size_t m_count;
};


///	Kernel which sums the values for all entries it was called on
/** For use with ug::Raster.
 * 
 * This class defines a default constructor, the type 'result_t',
 * and the method 'result_t result() const' and can thus be used like this:
 *
 * \code
 * Raster<T,TDIM> raster;
 * //...
 * size_t sumAll = raster.run_on_all<Sum<T,TDIM> > ();
 * size_t sumNbrs = raster.run_on_nbrs<Sum<T,TDIM> > (someMultiIndex);
 * \endcode
 *
 * \note	'no_data_values' are ignored and will not be summed.
 */
template <typename T, int TDIM>
class Sum {
public:
	using result_t = T;

	Sum () :
		m_sum (0)
	{}

	inline void operator () (Raster<T, TDIM>& raster,
							 const typename Raster<T, TDIM>::MultiIndex& cur)
	{
		if(raster.node_value(cur) != raster.no_data_value())
			m_sum += raster.node_value(cur);
	}

	inline T result() const
	{return m_sum;}

private:
	T m_sum;
};


///	Kernel which blurs all values of a raster it was called on
/** \note	This class does not feature a default constructor and thus
 *			must be constructed before being run on a raster.
 *
 * Use, e.g., like this:
 * \code
 * Raster<T,TDIM> raster;
 * //...
 * Blur<T,TDIM> blurKernel(0.1);
 * raster.run_on_all (blurKernel);
 * \endcode
 *
 * \note 'no_data_values' will not be affected by the blur operation
 */
template <typename T, int TDIM>
class Blur {
public:
	Blur (T alpha) :
		m_alpha (alpha)
	{}

	void operator () (Raster<T, TDIM>& raster,
	            	  const typename Raster<T, TDIM>::MultiIndex& cur)
	{
		if(raster.node_value(cur) != raster.no_data_value()) {
			T numNbrs = (T)raster.template run_on_nbrs<Count<T, TDIM> >(cur);
			if(numNbrs) {
				T nbrVal = raster.template run_on_nbrs<Sum<T, TDIM> >(cur);
				
				raster.node_value(cur) = raster.node_value(cur) * (1. - m_alpha)
									   + (nbrVal * m_alpha / numNbrs);
			}
		}
	}

private:
	T m_alpha;
};

}//	end of namespace
}//	end of namespace

#endif