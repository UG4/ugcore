/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_TRAITS__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_TRAITS__

#include "common/common.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Linker Traits to allow generic programming
////////////////////////////////////////////////////////////////////////////////

/// Linker Traits
template <typename TData, typename TDataIn, typename TRet = TData>
struct linker_traits
{
///	computes out += s * in1 (with appropriate '*')
	static void mult_add(TData& out, const TData& in1, const TDataIn& s);
};

template <>
struct linker_traits<number, number>
{
	static void mult_add(number& out, const number& in1, const number& s)
	{
		out += in1 * s;
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>, number >
{
	static void mult_add(MathVector<dim>& out,
	                     const MathVector<dim>& in1,
						 const number& s)
	{
		VecScaleAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>, MathMatrix<dim,dim> >
{
	static void mult_add(MathVector<dim>& out,
	                     const MathVector<dim>& in1,
						 const MathMatrix<dim,dim>& s)
	{
		MatVecMultAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>,  MathVector<dim>, number >
{
	static void mult_add(number& out,
	                     const MathVector<dim>& in1,
						 const MathVector<dim>& s)
	{
		out = VecDot(s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathMatrix<dim,dim>, number >
{
	static void mult_add(MathMatrix<dim,dim>& out,
	                     const MathMatrix<dim,dim>& in1,
						 const number& s)
	{
		MatScaleAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathTensor<4,dim>, number >
{
	static void mult_add(MathTensor<4,dim>& out,
	                     const MathTensor<4, dim>& in1,
						 const number& s)
	{
		out = in1;
		UG_THROW("linker_traits for MathTensor4 not implemented");
	}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_TRAITS__ */
