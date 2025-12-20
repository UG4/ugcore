/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__PIECEWISE_CONSTANT__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__PIECEWISE_CONSTANT__

#include "common/util/provider.h"
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
// #include "lib_disc/reference_element/reference_element_util.h"

namespace ug {

/// Elementwise constant shape functions
template <typename TRefElem>
class PiecewiseConstantLSFS
	: public BaseLSFS<PiecewiseConstantLSFS<TRefElem>, TRefElem::dim>
{
	public:
	///	Dimension, where shape functions are defined
		static constexpr int dim = TRefElem::dim;

	public:
	///	Constructor
		PiecewiseConstantLSFS()
		{
			const TRefElem& rRef = Provider<TRefElem>::get();

			bary = rRef.corner(0);
            for (size_t j=1; j < rRef.num(0); ++j){
            	bary += rRef.corner(j);
            }
            bary *= 1./rRef.num(0);

            m_vLocalDoF = LocalDoF(dim, 0, 0);
		}

	public:
	///	\copydoc ug::DimLocalDoFSet::roid()
		[[nodiscard]] ReferenceObjectID_t roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	\copydoc ug::DimLocalDoFSet::num_dof()
		[[nodiscard]] size_t num_dof() const {return 1;};

	///	\copydoc ug::DimLocalDoFSet::num_dof()
		[[nodiscard]] size_t num_dof(ReferenceObjectID_t type) const
		{
			if (type == TRefElem::REFERENCE_OBJECT_ID)   return 1;
			else return 0;
		}

	///	\copydoc ug::DimLocalDoFSet::local_dof()
		[[nodiscard]] const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF;}

	///	\copydoc ug::DimLocalDoFSet::position()
		[[nodiscard]] inline bool position(size_t i, MathVector<dim>& pos) const
		{
			pos = bary; return true;
		}

	///	\copydoc ug::DimLocalDoFSet::exact_position_available()
		[[nodiscard]] bool exact_position_available() const {return true;};

	public:
	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		[[nodiscard]] bool continuous() const {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		[[nodiscard]] size_t num_sh() const {return 1;}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		[[nodiscard]] number shape(const size_t i, const MathVector<dim>& x) const
		{
			return 1;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			TRefElem::check_position(x);
			VecSet(g, 0.0);
		}

	protected:
		MathVector<dim> bary; ///< Barycenter
		LocalDoF m_vLocalDoF; ///< association to elements
};

} //namespace ug

#endif