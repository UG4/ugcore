/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__ENERGY_CONVERGENCE_CHECK__
#define __H__LIB_ALGEBRA__OPERATOR__ENERGY_CONVERGENCE_CHECK__

#include "convergence_check.h"

namespace ug{


template <typename TVector>
class EnergyConvCheck : public StdConvCheck<TVector>
{
	using base_type = StdConvCheck<TVector>;
	public:
	EnergyConvCheck() : base_type() {}
	EnergyConvCheck(int maxSteps, number minDefect, number relReduction) : base_type(maxSteps, minDefect, relReduction) {}
	EnergyConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose)
	 	 : base_type(maxSteps, minDefect, relReduction, verbose) {}
	EnergyConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose,bool suppressUnsuccessful)
	 	 : base_type(maxSteps, minDefect, relReduction, verbose, suppressUnsuccessful) {}

	~EnergyConvCheck() override = default;

	SmartPtr<IConvergenceCheck<TVector> > clone() override
		{
		SmartPtr<EnergyConvCheck > newInst(new EnergyConvCheck);
		// use std assignment (implicit member-wise is fine here)
		*newInst = *this;
		return newInst;
	}

	void start(const TVector& d) override
	{
		base_type::start_defect(energy_norm(d));
	}
	void update(const TVector& d) override
	{
		base_type::update_defect(energy_norm(d));
	}

	double energy_norm(const TVector &d)
	{
		if(tmp.valid() == false || tmp->size() != d.size())
		{
			tmp = d.clone_without_values();
			tmp2 = d.clone_without_values();
		}
		TVector &t = *tmp;
		TVector &t2 = *tmp2;
		t = d;
#ifdef UG_PARALLEL
		t.change_storage_type(PST_CONSISTENT);
#endif
		m_op->apply(t2, t);
		return sqrt(VecProd(t, t2));
	}

	std::string config_string() const override
		{
		std::stringstream ss;
		ss << "EnergyConvCheck( max steps = " << base_type::m_maxSteps << ", min defect = " << base_type::m_minDefect <<
				", relative reduction = " << base_type::m_relReduction << ")";
		return ss.str();
	}

	void set_linear_operator(SmartPtr<ILinearOperator<TVector> > op)
	{
		m_op = op;
	}

private:
	SmartPtr<TVector> tmp, tmp2;
	SmartPtr<ILinearOperator<TVector> > m_op;

};


}

#endif