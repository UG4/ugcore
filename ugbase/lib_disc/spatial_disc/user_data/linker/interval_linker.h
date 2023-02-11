/*
 * Copyright (c) 2013-2022:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/*
 * A linker that cuts out values of a given userdata objection in some interval.
 */
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__INTERVAL_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__INTERVAL_LINKER__

/* ug4 headers */
#include "common/common.h"

#include "linker.h"

namespace ug {

/**
 * User data filtering out values of a given argument beyond a given interval. Out of the
 * interval, the return value is "def_val" (0 by default)
 */
template <int dim>
class IntervalNumberLinker
	: public StdDataLinker<IntervalNumberLinker<dim>, number, dim>
{
///	Base class type
	typedef StdDataLinker<IntervalNumberLinker<dim>, number, dim> base_type;
	
public:

///	Constructor
	IntervalNumberLinker
	(
		SmartPtr<CplUserData<number, dim> > spData, ///< data to filter
		MathVector<dim>& minCoord, ///< left boundaries of the interval
		MathVector<dim>& maxCoord ///< right boundaries of the interval
	)
	{
		init (spData, minCoord, maxCoord);
	}
	
///	Constructor
	IntervalNumberLinker
	(
		SmartPtr<CplUserData<number, dim> > spData, ///< data to filter
		std::vector<number> min_coord, ///< left boundaries of the interval
		std::vector<number> max_coord ///< right boundaries of the interval
	)
	{
		if (min_coord.size () != dim || max_coord.size () != dim)
			UG_THROW ("IntervalNumberLinker: Illegal sizes of the boundar arrays!");
		MathVector<dim> minCoord, maxCoord;
		for (size_t i = 0; i < dim; i++)
		{
			minCoord[i] = min_coord[i]; maxCoord[i] = max_coord[i];
		}
		
		init (spData, minCoord, maxCoord);
	}
	
///	sets the default values out of the interval
	void set_default (number v) {def_val = v;}
	
///	Returns true because without a grid function, we do not get the element
	virtual bool requires_grid_fct() const {return true;}

///	Evaluation for the global coordinates
	inline void evaluate
	(
		number& value,
		const MathVector<dim>& glob_ip,
		number time,
		int si
	) const
	{
		if (is_in (glob_ip))
			(* m_spData) (value, glob_ip, time, si);
		else
			value = def_val;
	}

///	Computation without the derivatives
	template <int refDim>
	inline void evaluate
	(
		number vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
	//	Compute all, then replace:
	
		(*m_spData) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		
		for (size_t ip = 0; ip < nip; ip++)
			if (! is_in (vGlobIP[ip]))
				vValue[ip] = def_val;
	}
	
///	Computation of the values and the derivatives
	template <int refDim>
	void eval_and_deriv
	(
		number vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		bool bDeriv,
		int s,
		std::vector<std::vector<number> > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		if (this->zero_derivative ())
			bDeriv = false;
		else
			this->set_zero (vvvDeriv, nip);
		
		const number* vValues = m_spData->values (s);
		
		for (size_t ip = 0; ip < nip; ip++)
			if (is_in (vGlobIP[ip]))
			{
				vValue[ip] = vValues[ip];
				if (bDeriv)
					for (size_t fct = 0; fct < m_spDData->num_fct (); fct++)
					{
						const number* vDValues = m_spDData->deriv (s, ip, fct);
						const size_t c_fct = this->input_common_fct (0, fct);
						for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
							vvvDeriv[ip][c_fct][sh] = vDValues [sh];
					}
			}
			else
				vValue[ip] = def_val;
				// The derivatives are all initialized with 0.
	}
	
private:

///	a general initializer (to call from a constructor)
	void init
	(
		SmartPtr<CplUserData<number, dim> > spData, ///< data to filter
		MathVector<dim>& minCoord, ///< left boundaries of the interval
		MathVector<dim>& maxCoord ///< right boundaries of the interval
	)
	{
		this->set_num_input (1);
		m_spData = spData;
		m_spDData = spData.template cast_dynamic<DependentUserData<number, dim> > ();
		this->set_input (0, spData, spData);
		m_minCoord = minCoord; m_maxCoord = maxCoord;
		
		def_val = 0;
	}

///	checks if the point is in the interval
	bool is_in
	(
		const MathVector<dim>& x // the point
	) const
	{
		for (size_t i = 0; i < dim; i++)
			if (x[i] < m_minCoord[i] || x[i] > m_maxCoord[i])
				return false;
		return true;
	}

///	data to filter
	SmartPtr<CplUserData<number, dim> > m_spData;
	SmartPtr<DependentUserData<number, dim> > m_spDData;
	
///	the interval
	MathVector<dim> m_minCoord, m_maxCoord;
	
//	the default value
	number def_val;
};

} // end namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__INTERVAL_LINKER__

/* End of File */
