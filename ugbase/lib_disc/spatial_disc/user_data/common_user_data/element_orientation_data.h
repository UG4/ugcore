/*
 * Copyright (c) 2021 - :  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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
 * User data of a subset indicator (1 in the subset, 0 everywhere else)
 */
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__LINE_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__LINE_USER_DATA__

#include <vector>

// ug4 headers
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"

namespace ug {

/// User data for the orientation of a line element
template <typename TDomain>
class EdgeOrientation
	: public StdUserData<EdgeOrientation<TDomain>, MathVector<TDomain::dim>, TDomain::dim>
{
public:
///	World dimension
	static constexpr int dim = TDomain::dim;

/// Return type
	using return_type = MathVector<dim>;

///	Type of domain
	using base_type = StdUserData<EdgeOrientation<TDomain>, MathVector<TDomain::dim>, TDomain::dim, void>;
	using cpl_user_data = CplUserData<MathVector<TDomain::dim>, TDomain::dim>;


	
public:

///	Constructor
	EdgeOrientation(SmartPtr<TDomain> domain) : m_spDomain(domain)
	{}

///	Indicator functions are discontinuous
	virtual bool continuous () const {return false;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return true;}

///	This function should not be used
	void operator() (return_type & vValue, const MathVector<dim> & globIP, number time, int si) const
	{ UG_THROW("SubsetIndicatorUserData: Element required for evaluation, but not passed. Cannot evaluate."); }

///	This function should not be used
	void operator() (return_type vValue [], const MathVector<dim> vGlobIP [], number time, int si, const size_t nip) const
	{ UG_THROW("SubsetIndicatorUserData: Element required for evaluation, but not passed. Cannot evaluate."); }


	inline void elem_evaluate(return_type &delta, GridObject * elem) const
	{
		auto& aaPos = m_spDomain->position_accessor();
		EdgeVertices *e = dynamic_cast<EdgeVertices*> (elem);
		UG_ASSERT(e != nullptr, "ERROR: No edge! ");

		VecSubtract(delta, aaPos[e->vertex(1)], aaPos[e->vertex(0)]); 	// v = e_1 - e_0
		delta /= VecLength(delta);										// normalize
	}



///	Evaluator
	template <int refDim>
	inline void evaluate
	(
		return_type vValue [],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		GridObject * elem,
		const MathVector<dim> vCornerCoords [],
		const MathVector<refDim> vLocIP [],
		const size_t nip,
		LocalVector * u,
		const MathMatrix<refDim, dim> * vJT = nullptr
	) const
	{

		MathVector<dim> delta;
		elem_evaluate(delta, elem);

		for (size_t i = 0; i < nip; i++)
		{
			vValue[i] = delta;
		}

	};
	
	///	implement as a UserData
			virtual void compute(LocalVector* u, GridObject* elem,
			                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
			{
				for(size_t s = 0; s < this->num_series(); ++s)
		   			for(size_t ip = 0; ip < this->num_ip(s); ++ip)
		   				elem_evaluate(this->value(s,ip), elem);
			}

		///	implement as a UserData
			virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
			                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
			{
				for(size_t s = 0; s < this->num_series(); ++s)
					for(size_t ip = 0; ip < this->num_ip(s); ++ip)
						elem_evaluate(this->value(s,ip), elem);
			}


protected:
	SmartPtr<TDomain> m_spDomain;
};

} // end namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__LINE_USER_DATA__

/* End of File */
