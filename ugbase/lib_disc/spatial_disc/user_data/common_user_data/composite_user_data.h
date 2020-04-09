/*
 * Copyright (c) 2019-2020:  G-CSC, Goethe University Frankfurt
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

#ifndef __LIB_DISC__COMPOSITE_USER_DATA_H_
#define __LIB_DISC__COMPOSITE_USER_DATA_H_

// UG4
#include "lib_disc/spatial_disc/user_data/linker/linker.h"

namespace ug{

//! This is a compositum for user data from different subsets.
/*! (This is handy, but the implementation is rather slow.) */
template <typename TData, int dim, typename TRet = void>
class CompositeUserData : public UserData<TData, dim, TRet>
{
public:
	typedef UserData<TData, dim, TRet> base_type;
	typedef SmartPtr<base_type> ref_type;
	typedef std::map<int, ref_type>  map_type;

	CompositeUserData(bool continuous) : m_continuous(continuous), m_bRequiresGridFunction(false)
	{}

	virtual ~CompositeUserData(){}

	///! Add 'UserData' for given subset index.
	void add(int si, ref_type ref)

	{
		// UG_ASSERT(ref->continuous() == m_continuous, "ERROR: Mixing continuous and discontinuous data!");
		m_map[si] = ref;
		m_continuous = m_continuous && ref->continuous();
		m_bRequiresGridFunction = m_bRequiresGridFunction || ref->requires_grid_fct();
	}

	// Implementing virtual functions

	virtual bool continuous() const
	{return m_continuous;}


	//! returns true, if at least one of the underlying UserData requires grid functions.
	 virtual bool requires_grid_fct() const
	 {
		 return m_bRequiresGridFunction;
	 }

		///	returns value for a global position
			virtual TRet operator() (TData& value,
									 const MathVector<dim>& globIP,
									 number time, int si) const
			{ return (*find(si)->second)(value, globIP, time, si); }

		///	returns values for global positions
			virtual void operator()(TData vValue[],
									const MathVector<dim> vGlobIP[],
									number time, int si, const size_t nip) const
			{ return (*find(si)->second)(vValue, vGlobIP, time, si, nip); }


			virtual void operator()(TData vValue[],
					                        const MathVector<dim> vGlobIP[],
					                        number time, int si,
					                        GridObject* elem,
					                        const MathVector<dim> vCornerCoords[],
					                        const MathVector<1> vLocIP[],
					                        const size_t nip,
					                        LocalVector* u,
					                        const MathMatrix<1, dim>* vJT = NULL) const
			{
				return (*find(si)->second)(vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			}

			virtual void operator()(TData vValue[],
					                        const MathVector<dim> vGlobIP[],
					                        number time, int si,
					                        GridObject* elem,
					                        const MathVector<dim> vCornerCoords[],
					                        const MathVector<2> vLocIP[],
					                        const size_t nip,
					                        LocalVector* u,
					                        const MathMatrix<2, dim>* vJT = NULL) const
			{
				return (*find(si)->second)(vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			}

			virtual void operator()(TData vValue[],
			                       const MathVector<dim> vGlobIP[],
			                        number time, int si,
			                        GridObject* elem,
			                        const MathVector<dim> vCornerCoords[],
			                        const MathVector<3> vLocIP[],
			                        const size_t nip,
			                        LocalVector* u,
			                        const MathMatrix<3, dim>* vJT = NULL) const
			{
				return (*find(si)->second)(vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			}



protected:
	typename map_type::const_iterator find(int si) const
	{
		typename map_type::const_iterator it = m_map.find(si);
		UG_ASSERT(it != m_map.end(), "ERROR:");
		return it;
	}

	map_type m_map;
	bool m_continuous;
	bool m_bRequiresGridFunction;
};


} // namespace ug

#endif /* __LIB_DISC__COMPOSITE_USER_DATA_H_ */
