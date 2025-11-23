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

#include <vector>
#include <string>

// UG4 headers
#include "lib_grid/tools/subset_group.h"
#include "lib_disc/spatial_disc/user_data/linker/linker.h"

namespace ug{

/// This is a compositum for user data defined on different subsets
/**
 * This combines user data objects defined on several subsets which is handy,
 * but may be slow.
 */
template <typename TData, int dim, typename TRet = void>
class CompositeUserData : public UserData<TData, dim, TRet>
{
protected:
	using TCplUserData = CplUserData<TData, dim, TRet>;
public:
	using base_type = UserData<TData, dim, TRet>;
	using ref_type = SmartPtr<base_type>; ///< the attached UserData objects should have the same type as this class (i.e. they are "remapped")

	CompositeUserData() : m_bContinuous(true), m_bRequiresGridFunction(false) {}

	explicit CompositeUserData(bool continuous) : m_bContinuous(continuous), m_bRequiresGridFunction(false) {}

	~CompositeUserData() override = default;

	/// Add 'UserData' object for given subset index.
	void add
	(
		int si, ///< the subset index
		SmartPtr<base_type> ref ///< pointer to the user-data object
	)
	{
		UG_ASSERT (si >= 0, "CompositeUserData: Non-existing subset index!");
		
		if ((size_t) si >= m_vData.size ())
			m_vData.resize (si + 1);
		m_vData[si] = ref;
		
		// UG_ASSERT(ref->continuous() == m_continuous, "CompositeUserData: Mixing continuous and discontinuous data!");
		m_bContinuous = m_bContinuous && ref->continuous();
		m_bRequiresGridFunction = m_bRequiresGridFunction || ref->requires_grid_fct();
	}
	
	/// Add 'UserData' object for all subsets in a given group
	void add
	(
		const SubsetGroup & ssg, ///< the subset group
		SmartPtr<base_type> ref ///< pointer to the user-data object
	)
	{
		for (size_t i = 0; i < ssg.size (); i++) add (ssg[i], ref);
	}
	
	/// Add 'UserData' object for all subsets by their names
	void add
	(
		ConstSmartPtr<ISubsetHandler> ssh, ///< subset handler of the domain
		const char * ss_names, ///< names of the subdomains
		SmartPtr<base_type> ref ///< pointer to the user-data object
	)
	{
		std::vector<std::string> v_ss_names;
		SubsetGroup ssg (ssh);
		
		TokenizeTrimString (std::string (ss_names), v_ss_names);
		ssg.add (v_ss_names);
		add (ssg, ref);
	}

	/// Checks if anything is assigned to a given subset index
	bool has(int si) const  { return si >= 0 && (size_t) si < m_vData.size () && m_vData[si].valid ();}
	
	SmartPtr<base_type> get(int si) const { check (si); return m_vData[si]; }
	
	bool is_coupled(int si) { return has(si) && m_vData[si].template is_of_type<TCplUserData>(); }

	SmartPtr<TCplUserData> get_coupled(int si) { return m_vData[si].template cast_dynamic<TCplUserData>(); }

	// Implementing virtual functions

	bool continuous() const override {return m_bContinuous;}

	/// returns true, if at least one of the underlying UserData requires grid functions.
	bool requires_grid_fct() const override {
		 return m_bRequiresGridFunction;
	 }

	///	returns value for a global position
	TRet operator () (TData& value,
	                  const MathVector<dim>& globIP,
	                  number time, int si) const override
	{ check (si); return (*m_vData[si]) (value, globIP, time, si); }

	///	returns values for global positions
	void operator ()(TData vValue[],
	                 const MathVector<dim> vGlobIP[],
	                 number time, int si, const size_t nip) const override
	{ check (si); return (*m_vData[si]) (vValue, vGlobIP, time, si, nip); }


	void operator ()(TData vValue[],
	                 const MathVector<dim> vGlobIP[],
	                 number time, int si,
	                 GridObject* elem,
	                 const MathVector<dim> vCornerCoords[],
	                 const MathVector<1> vLocIP[],
	                 const size_t nip,
	                 LocalVector* u,
	                 const MathMatrix<1, dim>* vJT = nullptr) const override {
		check (si); return (*m_vData[si]) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
	}

	void operator () (TData vValue[],
	                  const MathVector<dim> vGlobIP[],
	                  number time, int si,
	                  GridObject* elem,
	                  const MathVector<dim> vCornerCoords[],
	                  const MathVector<2> vLocIP[],
	                  const size_t nip,
	                  LocalVector* u,
	                  const MathMatrix<2, dim>* vJT = nullptr) const override {
		check (si); return (*m_vData[si]) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
	}

	void operator () (TData vValue[],
	                  const MathVector<dim> vGlobIP[],
	                  number time, int si,
	                  GridObject* elem,
	                  const MathVector<dim> vCornerCoords[],
	                  const MathVector<3> vLocIP[],
	                  const size_t nip,
	                  LocalVector* u,
	                  const MathMatrix<3, dim>* vJT = nullptr) const override {
		check (si); return (*m_vData[si]) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
	}

private:

	// checks if the subset si is present in the list
	void check (int si) const
	{
		if (! has (si))
		{
			UG_THROW ("CompositeUserData: No data for subset " << si);
		}
	}

private:

	std::vector<SmartPtr<base_type> > m_vData;
	bool m_bContinuous;
	bool m_bRequiresGridFunction;
};


} // namespace ug

#endif
