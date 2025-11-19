/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_BASE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_BASE__

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

template<typename TDomain>
class NeumannBoundaryBase
	: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		using base_type = IElemDisc<TDomain>;

	///	Base class type
		using this_type = NeumannBoundaryBase;

	public:
	///	World dimension
		static constexpr int dim = base_type::dim;

	public:
	///	default constructor
		NeumannBoundaryBase(const char* function);

	///	add a boundary value
	///	\{
		virtual void add(SmartPtr<CplUserData<number, dim> > data, const char* BndSubsets, const char* InnerSubsets) = 0;
		void add(SmartPtr<CplUserData<number, dim> > data, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);
		virtual void add(SmartPtr<CplUserData<number, dim, bool> > user, const char* BndSubsets, const char* InnerSubsets) = 0;
		void add(SmartPtr<CplUserData<number, dim, bool> > user, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);
		virtual void add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* BndSubsets, const char* InnerSubsets) = 0;
		void add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);

		void add(number val, const char* BndSubsets, const char* InnerSubsets);
		void add(number val, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);
		void add(const std::vector<number>& val, const char* BndSubsets, const char* InnerSubsets);
		void add(const std::vector<number>& val, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);
#ifdef UG_FOR_LUA
		void add(const char* name, const char* BndSubsets, const char* InnerSubsets);
		void add(const char* name, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);
		void add(LuaFunctionHandle fct, const char* BndSubsets, const char* InnerSubsets);
		void add(LuaFunctionHandle fct, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets);
#endif
	/// \}

	protected:
	///	base class for user data
		struct Data
		{
			Data(std::string BndSubsets_, std::string InnerSubsets_)
							: BndSubsetNames(BndSubsets_), InnerSubsetNames(InnerSubsets_) {}
			SubsetGroup BndSSGrp;
			std::string BndSubsetNames;
			SubsetGroup InnerSSGrp;
			std::string InnerSubsetNames;
		};

	///	method used to extract subsets id
		void update_subset_groups(Data& userData);

	///	adds subsets to the looped inner subsets
		void add_inner_subsets(const char* InnerSubsets);

	public:
	///	 returns the type of elem disc
		virtual int type() const {return EDT_BND;}

	protected:
	///	dummy add methods
	///	\{
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template<typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template<typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
	/// \}
};

} // end namespac ug

#endif