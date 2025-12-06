/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#include "neumann_boundary_base.h"
//ø #include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NeumannBoundaryBase<TDomain>::NeumannBoundaryBase(const char* function)
 :IElemDisc<TDomain>(function, "")
{
	if(this->num_fct() != 1)
		UG_THROW("NeumannBoundaryBase: needed exactly one function.");
}

////////////////////////////////////////////////////////////////////////////////
// User Data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
update_subset_groups(Data& userData)
{
//	create Function Group and Subset Group
	FunctionGroup functionGroup;

//	convert strings
	try{
		userData.InnerSSGrp = this->approx_space()->subset_grp_by_name(userData.InnerSubsetNames.c_str());
	}UG_CATCH_THROW("NeumannBoundaryBase:"
					" Subsets '"<<userData.InnerSubsetNames<<"' not"
					" all contained in ApproximationSpace.");
	try{
		userData.BndSSGrp = this->approx_space()->subset_grp_by_name(userData.BndSubsetNames.c_str());
	}UG_CATCH_THROW("NeumannBoundaryBase:"
					" Subsets '"<<userData.BndSubsetNames<<"' not"
					" all contained in ApproximationSpace.");
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add_inner_subsets(const char* InnerSubsets)
{
	std::vector<std::string> vSubsets = this->symb_subsets();
	std::vector<std::string> vNew = TokenizeTrimString(InnerSubsets);
	for(size_t i = 0; i < vNew.size(); ++i)
		if(std::find(vSubsets.begin(), vSubsets.end(), vNew[i]) == vSubsets.end())
			vSubsets.push_back(vNew[i]);
	this->set_subsets(vSubsets);

}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(SmartPtr<CplUserData<number, dim> > data, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}

	add(data, bnd.c_str(), inner.c_str());
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(SmartPtr<CplUserData<number, dim, bool> > data, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}

	add(data, bnd.c_str(), inner.c_str());
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(SmartPtr<CplUserData<MathVector<dim>, dim> > data, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}

	add(data, bnd.c_str(), inner.c_str());
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(number val, const char* function, const char* subsets)
{
	SmartPtr<CplUserData<number, dim> > sp = make_sp(new ConstUserNumber<dim>(val));
	add(sp, function, subsets);
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(number val, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}

	add(val, bnd.c_str(), inner.c_str());
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(const std::vector<number>& val, const char* function, const char* subsets)
{
	SmartPtr<CplUserData<MathVector<dim>, dim> > sp = make_sp(new ConstUserVector<dim>(val));
	add(sp, function, subsets);
}

template<typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(const std::vector<number>& val, const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}

	add(val, bnd.c_str(), inner.c_str());
}


#ifdef UG_FOR_LUA
template <typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(const char* name, const char* function, const char* subsets)
{
	if(LuaUserData<number, dim>::check_callback_returns(name)){
		SmartPtr<CplUserData<number, dim> > sp =
							LuaUserDataFactory<number, dim>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(name)){
		SmartPtr<CplUserData<number, dim, bool> > sp =
				LuaUserDataFactory<number, dim, bool>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(name)){
		SmartPtr<CplUserData<MathVector<dim>, dim> > sp =
				LuaUserDataFactory<MathVector<dim>, dim>::create(name);
		add(sp, function, subsets);
		return;
	}

//	no match found
	if(!CheckLuaCallbackName(name))
		UG_THROW("NeumannBoundaryBase: Lua-Callback with name '"<<name<<
		               "' does not exist.");

//	name exists but wrong signature
	UG_THROW("NeumannBoundaryBase: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}

template <typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(const char* name,  const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}
	add(name, bnd.c_str(), inner.c_str());
}

template <typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(LuaFunctionHandle fct, const char* function, const char* subsets)
{
	if(LuaUserData<number, dim>::check_callback_returns(fct)){
		SmartPtr<CplUserData<number, dim> > sp =
							make_sp(new LuaUserData<number, dim>(fct));
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(fct)){
		SmartPtr<CplUserData<number, dim, bool> > sp =
							make_sp(new LuaUserData<number, dim, bool>(fct));
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(fct)){
		SmartPtr<CplUserData<MathVector<dim>, dim> > sp =
							make_sp(new LuaUserData<MathVector<dim>, dim>(fct));
		add(sp, function, subsets);
		return;
	}

//	name exists but wrong signature
	UG_THROW("NeumannBoundaryBase: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}

template <typename TDomain>
void NeumannBoundaryBase<TDomain>::
add(LuaFunctionHandle fct,  const std::vector<std::string>& BndSubsets, const std::vector<std::string>& InnerSubsets)
{
	std::string bnd;
	for(size_t i = 0; i < BndSubsets.size(); ++i){
		if(i > 0) bnd.append(",");
		bnd.append(BndSubsets[i]);
	}
	std::string inner;
	for(size_t i = 0; i < InnerSubsets.size(); ++i){
		if(i > 0) inner.append(",");
		inner.append(InnerSubsets[i]);
	}
	add(fct, bnd.c_str(), inner.c_str());
}

#endif

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NeumannBoundaryBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NeumannBoundaryBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NeumannBoundaryBase<Domain3d>;
#endif

} // namespace ug

