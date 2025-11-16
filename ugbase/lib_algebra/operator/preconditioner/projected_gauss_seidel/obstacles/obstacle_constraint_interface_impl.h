/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__

#include "obstacle_constraint_interface.h"
#include "lib_disc/function_spaces/dof_position_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
init()
{
	if(m_spDomain.invalid())
		UG_THROW("No domain set in 'IObstacleConstraint::init' \n");

	if(m_spDD.invalid())
		UG_THROW("DofDistribution not set in 'IObstacleConstraint'.");

	//	build up a map of obstacle dofs and its corresponding obstacle values
	init_obstacle_dofs_with_values(1.0);

	UG_LOG("In IObstacleConstraint::init: "<<m_mObstacleValues.size()<< " obstacleDoFs tagged \n");
	UG_LOG("\n");
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim, bool> > func, const char* function)
{
	m_vCondNumberData.push_back(CondNumberData(func, function));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim, bool> > func, const char* function, const char* subsets)
{
	m_vCondNumberData.push_back(CondNumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim> > func, const char* function)
{
	m_vNumberData.push_back(NumberData(func, function));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim> > func, const char* function, const char* subsets)
{
	m_vNumberData.push_back(NumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(number value, const char* function)
{
	m_vConstNumberData.push_back(ConstNumberData(value, function));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(number value, const char* function, const char* subsets)
{
	m_vConstNumberData.push_back(ConstNumberData(value, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions)
{
	m_vVectorData.push_back(VectorData(func, functions));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions, const char* subsets)
{
	m_vVectorData.push_back(VectorData(func, functions, subsets));
}

#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(const char* name, const char* function)
{
	if(LuaUserData<number, dim>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim> > sp =
							LuaUserDataFactory<number, dim>::create(name);
		add(sp, function);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim, bool> > sp =
						LuaUserDataFactory<number, dim, bool>::create(name);
		add(sp, function);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(name)){
		SmartPtr<UserData<MathVector<dim>, dim> > sp =
				LuaUserDataFactory<MathVector<dim>, dim>::create(name);
		add(sp, function);
		return;
	}

//	no match found
	if(!CheckLuaCallbackName(name))
		UG_THROW("IObstacleConstraint::add: Lua-Callback with name '"<<name<<
		               "' does not exist.");

//	name exists but wrong signature
	UG_THROW("IObstacleConstraint::add: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(const char* name, const char* function, const char* subsets)
{
	if(LuaUserData<number, dim>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim> > sp =
							LuaUserDataFactory<number, dim>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim, bool> > sp =
						LuaUserDataFactory<number, dim, bool>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(name)){
		SmartPtr<UserData<MathVector<dim>, dim> > sp =
				LuaUserDataFactory<MathVector<dim>, dim>::create(name);
		add(sp, function, subsets);
		return;
	}

//	no match found
	if(!CheckLuaCallbackName(name))
		UG_THROW("IObstacleConstraint::add: Lua-Callback with name '"<<name<<
		               "' does not exist.");

//	name exists but wrong signature
	UG_THROW("IObstacleConstraint::add: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}
#endif


template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
clear()
{
	m_vCondNumberData.clear();
	m_vNumberData.clear();
	m_vConstNumberData.clear();
	m_vVectorData.clear();

	m_mObstacleValues.clear();
	m_vActiveDofs.clear();
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
check_functions_and_subsets(FunctionGroup& functionGroup, SubsetGroup& subsetGroup,
		size_t numFct) const
{
//	only number of functions allowed
	if(functionGroup.size() != numFct)
		UG_THROW("IObstacleConstraint:extract_data:"
					" Only "<<numFct<<" function(s) allowed in specification of a"
					" Obstacle Value, but the following functions given:"
					<<functionGroup);

//	get subsethandler
	ConstSmartPtr<ISubsetHandler> pSH = m_spDD->subset_handler();

// 	loop subsets
	for(size_t si = 0; si < subsetGroup.size(); ++si)
	{
	//	get subset index
		const int subsetIndex = subsetGroup[si];

	//	check that subsetIndex is valid
		if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
			UG_THROW("IObstacleConstraint:extract_data:"
							" Invalid Subset Index " << subsetIndex << ". (Valid is"
							" 0, .. , " << pSH->num_subsets() <<").");

	//	check all functions
		for(size_t i=0; i < functionGroup.size(); ++i)
		{
			const size_t fct = functionGroup[i];

		// 	check if function exist
			if(fct >= m_spDD->num_fct())
				UG_THROW("IObstacleConstraint:extract_data:"
							" Function "<< fct << " does not exist in pattern.");

		// 	check that function is defined for segment
			if(!m_spDD->is_def_in_subset(fct, subsetIndex))
				UG_THROW("IObstacleConstraint:extract_data:"
								" Function "<<fct<<" not defined on subset "<<subsetIndex);
		}
	}

//	everything ok
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TScheduledUserData>
void IObstacleConstraint<TDomain, TAlgebra>::
extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataObsSegment,
             std::vector<TScheduledUserData>& vUserData)
{
//	clear the extracted data
	mvUserDataObsSegment.clear();

	for(size_t i = 0; i < vUserData.size(); ++i)
	{
	//	create Function Group and Subset Group
		if (vUserData[i].bWholeDomain == false){
			try{
				vUserData[i].ssGrp = m_spDD->subset_grp_by_name(vUserData[i].ssName.c_str());
			}UG_CATCH_THROW(" Subsets '"<<vUserData[i].ssName<<"' not"
							" all contained in DoFDistribution.");
		}
		else{
			SubsetGroup ssGrp = SubsetGroup(m_spDD->subset_handler());
			ssGrp.add_all();
			vUserData[i].ssGrp = ssGrp;
		}

		FunctionGroup fctGrp;
		try{
			fctGrp = m_spDD->fct_grp_by_name(vUserData[i].fctName.c_str());
		}UG_CATCH_THROW(" Functions '"<<vUserData[i].fctName<<"' not"
		                " all contained in DoFDistribution.");

	//	check functions and subsets
		check_functions_and_subsets(fctGrp, vUserData[i].ssGrp, TUserData::numFct);

	//	set functions
		if(fctGrp.size() != TUserData::numFct)
			UG_THROW("IObstacleConstraint: wrong number of fct");

		for(size_t fct = 0; fct < TUserData::numFct; ++fct)
		{
			vUserData[i].fct[fct] = fctGrp[fct];
		}

	// 	loop subsets
		for(size_t si = 0; si < vUserData[i].ssGrp.size(); ++si)
		{
			//	get subset index and function
			const int subsetIndex = vUserData[i].ssGrp[si];

			//	remember functor and function
			mvUserDataObsSegment[subsetIndex].push_back(&vUserData[i]);
		}
	}
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
extract_data()
{
	extract_data(m_mNumberObsSegment, m_vNumberData);
	extract_data(m_mCondNumberObsSegment, m_vCondNumberData);
	extract_data(m_mConstNumberObsSegment, m_vConstNumberData);
	extract_data(m_mVectorObsSegment, m_vVectorData);
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
init_obstacle_dofs_with_values(number time)
{
	extract_data();

	//	reset map of obstacle values and vector of obstacle subset-indices
	m_mObstacleValues.clear();
	m_vObsSubsets.resize(0);

	init_obstacle_dofs_with_values<CondNumberData>(m_mCondNumberObsSegment, time);
	init_obstacle_dofs_with_values<NumberData>(m_mNumberObsSegment, time);
	init_obstacle_dofs_with_values<ConstNumberData>(m_mConstNumberObsSegment, time);
	init_obstacle_dofs_with_values<VectorData>(m_mVectorObsSegment, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void IObstacleConstraint<TDomain, TAlgebra>::
init_obstacle_dofs_with_values(const std::map<int, std::vector<TUserData*> >& mvUserData, number time)
{
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	store obstacle subsets
		m_vObsSubsets.push_back(si);

	//	get vector of scheduled obstacle data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	gets obstacle values in each base element type
		try
		{
		if(m_spDD->max_dofs(VERTEX))
			init_obstacle_dofs_with_values<RegularVertex, TUserData>(vUserData, si, time);
		if(m_spDD->max_dofs(EDGE))
			init_obstacle_dofs_with_values<Edge, TUserData>(vUserData, si, time);
		if(m_spDD->max_dofs(FACE))
			init_obstacle_dofs_with_values<Face, TUserData>(vUserData, si, time);
		if(m_spDD->max_dofs(VOLUME))
			init_obstacle_dofs_with_values<Volume, TUserData>(vUserData, si, time);
		}
		UG_CATCH_THROW("IObstacleConstraint::init_obstacle_dofs_with_values:"
						" While calling 'obstacle_value' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void IObstacleConstraint<TDomain, TAlgebra>::
init_obstacle_dofs_with_values(const std::vector<TUserData*>& vUserData, int si, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

	//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	using iter_type = typename DoFDistribution::traits<TBaseElem>::const_iterator;
	iter_type iter = m_spDD->begin<TBaseElem>(si);
	iter_type iterEnd = m_spDD->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get baseElem
		TBaseElem* elem = *iter;

	//	loop obstacle functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = m_spDD->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				m_spDD->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(),
						  "Mismatch: numInd="<<multInd.size()<<", numPos="
						  <<vPos.size()<<" on "<<elem->reference_object_id());

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is an obstacle fct and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					//	deposit obstacle values in a map
					m_mObstacleValues[multInd[j]] = val[f];
				}
			}
		}
	}
}

template <typename TDomain, typename TAlgebra>
bool
IObstacleConstraint<TDomain,TAlgebra>::
is_obs_dof(const DoFIndex& dof)
{
	if (m_mObstacleValues.find(dof) == m_mObstacleValues.end()){return false;}
	else {return true;}
}

template <typename TDomain, typename TAlgebra>
void
IObstacleConstraint<TDomain,TAlgebra>::
adjust_restriction(matrix_type& R, ConstSmartPtr<DoFDistribution> ddCoarse,
	ConstSmartPtr<DoFDistribution> ddFine, int type, number time)
{
	UG_LOG("IObstacleConstraint<TDomain,TAlgebra>::adjust_restrictionR \n");

	R.print();

	using iter_type = vector<DoFIndex>::iterator;
	iter_type dofIter = m_vActiveDofs.begin();
	iter_type dofIterEnd = m_vActiveDofs.end();
	for( ; dofIter != dofIterEnd; dofIter++)
	{
		UG_LOG("IObstacleConstraint<TDomain,TAlgebra>::"
				"adjust_restrictionR::activeDof : " <<*dofIter<< "\n");
		SetCol(R, (*dofIter)[0], (*dofIter)[1]);
	}

	if (m_vActiveDofs.size() > 0)
	{
		UG_LOG("#OfActiveDofs: " <<m_vActiveDofs.size()<< "\n");
		R.print();
	}
	UG_LOG("IObstacleConstraint::adjust_restrictionR() \n");
};

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__ */
