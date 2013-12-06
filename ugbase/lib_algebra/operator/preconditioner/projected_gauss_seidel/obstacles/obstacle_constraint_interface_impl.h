/*
 * obstacle_constraint_interface_impl.h
 *
 *  Created on: 26.11.2013
 *      Author: raphaelprohl
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
initObsValues()
{
	clear();

	if(m_spDomain.invalid())
		UG_THROW("No domain set in 'IObstacleConstraint::init' \n");

	if(m_spDD.invalid())
		UG_THROW("DofDistribution not set in 'IObstacleConstraint'.");

	obstacle_value(1.0);
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
clear()
{
	m_vCondNumberData.clear();
	m_vNumberData.clear();
	m_vConstNumberData.clear();
	m_vVectorData.clear();

	m_vDoFsOnObsSubset.clear();
	m_mObsValues.clear();
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim, bool> > func, const char* function, const char* subsets)
{
	m_vCondNumberData.push_back(CondNumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim> > func, const char* function, const char* subsets)
{
	m_vNumberData.push_back(NumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void IObstacleConstraint<TDomain, TAlgebra>::
add(number value, const char* function, const char* subsets)
{
	m_vConstNumberData.push_back(ConstNumberData(value, function, subsets));
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
		try{
			vUserData[i].ssGrp = m_spDD->subset_grp_by_name(vUserData[i].ssName.c_str());
		}UG_CATCH_THROW(" Subsets '"<<vUserData[i].ssName<<"' not"
		                " all contained in DoFDistribution.");

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
obstacle_value(number time)
{
	extract_data();

	//	reset vector of indices in obstacle-subsets
	m_vDoFsOnObsSubset.resize(0);

	obstacle_value<CondNumberData>(m_mCondNumberObsSegment, time);
	obstacle_value<NumberData>(m_mNumberObsSegment, time);
	obstacle_value<ConstNumberData>(m_mConstNumberObsSegment, time);
	obstacle_value<VectorData>(m_mVectorObsSegment, time);

	//	sort m_vDoFsOnObsSubset
	sort(m_vDoFsOnObsSubset.begin(), m_vDoFsOnObsSubset.end());

	for (vector<DoFIndex>::iterator itObsInd = m_vDoFsOnObsSubset.begin();
		itObsInd < m_vDoFsOnObsSubset.end(); ++itObsInd)
	{
		UG_LOG("dofIndex: "<<(*itObsInd)<<"\n");
		UG_LOG("\n");
		UG_LOG("\n");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void IObstacleConstraint<TDomain, TAlgebra>::
obstacle_value(const std::map<int, std::vector<TUserData*> >& mvUserData, number time)
{
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled obstacle data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	gets obstacle values in each base element type
		try
		{
		if(m_spDD->max_dofs(VERTEX))
			obstacle_value<Vertex, TUserData>(vUserData, si, time);
		if(m_spDD->max_dofs(EDGE))
			obstacle_value<EdgeBase, TUserData>(vUserData, si, time);
		if(m_spDD->max_dofs(FACE))
			obstacle_value<Face, TUserData>(vUserData, si, time);
		if(m_spDD->max_dofs(VOLUME))
			obstacle_value<Volume, TUserData>(vUserData, si, time);
		}
		UG_CATCH_THROW("IObstacleConstraint::obstacle_value:"
						" While calling 'obstacle_value' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void IObstacleConstraint<TDomain, TAlgebra>::
obstacle_value(const std::vector<TUserData*>& vUserData, int si, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

	//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator iter_type;
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
			//	store the multi indices in a vector
				m_spDD->inner_dof_indices(elem, fct, m_vDoFsOnObsSubset, false);

				UG_ASSERT(multInd.size() == vPos.size(),
						  "Mismatch: numInd="<<multInd.size()<<", numPos="
						  <<vPos.size()<<" on "<<elem->reference_object_id());

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is an obstacle fct and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					//	deposit obstacle values in a map
					m_mObsValues[multInd[j]] = val[f];
				}
			}
		}
	}
}

template <typename TDomain, typename TAlgebra>
bool
IObstacleConstraint<TDomain,TAlgebra>::
dof_lies_on_obs_subset(const DoFIndex& dof, const vector<DoFIndex>& vOfDofs)
{
	//	TODO: it would be more efficient to attach a flag to every dof
	//	which indicates whether the dof is a dof in an obstacle subset or not
	typedef vector<DoFIndex>::const_iterator iter_type;
	iter_type iter = vOfDofs.begin();
	iter_type iterEnd =	vOfDofs.end();

	for( ; iter != iterEnd; iter++)
	{
		if (*iter == dof)
			return true;
	}

	return false;
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__ */
