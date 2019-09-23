/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__

#include "lagrange_dirichlet_boundary.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/dof_position_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
//	setup
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	base_type::set_approximation_space(approxSpace);
	m_spApproxSpace = approxSpace;
	m_spDomain = approxSpace->domain();
	m_aaPos = m_spDomain->position_accessor();
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
clear()
{
	m_vBNDNumberData.clear();
	m_vNumberData.clear();
	m_vConstNumberData.clear();
	m_vVectorData.clear();
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim, bool> > func, const char* function, const char* subsets)
{
	m_vBNDNumberData.push_back(CondNumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim, bool> > func, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(func, function.c_str(), subsets.c_str());
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim> > func, const char* function, const char* subsets)
{
	m_vNumberData.push_back(NumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim> > func, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(func, function.c_str(), subsets.c_str());
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(number value, const char* function, const char* subsets)
{
	m_vConstNumberData.push_back(ConstNumberData(value, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(number value, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(value, function.c_str(), subsets.c_str());
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions, const char* subsets)
{
	m_vVectorData.push_back(VectorData(func, functions, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > func, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(func, function.c_str(), subsets.c_str());
}

#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
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
		UG_THROW("LagrangeDirichlet::add: Lua-Callback with name '"<<name<<
		               "' does not exist.");

//	name exists but wrong signature
	UG_THROW("LagrangeDirichlet::add: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(const char* name, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(name, function.c_str(), subsets.c_str());
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(LuaFunctionHandle fct, const char* function, const char* subsets)
{
	if(LuaUserData<number, dim>::check_callback_returns(fct)){
		SmartPtr<UserData<number, dim> > sp =
				make_sp(new LuaUserData<number, dim>(fct));
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(fct)){
		SmartPtr<UserData<number, dim, bool> > sp =
				make_sp(new LuaUserData<number, dim, bool>(fct));
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(fct)){
		SmartPtr<UserData<MathVector<dim>, dim> > sp =
				make_sp(new LuaUserData<MathVector<dim>, dim>(fct));
		add(sp, function, subsets);
		return;
	}

//	name exists but wrong signature
	UG_THROW("LagrangeDirichlet::add: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(LuaFunctionHandle fct, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(fct, function.c_str(), subsets.c_str());
}
#endif

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(const char* functions, const char* subsets)
{
	m_vOldNumberData.push_back(OldNumberData(functions, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets)
{
	std::string function;
	for(size_t i = 0; i < Fcts.size(); ++i){
		if(i > 0) function.append(",");
		function.append(Fcts[i]);
	}
	std::string subsets;
	for(size_t i = 0; i < Subsets.size(); ++i){
		if(i > 0) subsets.append(",");
		subsets.append(Subsets[i]);
	}

	add(function.c_str(), subsets.c_str());
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
check_functions_and_subsets(FunctionGroup& functionGroup, SubsetGroup& subsetGroup, size_t numFct) const
{
//	only number of functions allowed
	if(functionGroup.size() != numFct)
		UG_THROW("DirichletBoundary:extract_data:"
					" Only "<<numFct<<" function(s) allowed in specification of a"
					" Dirichlet Value, but the following functions given:"
					<<functionGroup);

//	get subsethandler
	ConstSmartPtr<ISubsetHandler> pSH = m_spApproxSpace->subset_handler();

// 	loop subsets
	for(size_t si = 0; si < subsetGroup.size(); ++si)
	{
	//	get subset index
		const int subsetIndex = subsetGroup[si];

	//	check that subsetIndex is valid
		if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
			UG_THROW("DirichletBoundary:extract_data:"
							" Invalid Subset Index " << subsetIndex << ". (Valid is"
							" 0, .. , " << pSH->num_subsets() <<").");

	//	check all functions
		for(size_t i=0; i < functionGroup.size(); ++i)
		{
			const size_t fct = functionGroup[i];

		// 	check if function exist
			if(fct >= m_spApproxSpace->num_fct())
				UG_THROW("DirichletBoundary:extract_data:"
							" Function "<< fct << " does not exist in pattern.");

		// 	check that function is defined for segment
			if(!m_spApproxSpace->is_def_in_subset(fct, subsetIndex))
				UG_THROW("DirichletBoundary:extract_data:"
								" Function "<<fct<<" not defined on subset "<<subsetIndex);
		}
	}

//	everything ok
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TScheduledUserData>
void DirichletBoundary<TDomain, TAlgebra>::
extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
             std::vector<TScheduledUserData>& vUserData)
{
//	clear the extracted data
	mvUserDataBndSegment.clear();

	for(size_t i = 0; i < vUserData.size(); ++i)
	{
	//	create Function Group and Subset Group
		try
		{
			if (! m_bInvertSubsetSelection)
				vUserData[i].ssGrp = m_spApproxSpace->subset_grp_by_name
					(vUserData[i].ssName.c_str());
			else
				vUserData[i].ssGrp = m_spApproxSpace->all_subsets_grp_except_for
					(vUserData[i].ssName.c_str());
		}
		UG_CATCH_THROW(" Subsets '"<<vUserData[i].ssName<<"' not"
		                " all contained in ApproximationSpace.");

		FunctionGroup fctGrp;
		try{
			fctGrp = m_spApproxSpace->fct_grp_by_name(vUserData[i].fctName.c_str());
		}UG_CATCH_THROW(" Functions '"<<vUserData[i].fctName<<"' not"
		                " all contained in ApproximationSpace.");

	//	check functions and subsets
		check_functions_and_subsets(fctGrp, vUserData[i].ssGrp, TUserData::numFct);

	//	set functions
		if(fctGrp.size() != TUserData::numFct)
			UG_THROW("LagrangeDirichletBoundary: wrong number of fct");

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
			mvUserDataBndSegment[subsetIndex].push_back(&vUserData[i]);
		}
	}
}


template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
extract_data()
{
//	check that function pattern exists
	if(!m_spApproxSpace.valid())
		UG_THROW("DirichletBoundary:extract_data: "
				" Approximation Space not set.");

	extract_data(m_mNumberBndSegment, m_vNumberData);
	extract_data(m_mBNDNumberBndSegment, m_vBNDNumberData);
	extract_data(m_mConstNumberBndSegment, m_vConstNumberData);
	extract_data(m_mVectorBndSegment, m_vVectorData);
	extract_data(m_mOldNumberBndSegment, m_vOldNumberData);
}

////////////////////////////////////////////////////////////////////////////////
//	assemble_dirichlet_rows
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
assemble_dirichlet_rows(matrix_type& mat, ConstSmartPtr<DoFDistribution> dd, number time)
{
	extract_data();

//	loop boundary subsets
	typename std::map<int, std::vector<CondNumberData*> >::const_iterator iter;
	for(iter = m_mBNDNumberBndSegment.begin(); iter != m_mBNDNumberBndSegment.end(); ++iter)
	{
		int si = (*iter).first;
		const std::vector<CondNumberData*>& userData = (*iter).second;

		DoFDistribution::traits<Vertex>::const_iterator iterBegin 	= dd->begin<Vertex>(si);
		DoFDistribution::traits<Vertex>::const_iterator iterEnd 	= dd->end<Vertex>(si);

	//	create Multiindex
		std::vector<DoFIndex> multInd;

	//	for readin
		MathVector<1> val;
		position_type corner;

	//	loop vertices
		for(DoFDistribution::traits<Vertex>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
		{
		//	get vertex
			Vertex* vertex = *iter;

		//	get corner position
			corner = m_aaPos[vertex];

		//	loop dirichlet functions on this segment
			for(size_t i = 0; i < userData.size(); ++i)
			{
			// 	check if function is dirichlet
				if(!(*userData[i])(val, corner, time, si)) continue;

			//	get function index
				const size_t fct = userData[i]->fct[0];

			//	get multi indices
				if(dd->inner_dof_indices(vertex, fct, multInd) != 1)
					return;

				this->m_spAssTuner->set_dirichlet_row(mat, multInd[0]);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust TRANSFER
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_prolongation(matrix_type& P,
                    ConstSmartPtr<DoFDistribution> ddFine,
                    ConstSmartPtr<DoFDistribution> ddCoarse,
					int type,
                    number time)
{
#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
	 if (!m_bAdjustTransfers)
	 {
		 std::cerr << "Avoiding  adjust_prolongation" << std::endl;
		 return;
	}
#endif
	extract_data();

	adjust_prolongation<CondNumberData>(m_mBNDNumberBndSegment, P, ddFine, ddCoarse, time);
	adjust_prolongation<NumberData>(m_mNumberBndSegment, P, ddFine, ddCoarse, time);
	adjust_prolongation<ConstNumberData>(m_mConstNumberBndSegment, P, ddFine, ddCoarse, time);

	adjust_prolongation<VectorData>(m_mVectorBndSegment, P, ddFine, ddCoarse, time);

	adjust_prolongation<OldNumberData>(m_mOldNumberBndSegment, P, ddFine, ddCoarse, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_prolongation(const std::map<int, std::vector<TUserData*> >& mvUserData,
                    matrix_type& P,
                    ConstSmartPtr<DoFDistribution> ddFine,
                    ConstSmartPtr<DoFDistribution> ddCoarse,
                    number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(ddFine->max_dofs(VERTEX)) adjust_prolongation<RegularVertex, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		if(ddFine->max_dofs(EDGE))   adjust_prolongation<Edge, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		if(ddFine->max_dofs(FACE))   adjust_prolongation<Face, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		if(ddFine->max_dofs(VOLUME)) adjust_prolongation<Volume, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_prolongation:"
						" While calling 'adapt_prolongation' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_prolongation(const std::vector<TUserData*>& vUserData, int si,
                    matrix_type& P,
                    ConstSmartPtr<DoFDistribution> ddFine,
                    ConstSmartPtr<DoFDistribution> ddCoarse,
                    number time)
{
//	create Multiindex
	std::vector<DoFIndex> vFineDoF, vCoarseDoF;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = ddFine->begin<TBaseElem>(si);
	iterEnd = ddFine->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;
		GridObject* parent = m_spDomain->grid()->get_parent(elem);
		if(!parent) continue;
		if(!ddCoarse->is_contained(parent)) continue;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = ddFine->local_finite_element_id(fct);

			//	get multi indices
				ddFine->inner_dof_indices(elem, fct, vFineDoF);
				ddCoarse->inner_dof_indices(parent, fct, vCoarseDoF);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
					UG_ASSERT(vFineDoF.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < vFineDoF.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					SetRow(P, vFineDoF[j], 0.0);
				}

				if(vFineDoF.size() > 0){
					for(size_t k = 0; k < vCoarseDoF.size(); ++k){
						DoFRef(P, vFineDoF[0], vCoarseDoF[k]) = 1.0;
					}
				}
			}
		}
	}
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_restriction(matrix_type& R,
					ConstSmartPtr<DoFDistribution> ddCoarse,
					ConstSmartPtr<DoFDistribution> ddFine,
					int type,
					number time)
{
#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
	if (!m_bAdjustTransfers)
	{
		std::cerr << "Avoiding adjust_restriction" << std::endl;
		return;
	}
#endif
	extract_data();

	adjust_restriction<CondNumberData>(m_mBNDNumberBndSegment, R, ddCoarse, ddFine, time);
	adjust_restriction<NumberData>(m_mNumberBndSegment, R, ddCoarse, ddFine, time);
	adjust_restriction<ConstNumberData>(m_mConstNumberBndSegment, R, ddCoarse, ddFine, time);

	adjust_restriction<VectorData>(m_mVectorBndSegment, R, ddCoarse, ddFine, time);

	adjust_restriction<OldNumberData>(m_mOldNumberBndSegment, R, ddCoarse, ddFine, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_restriction(const std::map<int, std::vector<TUserData*> >& mvUserData,
                   matrix_type& R,
                   ConstSmartPtr<DoFDistribution> ddCoarse,
                   ConstSmartPtr<DoFDistribution> ddFine,
                   number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(ddFine->max_dofs(VERTEX)) adjust_restriction<RegularVertex, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		if(ddFine->max_dofs(EDGE))   adjust_restriction<Edge, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		if(ddFine->max_dofs(FACE))   adjust_restriction<Face, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		if(ddFine->max_dofs(VOLUME)) adjust_restriction<Volume, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_restriction:"
						" While calling 'adjust_restriction' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_restriction(const std::vector<TUserData*>& vUserData, int si,
                   matrix_type& R,
                   ConstSmartPtr<DoFDistribution> ddCoarse,
                   ConstSmartPtr<DoFDistribution> ddFine,
                   number time)
{
//	create Multiindex
	std::vector<DoFIndex> vFineDoF, vCoarseDoF;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = ddFine->begin<TBaseElem>(si);
	iterEnd = ddFine->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;
		GridObject* parent = m_spDomain->grid()->get_parent(elem);
		if(!parent) continue;
		if(!ddCoarse->is_contained(parent)) continue;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = ddFine->local_finite_element_id(fct);

			//	get multi indices
				ddFine->inner_dof_indices(elem, fct, vFineDoF);
				ddCoarse->inner_dof_indices(parent, fct, vCoarseDoF);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, parent, *m_spDomain, lfeID);
					UG_ASSERT(vCoarseDoF.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < vCoarseDoF.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					SetRow(R, vCoarseDoF[j], 0.0);
				}

				if(vFineDoF.size() > 0){
					for(size_t k = 0; k < vCoarseDoF.size(); ++k){
						DoFRef(R, vCoarseDoF[k], vFineDoF[0]) = 1.0;
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust JACOBIAN
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
		ConstSmartPtr<DoFDistribution> dd, int type, number time,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		const number s_a0)
{
	extract_data();

	adjust_jacobian<CondNumberData>(m_mBNDNumberBndSegment, J, u, dd, time);
	adjust_jacobian<NumberData>(m_mNumberBndSegment, J, u, dd, time);
	adjust_jacobian<ConstNumberData>(m_mConstNumberBndSegment, J, u, dd, time);

	adjust_jacobian<VectorData>(m_mVectorBndSegment, J, u, dd, time);

	adjust_jacobian<OldNumberData>(m_mOldNumberBndSegment, J, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(const std::map<int, std::vector<TUserData*> >& mvUserData,
                matrix_type& J, const vector_type& u,
           	    ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_jacobian<RegularVertex, TUserData>(vUserData, si, J, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_jacobian<Edge, TUserData>(vUserData, si, J, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_jacobian<Face, TUserData>(vUserData, si, J, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_jacobian<Volume, TUserData>(vUserData, si, J, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_jacobian:"
						" While calling 'adapt_jacobian' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(const std::vector<TUserData*>& vUserData, int si,
                matrix_type& J, const vector_type& u,
           	    ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

// 	save all dirichlet degree of freedom indices.
	std::set<size_t> dirichletDoFIndices;


//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					this->m_spAssTuner->set_dirichlet_row(J, multInd[j]);
					if(m_bDirichletColumns)
						dirichletDoFIndices.insert(multInd[j][0]);
				}
			}
		}
	}


	if(m_bDirichletColumns){
	//	UG_LOG("adjust jacobian\n")

		// number of rows
		size_t nr = J.num_rows();

		// run over all rows of the local matrix J and save the colums
		// entries for the Dirichlet indices in the map

		typename std::set<size_t>::iterator currentDIndex;

		for(size_t i = 0; i<nr; i++)
		{
			for(typename matrix_type::row_iterator it = J.begin_row(i); it!=J.end_row(i); ++it){

				// look if the current index is a dirichlet index
				// if it.index is a dirichlet index
				// the iterator currentDIndex is delivered otherwise set::end()
				currentDIndex = dirichletDoFIndices.find(it.index());

				// fill dirichletMap & set corresponding entry to zero
				if(currentDIndex != dirichletDoFIndices.end()){
					// the dirichlet-dof-index it.index is assigned
					// the row i and the matrix entry it.value().
					// if necessary for defect remove comment

						//	m_dirichletMap[it.index()][i] = it.value();

					// the corresponding entry at column it.index() is set zero
					// this corresponds to a dirichlet column.
					// diagonal stays unchanged
					if(i!=it.index())
						it.value() = 0.0;
				}

			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
//	adjust DEFECT
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, int type, number time,
              ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			  const std::vector<number>* vScaleMass,
			  const std::vector<number>* vScaleStiff)
{
	extract_data();

	adjust_defect<CondNumberData>(m_mBNDNumberBndSegment, d, u, dd, time);
	adjust_defect<NumberData>(m_mNumberBndSegment, d, u, dd, time);
	adjust_defect<ConstNumberData>(m_mConstNumberBndSegment, d, u, dd, time);

	adjust_defect<VectorData>(m_mVectorBndSegment, d, u, dd, time);

	adjust_defect<OldNumberData>(m_mOldNumberBndSegment, d, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(const std::map<int, std::vector<TUserData*> >& mvUserData,
               vector_type& d, const vector_type& u,
               ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_defect<RegularVertex, TUserData>(vUserData, si, d, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_defect<Edge, TUserData>(vUserData, si, d, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_defect<Face, TUserData>(vUserData, si, d, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_defect<Volume, TUserData>(vUserData, si, d, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_defect:"
						" While calling 'adjust_defect' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(const std::vector<TUserData*>& vUserData, int si,
              vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch. (multInd.size()="<<
					          multInd.size()<<", vPos.size()="<<vPos.size()<<")");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					//	set zero for dirichlet values
					this->m_spAssTuner->set_dirichlet_val(d, multInd[j], 0.0);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust SOLUTION
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	extract_data();

	adjust_solution<CondNumberData>(m_mBNDNumberBndSegment, u, dd, time);
	adjust_solution<NumberData>(m_mNumberBndSegment, u, dd, time);
	adjust_solution<ConstNumberData>(m_mConstNumberBndSegment, u, dd, time);

	adjust_solution<VectorData>(m_mVectorBndSegment, u, dd, time);

	adjust_solution<OldNumberData>(m_mOldNumberBndSegment, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(const std::map<int, std::vector<TUserData*> >& mvUserData,
                vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_solution<RegularVertex, TUserData>(vUserData, si, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_solution<Edge, TUserData>(vUserData, si, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_solution<Face, TUserData>(vUserData, si, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_solution<Volume, TUserData>(vUserData, si, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_solution:"
						" While calling 'adjust_solution' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(const std::vector<TUserData*>& vUserData, int si,
                vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time)
{
//	check if the solution is to be adjusted
	if (! TUserData::setSolValue)
		return;
	
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	value readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				//  get dirichlet value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					this->m_spAssTuner->set_dirichlet_val(u, multInd[j], val[f]);
				}
			}
		}
	}
}



////////////////////////////////////////////////////////////////////////////////
//	adjust CORRECTION
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::adjust_correction
(
	vector_type& c,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	extract_data();

	adjust_correction<CondNumberData>(m_mBNDNumberBndSegment, c, dd, time);
	adjust_correction<NumberData>(m_mNumberBndSegment, c, dd, time);
	adjust_correction<ConstNumberData>(m_mConstNumberBndSegment, c, dd, time);

	adjust_correction<VectorData>(m_mVectorBndSegment, c, dd, time);

	adjust_correction<OldNumberData>(m_mOldNumberBndSegment, c, dd, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_correction(const std::map<int, std::vector<TUserData*> >& mvUserData,
                vector_type& c, ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt correction for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_correction<RegularVertex, TUserData>(vUserData, si, c, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_correction<Edge, TUserData>(vUserData, si,c, dd, time);
		if(dd->max_dofs(FACE))
			adjust_correction<Face, TUserData>(vUserData, si, c, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_correction<Volume, TUserData>(vUserData, si, c, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_correction:"
						" While calling 'adjust_correction' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_correction(const std::vector<TUserData*>& vUserData, int si,
                vector_type& c, ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	value readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				//  find out whether to use dirichlet value; concrete value is of no consequence
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					this->m_spAssTuner->set_dirichlet_val(c, multInd[j], 0.0);
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	adjust LINEAR
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(matrix_type& A, vector_type& b,
              ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	extract_data();

	adjust_linear<CondNumberData>(m_mBNDNumberBndSegment, A, b, dd, time);
	adjust_linear<NumberData>(m_mNumberBndSegment, A, b, dd, time);
	adjust_linear<ConstNumberData>(m_mConstNumberBndSegment, A, b, dd, time);

	adjust_linear<VectorData>(m_mVectorBndSegment, A, b, dd, time);

	adjust_linear<OldNumberData>(m_mOldNumberBndSegment, A, b, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(const std::map<int, std::vector<TUserData*> >& mvUserData,
              matrix_type& A, vector_type& b,
           	  ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_linear<RegularVertex, TUserData>(vUserData, si, A, b, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_linear<Edge, TUserData>(vUserData, si, A, b, dd, time);
		if(dd->max_dofs(FACE))
			adjust_linear<Face, TUserData>(vUserData, si, A, b, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_linear<Volume, TUserData>(vUserData, si, A, b, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_linear:"
						" While calling 'adjust_linear' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(const std::vector<TUserData*>& vUserData, int si,
              matrix_type& A, vector_type& b,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	readin value
	typename TUserData::value_type val;

// 	save all dirichlet degree of freedom indices.
	std::set<size_t> dirichletDoFIndices;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(),
						  "Mismatch: numInd="<<multInd.size()<<", numPos="
						  <<vPos.size()<<" on "<<elem->reference_object_id());

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					this->m_spAssTuner->set_dirichlet_row(A, multInd[j]);

					if(m_bDirichletColumns)
					{
						// FIXME: Beware, this is dangerous!
						//        It will not work for blocked algebras.
						UG_COND_THROW(multInd[j][1] != 0,
							"adjust_linear() is not implemented for block matrices and the symmetric case!");
						dirichletDoFIndices.insert(multInd[j][0]);
					}

					if (TUserData::setSolValue)
						this->m_spAssTuner->set_dirichlet_val(b, multInd[j], val[f]);
				}
			}
		}
	}



	if(m_bDirichletColumns){
//		UG_LOG("adjust linear\n")
		m_A = &A;
		// number of rows
		size_t nr = A.num_rows();

		typename std::set<size_t>::iterator currentDIndex;

		// run over all rows of the local matrix J and save the column
		// entries for the Dirichlet indices in the map

		for(size_t i = 0; i<nr; i++)
		{
			// do not save any entries in Dirichlet rows!
			if (dirichletDoFIndices.find(i) != dirichletDoFIndices.end())
				continue;

			for(typename matrix_type::row_iterator it = A.begin_row(i); it!=A.end_row(i); ++it)
			{
				currentDIndex = dirichletDoFIndices.find(it.index());

				// fill dirichletMap & set corresponding entry to zero
				if(currentDIndex != dirichletDoFIndices.end()){

					// the dirichlet-dof-index it.index is assigned
					// the row i and the matrix entry it.value().
					m_dirichletMap[si][i][it.index()] = it.value();
					it.value() = 0.0;
				}
			}
		}

		// TODO: for better efficiency use vectors instead of maps (rows and columns are ordered!)
		typename std::map<int, std::map<int, value_type> >::iterator itdirichletMap;
		typename std::map<int, value_type>::iterator itdirichletRowMap;
		for(size_t i = 0; i<nr; i++)
		{
			// step over if this row did not contain any connections to Dirichlet nodes
			if ((itdirichletMap = m_dirichletMap[si].find(i)) == m_dirichletMap[si].end())
				continue;

			for(typename matrix_type::row_iterator it = m_A->begin_row(i); it!=m_A->end_row(i); ++it){

				// current column index is a dirichlet index
				if ((itdirichletRowMap = itdirichletMap->second.find(it.index())) != itdirichletMap->second.end())
					b[i] -= itdirichletRowMap->second*b[it.index()];
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust RHS
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(vector_type& b, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	extract_data();

	adjust_rhs<CondNumberData>(m_mBNDNumberBndSegment, b, u, dd, time);
	adjust_rhs<NumberData>(m_mNumberBndSegment, b, u, dd, time);
	adjust_rhs<ConstNumberData>(m_mConstNumberBndSegment, b, u, dd, time);

	adjust_rhs<VectorData>(m_mVectorBndSegment, b, u, dd, time);

	adjust_rhs<OldNumberData>(m_mOldNumberBndSegment, b, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(const std::map<int, std::vector<TUserData*> >& mvUserData,
           vector_type& b, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_rhs<RegularVertex, TUserData>(vUserData, si, b, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_rhs<Edge, TUserData>(vUserData, si, b, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_rhs<Face, TUserData>(vUserData, si, b, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_rhs<Volume, TUserData>(vUserData, si, b, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_rhs:"
						" While calling 'adjust_rhs' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(const std::vector<TUserData*>& vUserData, int si,
           vector_type& b, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					if (TUserData::setSolValue)
						this->m_spAssTuner->set_dirichlet_val(b, multInd[j], val[f]);
					else
						this->m_spAssTuner->set_dirichlet_val(b, multInd[j], DoFRef(u, multInd[j]));
				}
			}
		}

	}


	// adjust the right hand side
	if(m_bDirichletColumns){
		typename std::map<int, std::map<int, value_type> >::iterator itdirichletMap;
		typename std::map<int, value_type>::iterator itdirichletRowMap;
		size_t nr = m_A->num_rows();
		for(size_t i = 0; i<nr; i++)
		{
			// step over if this row did not contain any connections to Dirichlet nodes
			if ((itdirichletMap = m_dirichletMap[si].find(i)) == m_dirichletMap[si].end())
				continue;

			for(typename matrix_type::row_iterator it = m_A->begin_row(i); it!=m_A->end_row(i); ++it){

				// current column index is a dirichlet index
				if ((itdirichletRowMap = itdirichletMap->second.find(it.index())) != itdirichletMap->second.end())
					b[i] -= itdirichletRowMap->second*b[it.index()];
			}
		}
	}
}


// //////////////////////////////////////////////////////////////////////////////
//	adjust error
// //////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_error
(
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const std::vector<number>* vScaleMass,
	const std::vector<number>* vScaleStiff
)
{
	//	get the error estimator data object and check that it is of the right type
	if (this->m_spErrEstData.get() == NULL)
	{
		UG_THROW("No ErrEstData object has been given to this constraint!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to MultipleSideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}


	extract_data();

	adjust_error<CondNumberData>(m_mBNDNumberBndSegment, u, dd, time);
	adjust_error<NumberData>(m_mNumberBndSegment, u, dd, time);
	adjust_error<ConstNumberData>(m_mConstNumberBndSegment, u, dd, time);
	
	adjust_error<VectorData>(m_mVectorBndSegment, u, dd, time);
	
	adjust_error<OldNumberData>(m_mOldNumberBndSegment, u, dd, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_error
(
	const std::map<int, std::vector<TUserData*> >& mvUserData,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
	// cast error estimator data object to the right type
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	typedef typename err_est_type::side_type side_type;

	// loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for (iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
		// get subset index
		const int si = (*iter).first;

		// get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

		try
		{
			// create multi-index
			std::vector<DoFIndex> multInd;

			// dummy for readin
			typename TUserData::value_type val;

			// position of dofs
			std::vector<position_type> vPos;

			// iterators
			typename DoFDistribution::traits<side_type>::const_iterator iter, iterEnd;
			iter = dd->begin<side_type>(si);
			iterEnd = dd->end<side_type>(si);

			// loop elements of side_type (only!)
			// elements of measure 0 in the boundary are ignored.
			for( ; iter != iterEnd; iter++)
			{
				// get vertex
				side_type* elem = *iter;

				// get reference object id
				ReferenceObjectID roid = elem->reference_object_id();

				// get corner coords (for later use in calculating global IPs)
				std::vector<typename TDomain::position_type> vCoCo;
				CollectCornerCoordinates(vCoCo, elem, *m_spDomain, false);

				// loop dirichlet functions on this segment
				for(size_t i = 0; i < vUserData.size(); ++i)
				{
					for(size_t f = 0; f < TUserData::numFct; ++f)
					{
						// get function index
						const size_t fct = vUserData[i]->fct[f];

						// get lfeID for function
						LFEID lfeID = dd->local_finite_element_id(fct);

						// get local and global IPs
						size_t numSideIPs;
						const MathVector<side_type::dim>* sideLocIPs;
						const MathVector<dim>* sideGlobIPs;

						try
						{
							numSideIPs = err_est_data->get(fct)->num_side_ips(roid);
							sideLocIPs = err_est_data->get(fct)->template side_local_ips<side_type::dim>(roid);
							sideGlobIPs = err_est_data->get(fct)->side_global_ips(elem, &vCoCo[0]);
						}
						UG_CATCH_THROW("Global integration points for error estimator cannot be determined.");

						// loop IPs
						for (size_t ip = 0; ip < numSideIPs; ++ip)
						{
							// get Dirichlet value (and do nothing, if conditional D value is false)
							if (!(*vUserData[i])(val, sideGlobIPs[ip], time, si)) continue;
							
							// check if we take the values directly from the solution
							if (! TUserData::setSolValue)
							{
								(*err_est_data->get(fct))(elem,ip) = 0;
								continue;
							}

							// evaluate shapes at ip
							LFEID new_lfeID(lfeID.type(), lfeID.dim()-1, lfeID.order());
							const LocalShapeFunctionSet<side_type::dim>& rTrialSpace =
								LocalFiniteElementProvider::get<side_type::dim>(roid, new_lfeID);
							std::vector<number> vShape;
							rTrialSpace.shapes(vShape, sideLocIPs[ip]);

							// get multiindices of element
							std::vector<DoFIndex> ind;
							dd->dof_indices(elem, fct, ind);

							// compute solution at integration point
							number uAtIP = 0.0;
							for (size_t sh = 0; sh < vShape.size(); ++sh)
								uAtIP += DoFRef(u, ind[sh]) * vShape[sh];

							// set error integrand value at IP
							(*err_est_data->get(fct))(elem,ip) = val[f] - uAtIP;
						}
					}
				}
			}
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_error:"
						" While calling 'adjust_error' for TUserData, aborting.");
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__ */
