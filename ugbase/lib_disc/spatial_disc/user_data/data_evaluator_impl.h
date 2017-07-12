/*
 * Copyright (c) 2011-2017:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Arne Naegel
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__

namespace ug{

template <typename TDomain, typename TElemDisc>
DataEvaluatorBase<TDomain, TElemDisc>::
DataEvaluatorBase(int discPart,
              const std::vector<TElemDisc*>& vElemDisc,
              ConstSmartPtr<FunctionPattern> fctPat,
              const bool bNonRegularGrid,
              LocalVectorTimeSeries* pLocTimeSeries,
              const std::vector<number>* pvScaleMass,
              const std::vector<number>* pvScaleStiff)
   : m_spFctPattern(fctPat)
{
// 	remember infos
	m_discPart = discPart;
	m_pLocTimeSeries = pLocTimeSeries;
	m_bNeedLocTimeSeries = false; // initially
	m_bUseHanging = false; // initially

	m_vElemDisc[PT_ALL] = vElemDisc;

//	create FunctionIndexMapping for each Disc
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get elem disc
		TElemDisc* disc = m_vElemDisc[PT_ALL][i];

	// 	handle time dependency
		if(pLocTimeSeries != NULL && pvScaleMass != NULL && pvScaleStiff != NULL){
			disc->set_time_dependent(*pLocTimeSeries, *pvScaleMass, *pvScaleStiff);
		}
		else if(pLocTimeSeries != NULL){
			disc->set_time_dependent(*pLocTimeSeries, std::vector<number>(), std::vector<number>());
		}
		else{
			disc->set_time_independent();
		}

	// 	checks
		disc->check_setup(bNonRegularGrid);

	//	cache time dependency
		m_bNeedLocTimeSeries |= disc->local_time_series_needed();

	//	cache grid type required
		if(bNonRegularGrid)
			m_bUseHanging |= disc->use_hanging();

	//	sort by process type
		ProcessType process;
		if(!disc->is_time_dependent()) process = PT_STATIONARY;
		else process = PT_INSTATIONARY;
		m_vElemDisc[process].push_back(disc);
	}
}


///////////////////////////////////////////////////////////////////////////////
// DataEvaluatorBase Setup
///////////////////////////////////////////////////////////////////////////////



template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::
add_data_to_eval_data(std::vector<SmartPtr<ICplUserData<dim> > >& vEvalData,
                      std::vector<SmartPtr<ICplUserData<dim> > >& vTryingToAdd)
{
//	if empty, we're done
	if(vTryingToAdd.empty()) return;

//	search for element in already scheduled data
	typename std::vector<SmartPtr<ICplUserData<dim> > >::iterator it, itEnd;
	it = find(vEvalData.begin(), vEvalData.end(), vTryingToAdd.back());

//	if found, skip this data
	if(it != vEvalData.end())
	{
		vTryingToAdd.pop_back();
		return;
	}

//	search if element already contained in list. Then, the element
//	did start the adding procedure before and a circle dependency
//	is found
	itEnd = vTryingToAdd.end(); itEnd--;
	it = find(vTryingToAdd.begin(), itEnd, *itEnd);

//	if found, return error of circle dependency
	if(it != itEnd)
		UG_THROW("DataEvaluatorBase::add_data_to_eval_data:"
						" Circle dependency of data detected for UserData.");

//	add all dependent datas
	SmartPtr<ICplUserData<dim> > data = vTryingToAdd.back();
	for(size_t i = 0; i < data->num_needed_data(); ++i)
	{
	//	add each data separately
		vTryingToAdd.push_back(data->needed_data(i));
		add_data_to_eval_data(vEvalData, vTryingToAdd);
	}

//	add this data to the evaluation list
	vEvalData.push_back(data);

//	pop last one, since now added to eval list
	if(!vTryingToAdd.empty())
		vTryingToAdd.pop_back();
}

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::extract_imports_and_userdata(int subsetIndex, int discPart)
{
	clear_extracted_data_and_mappings();

//	queue for all user data needed
	std::vector<SmartPtr<ICplUserData<dim> > > vEvalData;
	std::vector<SmartPtr<ICplUserData<dim> > > vTryingToAdd;

//	In the next loop we extract all needed UserData:
//	We only process the DataImport if there has been set data to the import
//	since otherwise no evaluation is needed.
//	If there is data given, we get the connected UserData and add it to the vector
//	of EvaluationData. This simply adds the UserData to the queue for UserData, if
//	the data does not depend on other Data. But if the UserData itself has
//	dependencies to other UserData, this data is added first (in a recursive
//	process). Of course, no circle dependency between UserData is allowed.

//	In the same loop over the data imports, we schedule the DataImports for
//	evaluation and compute the correct FunctionMapping for the linearization
//	of the defect and the Data, the Import is connected to:
//	If the UserData does not depend on the primary unknowns, we're done. Else
//	we have to setup the Function mappings between the common function group
//	and the DataImport-FunctionGroup. This is simply the same function map as
//	for the element discretization, since the DataImport depends by definition
//	from and only from the primary variables of its associated IElemDisc.

//	loop elem discs
	for(size_t d = 0; d < m_vElemDisc[PT_ALL].size(); ++d)
	{
		TElemDisc* disc = m_vElemDisc[PT_ALL][d];

	//	loop imports
		for(size_t i = 0; i < disc->num_imports(); ++i)
		{
		//	get import
			IDataImport<dim>* iimp = &(disc->get_import(i));

		//	skip non-given data (no need for evaluation)
			if(!iimp->data_given()) continue;

		//	check part
			if( !(iimp->part() & discPart) ) continue;

		//	check correct process type
			if(iimp->part() == MASS)
				if(!disc->is_time_dependent()) continue;

		//	push export on stack of needed data
			vTryingToAdd.push_back(iimp->data());

		//	add data and all dependency to evaluation list
			try{
				add_data_to_eval_data(vEvalData, vTryingToAdd);
			}
			UG_CATCH_THROW("DataEvaluatorBase:"
						" Circle dependency of data detected for UserData.");

		//	check that queue is empty now, else some internal error occured
			if(!vTryingToAdd.empty())
				UG_THROW("DataEvaluatorBase:"
						" Internal Error, UserData queue not empty after adding.");

		//	done if and only if zero-derivative
			if(iimp->zero_derivative()) continue;

		//	remember Import
			ProcessType process;
			if(!disc->is_time_dependent()) process = PT_STATIONARY;
			else process = PT_INSTATIONARY;

			m_vImport[PT_ALL][iimp->part()].push_back(iimp);
			m_vImport[process][iimp->part()].push_back(iimp);
		}
	}

//	Now, we have processed all imports, that must be evaluated and have a long
//	vector of UserData that is connected to those imports. The UserData is already
//	sorted in this way: Data that depends on other data appears after the data
//	it depends on. This is important since we will schedule now the data for
//	evaluation and the data, that is needed by other data, will be computed
//	first. In addition, the data linker have to update their FunctionGroup and
//	must be sure that the data they depend on has already a correct FunctionGroup
//	set. This all is ensured by the (already produced) correct ordering.
//
//	In the next loop we process all UserData, that will be evaluated during
//	assembling (i.e. is connected to an Import). First, we check if the data
//	is constant. If so simply add it the the Constant Data vector; nothing more
//	has to be done here. Else we check if the data depends on the primary
//	unknowns. If this is not the case, the UserData must be a Position-dependent
//	data, but not constant. Thus, schedule it at the Position Data vector.
//	If the data depends on the primary unknowns we must proceed as follows:
//	First, we update the FunctionGroup of the Data, since it could be a linker
//	and having an incorrect FunctionGroup (iff the FunctionGroup of the data
//	the linker depends on has been changed). Then we create the function
//	mapping between the functions the linker depends on and the common Function
//	Group.

//	loop all needed user data and group it
	for(size_t i = 0; i < vEvalData.size(); ++i)
	{
	//	get the user data
		SmartPtr<ICplUserData<dim> > ipData = vEvalData[i];

	//	update function pattern (this will update functionGroups and Map of Data)
		try{
			ipData->set_function_pattern(m_spFctPattern);
		}
		UG_CATCH_THROW("DataEvaluatorBase: Cannot set FunctionPattern to UserData.");

	//	sort data into const and non-solution dependent
		if(ipData->constant()) {m_vConstData.push_back(ipData); continue;}
		if(ipData->zero_derivative()){m_vPosData.push_back(ipData); continue;}

	//	save as dependent data
		m_vDependentData.push_back(ipData);
	}

// 	Handle time dependency
// 	NOTE: constant data is not processed.
	if(m_pLocTimeSeries != NULL){
		for(size_t i = 0; i < m_vPosData.size(); ++i)
			m_vPosData[i]->set_times(m_pLocTimeSeries->times());
		for(size_t i = 0; i < m_vDependentData.size(); ++i)
			m_vDependentData[i]->set_times(m_pLocTimeSeries->times());
	}

// 	NOTE: constant data is not processed, since constant == independent of si
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->set_subset(subsetIndex);
	for(size_t i = 0; i < m_vDependentData.size(); ++i)
		m_vDependentData[i]->set_subset(subsetIndex);
}

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::clear_extracted_data_and_mappings()
{
	for(int type = 0; type < MAX_PROCESS; ++type){
		m_vImport[type][MASS].clear();
		m_vImport[type][STIFF].clear();
		m_vImport[type][RHS].clear();
	}

	m_vConstData.clear();
	m_vPosData.clear();
	m_vDependentData.clear();
}



template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::clear_positions_in_user_data()
{
//	remove ip series for all used UserData
	for(size_t i = 0; i < m_vConstData.size(); ++i) m_vConstData[i]->clear();
	for(size_t i = 0; i < m_vPosData.size(); ++i)   m_vPosData[i]->clear();
	for(size_t i = 0; i < m_vDependentData.size(); ++i) m_vDependentData[i]->clear();
}


///////////////////////////////////////////////////////////////////////////////
// prepare / finish
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::set_time_point(const size_t timePoint)
{
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
		m_vElemDisc[PT_ALL][i]->set_time_point(timePoint);

	// NOTE: constant data is not processed.
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->set_time_point(timePoint);
	for(size_t i = 0; i < m_vDependentData.size(); ++i)
		m_vDependentData[i]->set_time_point(timePoint);
}



///////////////////////////////////////////////////////////////////////////////
// Error estimator's routines
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::
prepare_err_est_elem_loop(const ReferenceObjectID id, int si)
{
// 	prepare loop (elem disc set local ip series here)
	try{
		for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_prep_err_est_elem_loop(id, si);
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_err_est_elem_loop: "
						"Cannot prepare element loop.");

//	extract data imports and user data
	try{
		extract_imports_and_userdata(si, m_discPart);
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_err_est_elem_loop: "
					"Cannot extract imports and userdata.");

//	check setup of imports
	try{
		for(size_t i = 0; i < m_vImport[PT_ALL][MASS].size(); ++i)
			m_vImport[PT_ALL][MASS][i]->check_setup();
		for(size_t i = 0; i < m_vImport[PT_ALL][STIFF].size(); ++i)
			m_vImport[PT_ALL][STIFF][i]->check_setup();
		for(size_t i = 0; i < m_vImport[PT_ALL][RHS].size(); ++i)
			m_vImport[PT_ALL][RHS][i]->check_setup();
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_err_est_elem_loop: Import not correctly implemented.");

//	prepare and check dependent data
	try{
		for(size_t i = 0; i < m_vDependentData.size(); ++i){
			m_vDependentData[i]->check_setup();
		}
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_err_est_elem_loop: Dependent UserData "
				   " (e.g. Linker or Export) is not ready for evaluation.");

//	evaluate constant data
	for(size_t i = 0; i < m_vConstData.size(); ++i)
		m_vConstData[i]->compute((LocalVector*)NULL, NULL, NULL, false);
}

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::finish_err_est_elem_loop()
{
//	finish each elem error estimator disc
	try{
		for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_fsh_err_est_elem_loop();
	}
	UG_CATCH_THROW("DataEvaluatorBase::finish_err_est_elem_loop: Cannot finish error estimator element loop");

//	clear positions at user data
	clear_positions_in_user_data();
}

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::prepare_err_est_elem
(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[],
const LocalIndices& ind, bool bDeriv)
{
	try
	{
		for (size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_prep_err_est_elem(u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::compute_elem_err_est: Cannot prepare element.");

//	evaluate position data
	for (size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->compute(&u, elem, vCornerCoords, false);

// 	process dependent data:
//	We can not simply compute exports first, then Linker, because an export
//	itself could depend on other data if implemented somehow in the IElemDisc
//	(e.g. using data from some DataImport). Thus, we have to loop the sorted
//	vector of all dependent data (that is correctly sorted the way that always
//	needed data has previously computed).

//	compute the data
	try
	{
		for (size_t i = 0; i < m_vDependentData.size(); ++i)
		{
			u.access_by_map(m_vDependentData[i]->map());
			m_vDependentData[i]->compute(&u, elem, vCornerCoords, false);
		}
	}
	UG_CATCH_THROW("DataEvaluatorBase::compute_elem_err_est: Cannot compute data for Export or Linker.");
}


template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::
compute_err_est_A_elem
(	LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[],
	const LocalIndices& ind,
	const number scaleMass,
	const number scaleStiff)
{
	UG_ASSERT(m_discPart & STIFF, "Using compute_err_est_A_elem, but not STIFF requested.");
	try
	{
		for (std::size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_compute_err_est_A_elem(u, elem, vCornerCoords, scaleStiff);
	}
	UG_CATCH_THROW("DataEvaluatorBase::compute_err_est_A_elem: Cannot assemble stiffness part of error estimator");
}

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::
compute_err_est_M_elem
(	LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[],
	const LocalIndices& ind,
	const number scaleMass,
	const number scaleStiff)
{
	UG_ASSERT(m_discPart & MASS, "Using compute_err_est_M_elem, but not MASS requested.");
	try
	{
		for (std::size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_compute_err_est_M_elem(u, elem, vCornerCoords, scaleMass);
	}
	UG_CATCH_THROW("DataEvaluatorBase::compute_err_est_A_elem: Cannot assemble stiffness part of error estimator");
}

template <typename TDomain, typename TElemDisc>
void DataEvaluatorBase<TDomain, TElemDisc>::
compute_err_est_rhs_elem
(	GridObject* elem,
	const MathVector<dim> vCornerCoords[],
	const LocalIndices& ind,
	const number scaleMass,
	const number scaleStiff)
{
	UG_ASSERT(m_discPart & RHS, "Using compute_err_est_rhs_elem, but not RHS requested.");
	try
	{
		for (std::size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_compute_err_est_rhs_elem(elem, vCornerCoords, scaleStiff);
	}
	UG_CATCH_THROW("DataEvaluatorBase::compute_err_est_rhs_elem: Cannot assemble rhs part of error estimator");
}


} // namespace ug


#endif
