/*
 * data_evaluator_impl.h
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include <sstream>

#include "data_evaluator.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// DataEvaluator Setup
///////////////////////////////////////////////////////////////////////////////

DataEvaluator::DataEvaluator(int discPart,
                             const std::vector<IElemDisc*>& vElemDisc,
                             const FunctionPattern& fctPat,
                             const int subset,
                             const bool bNonRegularGrid,
                             LocalVectorTimeSeries* pLocTimeSeries,
                             const std::vector<number>* pvScaleMass,
                             const std::vector<number>* pvScaleStiff)
{
//	remember needed disc parts
	m_discPart = discPart;

//	currently only fast assembles allowed
	for(size_t i = 0; i < vElemDisc.size(); ++i){
		if(!vElemDisc[i]->fast_ass_elem_enabled()){
			UG_THROW("DataEvaluator: currently only fast assemble allowed."
					 " Please use enable_fast_ass_elem in all IElemDisc.");
		}
	}

// handle time dependency
	if(pLocTimeSeries != NULL && pvScaleMass != NULL && pvScaleStiff != NULL){
		for(size_t i = 0; i < m_vElemDisc.size(); ++i)
			vElemDisc[i]->set_time_dependent(*pLocTimeSeries, *pvScaleMass, *pvScaleStiff);
	}
	else if(pLocTimeSeries != NULL){
		for(size_t i = 0; i < m_vElemDisc.size(); ++i)
			vElemDisc[i]->set_time_dependent(*pLocTimeSeries, std::vector<number>(), std::vector<number>());
	}
	else{
		for(size_t i = 0; i < m_vElemDisc.size(); ++i)
			vElemDisc[i]->set_time_independent();
	}

// remember time series
	m_pLocTimeSeries = pLocTimeSeries;
	m_bNeedLocTimeSeries = false; // initially
	m_bUseHanging = false; // initially
	m_subset = subset;

//	set function pattern to common function group
	m_commonFctGroup.set_function_pattern(fctPat);
	m_commonFctGroup.add_all();

	m_vElemDisc.clear();
	m_vElemDisc.resize(vElemDisc.size());

//	create FunctionIndexMapping for each Disc
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
	//	store elem disc
		m_vElemDisc[i].elemDisc = vElemDisc[i];

	//	create function group of this elem disc
		try{
			m_vElemDisc[i].fctGrp.set_function_pattern(fctPat);
			m_vElemDisc[i].fctGrp.add(m_vElemDisc[i].elemDisc->symb_fcts());
		}UG_CATCH_THROW("'DataEvaluator::set_elem_discs': Cannot find "
					"some symbolic Function Name for disc "<<i<<".");

	//	create a mapping between all functions and the function group of this
	//	element disc.
		try{
			CreateFunctionIndexMapping(m_vElemDisc[i].map, m_vElemDisc[i].fctGrp,
			                           m_commonFctGroup);
		}UG_CATCH_THROW("'DataEvaluator::set_elem_discs': Cannot create "
						"Function Index Mapping for disc "<<i<<".");

	////////////////////////
	// checks
	////////////////////////

	//	check that all functions are defined on chosen subsets
		SubsetGroup discSubsetGrp(fctPat.subset_handler(), m_vElemDisc[i].elemDisc->symb_subsets());

	//	check that all functions are defined on chosen subsets
		for(size_t fct = 0; fct < m_vElemDisc[i].fctGrp.size(); ++fct)
		{
			for(size_t si = 0; si < discSubsetGrp.size(); ++si)
			{
				if(!fctPat.is_def_in_subset(m_vElemDisc[i].fctGrp[fct], discSubsetGrp[si])){
					UG_LOG("WARNING in 'DataEvaluator::set_elem_discs': On disc "<<i<<
						   ": symbolic Function "<< m_vElemDisc[i].elemDisc->symb_fcts()[fct]
					 << " is not defined on subset "<<m_vElemDisc[i].elemDisc->symb_subsets()[si]
					 << ". This may be senseful only in particular cases.\n");
				}
			}
		}

	//	check correct number of functions
		if(m_vElemDisc[i].fctGrp.size() != m_vElemDisc[i].elemDisc->num_fct())
		{
			std::stringstream ss;
			ss << "DataEvaluator::set_elem_discs: Elem Disc "<<i<<
					" requires "<<m_vElemDisc[i].elemDisc->num_fct()<<" symbolic "
					"Function Name, but "<<m_vElemDisc[i].fctGrp.size()<<" Functions "
					" specified: ";
			for(size_t f=0; f < m_vElemDisc[i].elemDisc->symb_fcts().size(); ++f)
			{
				if(f > 0) ss << ", ";
				ss << vElemDisc[i]->symb_fcts()[f];
			}
			UG_THROW(ss.str());
		}

	//	request assembling for local finite element id
		std::vector<LFEID> vLfeID(m_vElemDisc[i].fctGrp.size());
		for(size_t f = 0; f < vLfeID.size(); ++f)
			vLfeID[f] = m_vElemDisc[i].fctGrp.local_finite_element_id(f);
		if(!(m_vElemDisc[i].elemDisc->request_finite_element_id(vLfeID)))
		{
			std::stringstream ss;
			ss << "DataEvaluator::set_elem_discs: Elem Disc "<<i<<
				" can not assemble the specified local finite element space set:";
			for(size_t f=0; f < m_vElemDisc[i].elemDisc->symb_fcts().size(); ++f)
			{
				ss << "  Fct "<<f<<": '"<<m_vElemDisc[i].elemDisc->symb_fcts()[f];
				ss << "' using "<< vLfeID[f];
			}
			UG_THROW(ss.str());
		}

	//	check if time dependent
		m_vElemDisc[i].needLocTimeSeries = m_vElemDisc[i].elemDisc->requests_local_time_series();
		m_bNeedLocTimeSeries |= m_vElemDisc[i].needLocTimeSeries;

	//  let disc use non-regular grid assemblings
		if(!m_vElemDisc[i].elemDisc->request_non_regular_grid(bNonRegularGrid))
		{
			UG_THROW("DataEvaluator::set_non_regular_grid: "
					" Elem Disc " << i << " does not support non-regular"
					" grids, but this is requested.\n");
		}

		if(bNonRegularGrid)
			m_bUseHanging |= m_vElemDisc[i].elemDisc->use_hanging();
	}

}

void DataEvaluator::clear_extracted_data_and_mappings()
{
	m_vImport[MASS].clear();
	m_vImport[STIFF].clear();
	m_vImport[RHS].clear();

	m_vConstData.clear();
	m_vPosData.clear();
	m_vDependentData.clear();
	m_vDependentMap.clear();
}

void DataEvaluator::add_data_to_eval_data(std::vector<SmartPtr<IUserData> >& vEvalData,
                                          std::vector<SmartPtr<IUserData> >& vTryingToAdd)
{
//	if empty, we're done
	if(vTryingToAdd.empty()) return;

//	search for element in already scheduled data
	std::vector<SmartPtr<IUserData> >::iterator it, itEnd;
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
		UG_THROW("DataEvaluator::add_data_to_eval_data:"
						" Circle dependency of data detected for UserData.");

//	add all dependent datas
	SmartPtr<IUserData> data = vTryingToAdd.back();
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

void DataEvaluator::extract_imports_and_userdata(int discPart)
{
	clear_extracted_data_and_mappings();

//	copy function group in import/export of element discs
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
		for(size_t imp = 0; imp < m_vElemDisc[i].elemDisc->num_imports(); ++imp)
			m_vElemDisc[i].elemDisc->get_import(imp).set_function_group(m_vElemDisc[i].fctGrp);

		for(size_t exp = 0; exp < m_vElemDisc[i].elemDisc->num_exports(); ++exp)
			m_vElemDisc[i].elemDisc->get_export(exp)->set_function_group(m_vElemDisc[i].fctGrp);
	}

//	queue for all user data needed
	std::vector<SmartPtr<IUserData> > vEvalData;
	std::vector<SmartPtr<IUserData> > vTryingToAdd;

//	In the next loop we extract all need UserData:
//	We only process the DataImport if there has been set data to the import
//	since otherwise no evaluation is needed.
//	If there is data given, we get the connected UserData and add it to the vector
//	of EvaluationData. This simply adds the UserData to the queue for UserData, if
//	the data does not depend on other Data. But if the UserData itself has
//	dependencies to other UserData, this data is added first (in a recursive
//	process). Of coarse, no circle dependency between UserData is allowed.

//	loop elem discs
	for(size_t d = 0; d < m_vElemDisc.size(); ++d)
	{
	//	check correct process type
		if(discPart & MASS)
			if(m_vElemDisc[d].elemDisc->is_stationary()) continue;

	//	loop imports
		for(size_t i = 0; i < m_vElemDisc[d].elemDisc->num_imports(); ++i)
		{
		//	get import
			IDataImport* iimp = &(m_vElemDisc[d].elemDisc->get_import(i));

		//	skip non-given data (no need for evaluation)
			if(!iimp->data_given()) continue;

		//	check part
			if( !(iimp->part() & discPart) ) continue;

		//	push export on stack of needed data
			vTryingToAdd.push_back(iimp->data());

		//	add data and all dependency to evaluation list
			try{
				add_data_to_eval_data(vEvalData, vTryingToAdd);
			}
			UG_CATCH_THROW("DataEvaluator:"
						" Circle dependency of data detected for UserData.");

		//	check that queue is empty now, else some internal error occured
			if(!vTryingToAdd.empty())
				UG_THROW("DataEvaluator:"
						" Internal Error, UserData queue not empty after adding.");
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
		SmartPtr<IUserData> ipData = vEvalData[i];

	//	sort data into const and non-solution dependent
		if(ipData->constant()) {m_vConstData.push_back(ipData); continue;}
		if(ipData->zero_derivative()){m_vPosData.push_back(ipData); continue;}

	//	update function group of dependent data
		try{ipData->update_function_group();}
		UG_CATCH_THROW("DataEvaluator: Cannot update FunctionGroup of IDependentData.");

	//	create FuncMap
		FunctionIndexMapping map;
		try{CreateFunctionIndexMapping(map, ipData->function_group(),
									   m_commonFctGroup);
		}UG_CATCH_THROW("DataEvaluator:Cannot create Function Index Mapping for IDependData.");

	//	save as dependent data
		m_vDependentData.push_back(ipData);
		m_vDependentMap.push_back(map);
	}

//	In a second loop over the data imports, we schedule the DataImports for
//	evaluation and compute the correct FunctionMapping for the linearization
//	of the defect and the Data, the Import is connected to:
//	If the UserData does not depend on the primary unknowns, we're done. Else
//	we have to setup the Function mappings between the common function group
//	and the DataImport-FunctionGroup. This is simply the same function map as
//	for the element discretization, since the DataImport depends by definition
//	from and only from the primary variables of its associated IElemDisc.

//	loop elem discs
	for(size_t d = 0; d < m_vElemDisc.size(); ++d)
	{
	//	check correct process type
		if(discPart & MASS)
			if(m_vElemDisc[d].elemDisc->is_stationary()) continue;

	//	loop imports
		for(size_t i = 0; i < m_vElemDisc[d].elemDisc->num_imports(); ++i)
		{
		//	get import
			IDataImport* iimp = &(m_vElemDisc[d].elemDisc->get_import(i));

		//	skip non-given data (no need for evaluation)
			if(!iimp->data_given()) continue;

		//	check part
			if( !(iimp->part() & discPart) ) continue;

		//	done if and only if zero-derivative
			if(iimp->zero_derivative()) continue;

		//	get and cast dependent data
			SmartPtr<IUserData> dependData = iimp->data();

		//	create FuncMap for data
		//	this is ok, since the function group has been updated in the
		//	previous loop over all needed data
			FunctionIndexMapping map;
			try{CreateFunctionIndexMapping(map, dependData->function_group(),
										   m_commonFctGroup);
			}UG_CATCH_THROW("DataEvaluator: Cannot create Function Index Mapping for DependentData.");

		//	remember Import
			m_vImport[iimp->part()].push_back(Import(iimp, m_vElemDisc[d].map, map));
		}
	}

	// Handle time dependency
	// NOTE: constant data is not processed.
	if(m_pLocTimeSeries != NULL){
		for(size_t i = 0; i < m_vPosData.size(); ++i)
			m_vPosData[i]->set_times(m_pLocTimeSeries->times());
		for(size_t i = 0; i < m_vDependentData.size(); ++i)
			m_vDependentData[i]->set_times(m_pLocTimeSeries->times());
	}

	// NOTE: constant data is not processed, since constant == independent of si
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->set_subset(m_subset);
	for(size_t i = 0; i < m_vDependentData.size(); ++i)
		m_vDependentData[i]->set_subset(m_subset);

}

void DataEvaluator::set_time_point(const size_t timePoint)
{
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
		m_vElemDisc[i].elemDisc->set_time_point(timePoint);

	// NOTE: constant data is not processed.
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->set_time_point(timePoint);
	for(size_t i = 0; i < m_vDependentData.size(); ++i)
		m_vDependentData[i]->set_time_point(timePoint);
}

///////////////////////////////////////////////////////////////////////////////
// Assemble routines
///////////////////////////////////////////////////////////////////////////////

void DataEvaluator::compute_elem_data(LocalVector& u, GeometricObject* elem, bool bDeriv)
{
//	evaluate position data
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->compute(&u, elem, false);

// 	process dependent data:
//	We can not simply compute exports first, then Linker, because an export
//	itself could depend on other data if implemented somehow in the IElemDisc
//	(e.g. using data from some DataImport). Thus, we have to loop the sorted
//	vector of all dependent data (that is correctly sorted the way that always
//	needed data has previously computed). We look up, if a Export is given, if
//	so compute it, else compute the linker

//	loop all dependent data
	for(size_t i = 0; i < m_vDependentData.size(); ++i)
	{
	//	access needed components
		u.access_by_map(m_vDependentMap[i]);

	//	compute the data
		try{
			m_vDependentData[i]->compute(&u, elem, bDeriv);
		}
		UG_CATCH_THROW("DataEvaluator::compute_elem_data:"
						"Cannot compute data for Export " << i);
	}
}

void DataEvaluator::add_JA_elem(LocalMatrix& A, LocalVector& u, GeometricObject* elem)
{
	UG_ASSERT(m_discPart & STIFF, "Using add_JA_elem, but not STIFF requested.")

	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		A.access_by_map(map(i));

		if(m_vElemDisc[i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	assemble JA
		try{
			m_vElemDisc[i].elemDisc->fast_ass_JA_elem(A, u);
		}
		UG_CATCH_THROW("DataEvaluator::ass_JA_elem: "
						"Cannot assemble Jacobian (A) for IElemDisc "<<i);
	}

	add_coupl_JA(A, u);
}

void DataEvaluator::add_JM_elem(LocalMatrix& M, LocalVector& u, GeometricObject* elem)
{
	UG_ASSERT(m_discPart & MASS, "Using add_JM_elem, but not MASS requested.")

	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		M.access_by_map(map(i));

		if(m_vElemDisc[i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	assemble JM
		try{
			if(!m_vElemDisc[i].elemDisc->is_stationary())
				m_vElemDisc[i].elemDisc->fast_ass_JM_elem(M, u);
		}
		UG_CATCH_THROW("DataEvaluator::ass_JM_elem: "
						"Cannot assemble Jacobian (M) for IElemDisc "<<i);
	}

	add_coupl_JM(M, u);
}

void DataEvaluator::add_dA_elem(LocalVector& d, LocalVector& u, GeometricObject* elem)
{
	UG_ASSERT(m_discPart & STIFF, "Using add_dA_elem, but not STIFF requested.")

	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		d.access_by_map(map(i));

		if(m_vElemDisc[i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	assemble dA
		try{
			m_vElemDisc[i].elemDisc->fast_ass_dA_elem(d, u);
		}
		UG_CATCH_THROW("DataEvaluator::ass_dA_elem: "
						"Cannot assemble Defect (A) for IElemDisc "<<i);
	}
}

void DataEvaluator::add_dM_elem(LocalVector& d, LocalVector& u, GeometricObject* elem)
{
	UG_ASSERT(m_discPart & MASS, "Using add_dM_elem, but not MASS requested.")

	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		d.access_by_map(map(i));

		if(m_vElemDisc[i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	assemble dM
		try{
			if(!m_vElemDisc[i].elemDisc->is_stationary())
				m_vElemDisc[i].elemDisc->fast_ass_dM_elem(d, u);
		}
		UG_CATCH_THROW("DataEvaluator::ass_dM_elem: "
						"Cannot assemble Defect (M) for IElemDisc "<<i);
	}
}

void DataEvaluator::add_rhs_elem(LocalVector& rhs, GeometricObject* elem)
{
	UG_ASSERT(m_discPart & RHS, "Using add_rhs_elem, but not RHS requested.")

	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	access disc functions
		rhs.access_by_map(map(i));

		if(m_vElemDisc[i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	assemble rhs
		try{
			m_vElemDisc[i].elemDisc->fast_ass_rhs_elem(rhs);
		}
		UG_CATCH_THROW("DataEvaluator::ass_rhs_elem: "
						"Cannot assemble rhs for IElemDisc "<<i);
	}
}

void DataEvaluator::finish_elem_loop()
{
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
		try{
			m_vElemDisc[i].elemDisc->fast_finish_elem_loop();
		}
		UG_CATCH_THROW("DataEvaluator::finish_element_loop: "
						"Cannot finish element loop for IElemDisc "<<i);
	}

	clear_positions_in_user_data();
}

void DataEvaluator::clear_positions_in_user_data()
{
//	remove ip series for all used UserData
	for(size_t i = 0; i < m_vConstData.size(); ++i) m_vConstData[i]->clear();
	for(size_t i = 0; i < m_vPosData.size(); ++i)   m_vPosData[i]->clear();
	for(size_t i = 0; i < m_vDependentData.size(); ++i) m_vDependentData[i]->clear();

//	we remove all ips, since they may have been set in prepare_elem
	for(size_t d = 0; d < m_vElemDisc.size(); ++d)
		for(size_t i = 0; i < m_vElemDisc[d].elemDisc->num_imports(); ++i)
			m_vElemDisc[d].elemDisc->get_import(i).clear_ips();
}

///////////////////////////////////////////////////////////////////////////////
// Coupling
///////////////////////////////////////////////////////////////////////////////

void DataEvaluator::add_coupl_JA(LocalMatrix& J, LocalVector& u)
{
//	compute linearized defect
	for(size_t i = 0; i < m_vImport[STIFF].size(); ++i)
	{
	//	set correct access for import
		u.access_by_map(m_vImport[STIFF][i].map);

	//	compute linearization of defect
		try{m_vImport[STIFF][i].import->compute_lin_defect(u);
		}UG_CATCH_THROW("DataEvaluator::add_coupl_JA: Cannot compute"
						" linearized defect for Import " << i <<" (Stiffness part).");
	}
//	compute linearized defect
	for(size_t i = 0; i < m_vImport[RHS].size(); ++i)
	{
	//	set correct access for import
		u.access_by_map(m_vImport[RHS][i].map);

	//	compute linearization of defect
		try{m_vImport[RHS][i].import->compute_lin_defect(u);
		}UG_CATCH_THROW("DataEvaluator::add_coupl_JA: Cannot compute"
						" linearized defect for Import " << i <<" (Rhs part).");
	}

//	loop all imports located in the stiffness part
	for(size_t i = 0; i < m_vImport[STIFF].size(); ++i)
	{
	//	rows are given by import, cols are given by connected data
		J.access_by_map(m_vImport[STIFF][i].map, m_vImport[STIFF][i].connMap);

	//	add off diagonal coupling
		try{m_vImport[STIFF][i].import->add_jacobian(J, 1.0);}
		UG_CATCH_THROW("DataEvaluator::add_coupl_JA: Cannot add couplings.");
	}

//	loop all imports located in the rhs part
	for(size_t i = 0; i < m_vImport[RHS].size(); ++i)
	{
	//	rows are given by import, cols are given by connected data
		J.access_by_map(m_vImport[RHS][i].map, m_vImport[RHS][i].connMap);

	//	add off diagonal coupling
		try{m_vImport[RHS][i].import->add_jacobian(J, -1.0);}
		UG_CATCH_THROW("DataEvaluator::add_coupl_JA: Cannot add couplings.");
	}
}

void DataEvaluator::add_coupl_JM(LocalMatrix& J, LocalVector& u)
{
//	compute linearized defect
	for(size_t i = 0; i < m_vImport[MASS].size(); ++i)
	{
	//	set correct access for import
		u.access_by_map(m_vImport[MASS][i].map);

	//	compute linearization of defect
		try{m_vImport[MASS][i].import->compute_lin_defect(u);
		}UG_CATCH_THROW("DataEvaluator::add_coupl_JM: Cannot compute"
						" linearized defect for Import " << i <<" (Mass part).");
	}

//	loop all imports located in the mass part
	for(size_t i = 0; i < m_vImport[MASS].size(); ++i)
	{
	//	rows are given by import, cols are given by connected data
		J.access_by_map(m_vImport[MASS][i].map, m_vImport[MASS][i].connMap);

	//	add off diagonal coupling
		try{
			m_vImport[MASS][i].import->add_jacobian(J, 1.0);
		}
		UG_CATCH_THROW("DataEvaluator::add_coupl_JM: Cannot add couplings.");
	}
}

} // end namespace ug

