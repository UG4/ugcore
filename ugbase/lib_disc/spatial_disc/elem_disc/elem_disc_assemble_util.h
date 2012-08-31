/*
 * assemble_elem_disc.h
 *
 *  Created on: 08.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// intern headers
#include "../../reference_element/reference_element.h"
#include "./elem_disc_interface.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/user_data/data_evaluator.h"

#define PROFILE_ELEM_LOOP
#ifdef PROFILE_ELEM_LOOP
	#define EL_PROFILE_FUNC()		PROFILE_FUNC()
	#define EL_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define EL_PROFILE_END()		PROFILE_END()
#else
	#define EL_PROFILE_FUNC()
	#define EL_PROFILE_BEGIN(name)
	#define EL_PROFILE_END()
#endif


namespace ug {

////////////////////////////////////////////////////////////////////////////////
// Assemble Stiffness Matrix
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Stiffness matrix the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Stiffness matrix
 * \param[in]		u				solution
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleStiffnessMatrix(	const std::vector<IElemDisc*>& vElemDisc,
                        	ConstSmartPtr<TDD> dd,
                        	int si, bool bNonRegularGrid,
                        	typename TAlgebra::matrix_type& A,
                        	const typename TAlgebra::vector_type& u,
                        	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locA;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(false);
		Eval.template prepare_elem_loop<TElem>();
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("AssembleStiffnessMatrix': Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locA.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, true);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, true);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot compute element data.");

	//	Evaluate lin defect A
		try
		{
			Eval.compute_lin_defect_JA(locU, elem);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot compute lin_defect_JA.");

	// 	Assemble JA
		locA = 0.0;
		try
		{
			Eval.ass_JA_elem(locA, locU, elem);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot compute Jacobian (A).");

	//	add couplings
		try{
			Eval.add_coupl_JA(locA);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot add couplings (A).");

	// 	send local to global matrix
		AddLocalMatrixToGlobal(A, locA);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble Mass Matrix
////////////////////////////////////////////////////////////////////////////////
/**
 * This function adds to the Mass matrix the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	M				Mass matrix
 * \param[in]		u				solution
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleMassMatrix(	const std::vector<IElemDisc*>& vElemDisc,
					ConstSmartPtr<TDD> dd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& M,
					const typename TAlgebra::vector_type& u,
                	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locM;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(false);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("AssembleMassMatrix: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locM.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, true, true);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, true);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot compute element data.");

	//	Evaluate lin defect M
		try
		{
			Eval.compute_lin_defect_JM(locU, elem);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot compute lin_defect_JM.");

	// 	Assemble JM
		locM = 0.0;
		try
		{
			Eval.ass_JM_elem(locM, locU, elem);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot compute Jacobian (M).");

	//	add couplings
		try
		{
			Eval.add_coupl_JM(locM);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot add couplings (M).");

	// 	send local to global matrix
		AddLocalMatrixToGlobal(M, locM);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("AssembleMassMatrix: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Prepare Timestep (instationary)
////////////////////////////////////////////////////////////////////////////////

/**
 * This function calls the function "prepare_timestep_elem" of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vSol			current and previous solutions
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
PrepareTimestep(const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	get current time and vector
	const number time = vSol->time(0);
	const typename TAlgebra::vector_type& u = *vSol->solution(0);

//	create data evaluator
	DataEvaluator Eval;
	LocalVectorTimeSeries locTimeSeries;
	LocalIndices ind; LocalVector locU;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(true, &locTimeSeries);
		Eval.set_time(time);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//	read local values of time series
		if(Eval.time_series_needed())
		{
			locTimeSeries.read_values(vSol, ind);
			locTimeSeries.read_times(vSol);
		}

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, true, true);
		}
		UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, true);
		}
		UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot compute element data.");

	// 	prepare timestep
		try
		{
			Eval.prepare_timestep_elem(elem, locU);
		}
		UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot prepare timestep.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Jacobian
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the jacobian the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		u				solution
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleJacobian(	const std::vector<IElemDisc*>& vElemDisc,
					ConstSmartPtr<TDD> dd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u,
                	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locJ;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(false);
		Eval.template prepare_elem_loop<TElem>();
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locJ.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, true);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, true);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot compute element data.");

	//	Evaluate lin defect A
		try
		{
			Eval.compute_lin_defect_JA(locU, elem);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot compute lin_defect_JA.");

	// 	Assemble JA
		locJ = 0.0;
		try
		{
			Eval.ass_JA_elem(locJ, locU, elem);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot compute Jacobian (A).");

	//	add couplings
		try
		{
			Eval.add_coupl_JA(locJ);
		}
		UG_CATCH_THROW("(stationary AssembleJacobian: Cannot add couplings (A).");

	// 	send local to global matrix
		AddLocalMatrixToGlobal(J, locJ);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Jacobian
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the jacobian the entries of one subset for all passed
 * element discretizations.
 * Note, that it is assumed	to have s_m0 == 1
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		vSol			current and previous solutions
 * \param[in]		s_a0			scaling factor for stiffness part
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleJacobian(	const std::vector<IElemDisc*>& vElemDisc,
					ConstSmartPtr<TDD> dd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
					number s_a0,
                	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	get current time and vector
	const number time = vSol->time(0);
	const typename TAlgebra::vector_type& u = *vSol->solution(0);

//	create data evaluator
	DataEvaluator Eval;
	LocalVectorTimeSeries locTimeSeries;
	LocalIndices ind; LocalVector locU; LocalMatrix locJ;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(true, &locTimeSeries);
		Eval.set_time(time);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locJ.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//	read local values of time series
		if(Eval.time_series_needed())
		{
			locTimeSeries.read_values(vSol, ind);
			locTimeSeries.read_times(vSol);
		}

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, true, true);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, true);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute element data.");

	//	Evaluate lin defect A
		try
		{
			Eval.compute_lin_defect_JA(locU, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute lin_defect_JA.");

	//	Evaluate lin defect M
		try
		{
			Eval.compute_lin_defect_JM(locU, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute lin_defect_JM.");

	// 	Assemble JA
		locJ = 0.0;
		try
		{
			Eval.ass_JA_elem(locJ, locU, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (A).");

	//	add couplings
		try
		{
			Eval.add_coupl_JA(locJ);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot add couplings (A).");

	//	scale stiffness part
		locJ *= s_a0;

	// 	Assemble JM
		try
		{
			Eval.ass_JM_elem(locJ, locU, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (M).");

	//	add couplings
		try
		{
			Eval.add_coupl_JM(locJ);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot add couplings (M).");

	// 	send local to global matrix
		AddLocalMatrixToGlobal(J, locJ);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Defect
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the defect the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		u				solution
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleDefect(	const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
               	const typename TAlgebra::vector_type& u,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locD, tmpLocD;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(false);
		Eval.template prepare_elem_loop<TElem>();
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(stationary) AssembleDefect: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locD.resize(ind); tmpLocD.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, false);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute element data.");

	// 	Assemble A
		locD = 0.0;
		try
		{
			Eval.ass_dA_elem(locD, locU, elem);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute Defect (A).");

	// 	Assemble rhs
		tmpLocD = 0.0;
		try
		{
			Eval.ass_rhs_elem(tmpLocD, elem);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute Rhs.");

		locD.scale_append(-1, tmpLocD);

	// 	send local to global rhs
		AddLocalVector(d, locD);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(stationary) AssembleDefect: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Defect
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the defect the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleDefect(	const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

//	current time
	LocalVectorTimeSeries locTimeSeries;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locD, tmpLocD;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(true, &locTimeSeries);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(instationary) AssembleDefect: Cannot prepare element loop.");

	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleDefect: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleDefect: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

//	read time points
	locTimeSeries.read_times(vSol);

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locD.resize(ind); tmpLocD.resize(ind);

	//	read local values of time series
		locTimeSeries.read_values(vSol, ind);

	//	reset contribution of this element
		locD = 0.0;

	//	loop all time points and assemble them
		for(size_t t = 0; t < vScaleStiff.size(); ++t)
		{
		//	get local solution at timepoint
			LocalVector& locU = locTimeSeries.solution(t);
			const number time = vSol->time(t);
			Eval.set_time(time);

		// 	prepare element
			try
			{
				Eval.prepare_elem(elem, locU, ind, false, true);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot prepare element.");

		//	Compute element data
			try
			{
				Eval.compute_elem_data(locU, elem, false);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute element data.");

		// 	Assemble M
			tmpLocD = 0.0;
			try
			{
				Eval.ass_dM_elem(tmpLocD, locU, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Defect (M).");

			locD.scale_append(vScaleMass[t], tmpLocD);

		// 	Assemble A
			tmpLocD = 0.0;
			try
			{
				Eval.ass_dA_elem(tmpLocD, locU, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Defect (A).");

			locD.scale_append(vScaleStiff[t], tmpLocD);

		// 	Assemble rhs
			tmpLocD = 0.0;
			try
			{
				Eval.ass_rhs_elem(tmpLocD, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Rhs.");

			locD.scale_append( -vScaleStiff[t], tmpLocD);
		}

	// 	send local to global rhs
		AddLocalVector(d, locD);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleDefect: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Linear
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Matrix and to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD,typename TAlgebra>
void
AssembleLinear(	const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locRhs; LocalMatrix locA;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(false);
		Eval.template prepare_elem_loop<TElem>();
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(stationary) AssembleLinear: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locRhs.resize(ind); locA.resize(ind);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locRhs, ind);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locRhs, elem, false);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute element data.");

	// 	Assemble JA
		locA = 0.0;
		try
		{
			Eval.ass_JA_elem(locA, locRhs, elem);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute Jacobian (A).");

	// 	Assemble rhs
		locRhs = 0.0;
		try
		{
			Eval.ass_rhs_elem(locRhs, elem);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute Rhs.");

	// 	send local to global matrix
		AddLocalMatrixToGlobal(A, locA);

	// 	send local to global rhs
		AddLocalVector(rhs, locRhs);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(stationary) AssembleLinear: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Linear
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Matrix and to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleLinear(	const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               	const std::vector<number>& vScaleMass,
               	const std::vector<number>& vScaleStiff,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

//	get current time
	LocalVectorTimeSeries locTimeSeries;
	LocalIndices ind; LocalVector locRhs, tmpLocRhs; LocalMatrix locA, tmpLocA;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(true, &locTimeSeries);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element loop.");

//	get time points
	locTimeSeries.read_times(vSol);

	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleLinear: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleLinear: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locRhs.resize(ind); tmpLocRhs.resize(ind);
		locA.resize(ind); tmpLocA.resize(ind);

	//	read local values of time series
		locTimeSeries.read_values(vSol, ind);
		number time = vSol->time(0);
		Eval.set_time(time);

	//	reset element contribution
		locA = 0.0; locRhs = 0.0;

	/////////////////////
	//	current time step

	//	get local solution at time point
		LocalVector& locU = locTimeSeries.solution(0);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, false, true);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, false);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute element data.");

	// 	Assemble JM
		tmpLocA = 0.0;
		try
		{
			Eval.ass_JM_elem(tmpLocA, locU, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (M).");

		locA.scale_append(vScaleMass[0], tmpLocA);

	// 	Assemble JA
		tmpLocA = 0.0;
		try
		{
			Eval.ass_JA_elem(tmpLocA, locU, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

		locA.scale_append(vScaleStiff[0], tmpLocA);

	// 	Assemble rhs
		tmpLocRhs = 0.0;
		try
		{
			Eval.ass_rhs_elem(tmpLocRhs, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");

		locRhs.scale_append(vScaleStiff[0], tmpLocRhs);

	///////////////////
	//	old time steps

	//	loop all old time points
		for(size_t t = 1; t < vScaleStiff.size(); ++t)
		{
		//	get local solution at time point
			LocalVector& locU = locTimeSeries.solution(t);
			number time = vSol->time(t);
			Eval.set_time(time);

		// 	prepare element
			try
			{
				Eval.prepare_elem(elem, locU, ind, false, true);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

		//	Compute element data
			try
			{
				Eval.compute_elem_data(locU, elem, false);

			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute element data");

		// 	Assemble dM
			tmpLocRhs = 0.0;
			try
			{
				Eval.ass_dM_elem(tmpLocRhs, locU, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (M).");

			locRhs.scale_append(-vScaleMass[t], tmpLocRhs);

		// 	Assemble dA
			tmpLocRhs = 0.0;
			try
			{
				Eval.ass_dA_elem(tmpLocRhs, locU, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

			locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);

		// 	Assemble rhs
			tmpLocRhs = 0.0;
			try
			{
				Eval.ass_rhs_elem(tmpLocRhs, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");

			locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
		}

	// 	send local to global matrix
		AddLocalMatrixToGlobal(A, locA);

	// 	send local to global rhs
		AddLocalVector(rhs, locRhs);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleLinear: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble Rhs
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		u				solution
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleRhs(	const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
               	const typename TAlgebra::vector_type& u,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locRhs;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(false);
		Eval.template prepare_elem_loop<TElem>();
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("AssembleRhs: Cannot prepare element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locRhs.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind);
		}
		UG_CATCH_THROW("AssembleRhs: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, false);
		}
		UG_CATCH_THROW("AssembleRhs: Cannot compute element data.");

	// 	Assemble rhs
		locRhs = 0.0;
		try
		{
			Eval.ass_rhs_elem(locRhs, elem);
		}
		UG_CATCH_THROW("AssembleRhs: Cannot compute Rhs.");

	// 	send local to global rhs
		AddLocalVector(rhs, locRhs);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("AssembleRhs: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Rhs
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
AssembleRhs(	const std::vector<IElemDisc*>& vElemDisc,
               	ConstSmartPtr<TDD> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               	const std::vector<number>& vScaleMass,
               	const std::vector<number>& vScaleStiff,
            	BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	create data evaluator
	DataEvaluator Eval;

//	get current time
	LocalVectorTimeSeries locTimeSeries;
	LocalIndices ind; LocalVector locRhs, tmpLocRhs;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(true, &locTimeSeries);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element loop.");

//	get time points
	locTimeSeries.read_times(vSol);

	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleRhs: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleRhs: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locRhs.resize(ind); tmpLocRhs.resize(ind);

	//	read local values of time series
		locTimeSeries.read_values(vSol, ind);
		number time = vSol->time(0);
		Eval.set_time(time);

	//	reset element contribution
		locRhs = 0.0;

	/////////////////////
	//	current time step

	//	get local solution at time point
		LocalVector& locU = locTimeSeries.solution(0);

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, false, true);
		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, false);
		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute element data.");

	// 	Assemble rhs
		tmpLocRhs = 0.0;
		try
		{
			Eval.ass_rhs_elem(tmpLocRhs, elem);
		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Rhs.");

		locRhs.scale_append(vScaleStiff[0], tmpLocRhs);

	///////////////////
	//	old time steps

	//	loop all old time points
		for(size_t t = 1; t < vScaleStiff.size(); ++t)
		{
		//	get local solution at time point
			LocalVector& locU = locTimeSeries.solution(t);
			number time = vSol->time(t);
			Eval.set_time(time);

		// 	prepare element
			try
			{
				Eval.prepare_elem(elem, locU, ind, false, true);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element.");

		//	Compute element data
			try
			{
				Eval.compute_elem_data(locU, elem, false);

			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute element data");

		// 	Assemble dM
			tmpLocRhs = 0.0;
			try
			{
				Eval.ass_dM_elem(tmpLocRhs, locU, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Jacobian (M).");

			locRhs.scale_append(-vScaleMass[t], tmpLocRhs);

		// 	Assemble dA
			tmpLocRhs = 0.0;
			try
			{
				Eval.ass_dA_elem(tmpLocRhs, locU, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Jacobian (A).");

			locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);

		// 	Assemble rhs
			tmpLocRhs = 0.0;
			try
			{
				Eval.ass_rhs_elem(tmpLocRhs, elem);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Rhs.");

			locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
		}

	// 	send local to global rhs
		AddLocalVector(rhs, locRhs);
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleRhs: Cannot finish element loop.");
}

////////////////////////////////////////////////////////////////////////////////
// Finish Timestep (instationary)
////////////////////////////////////////////////////////////////////////////////

/**
 * This function calls the function "finish_timestep_elem" of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vSol			current and previous solutions
 * \param[in]		sel				Selector
 */
template <typename TElem, typename TDD, typename TAlgebra>
void
FinishTimestep(const std::vector<IElemDisc*>& vElemDisc,
               ConstSmartPtr<TDD> dd,
               int si, bool bNonRegularGrid,
               ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               BoolMarker* sel = NULL)
{
//	get element iterator
	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dd->template begin<TElem>(si);
	iterEnd = dd->template end<TElem>(si);

// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	get current time and vector
	const number time = vSol->time(0);
	const typename TAlgebra::vector_type& u = *vSol->solution(0);

//	create data evaluator
	DataEvaluator Eval;
	LocalVectorTimeSeries locTimeSeries;
	LocalIndices ind; LocalVector locU;

//	prepare for given elem discs
	try
	{
		Eval.set_elem_discs(vElemDisc, dd->function_pattern());
		Eval.set_non_regular_grid(bNonRegularGrid);
		Eval.set_time_dependent(true, &locTimeSeries);
		Eval.set_time(time);
		Eval.template prepare_elem_loop<TElem>(true);
		Eval.set_subset(si);
	}
	UG_CATCH_THROW("(instationary) FinishTimestep: Cannot finish element loop.");

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel) if(!sel->is_marked(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//	read local values of time series
		if(Eval.time_series_needed())
		{
			locTimeSeries.read_values(vSol, ind);
			locTimeSeries.read_times(vSol);
		}

	// 	prepare element
		try
		{
			Eval.prepare_elem(elem, locU, ind, true, true);
		}
		UG_CATCH_THROW("(instationary) FinishTimestep: Cannot prepare element.");

	//	Compute element data
		try
		{
			Eval.compute_elem_data(locU, elem, true);
		}
		UG_CATCH_THROW("(instationary) FinishTimestep: Cannot compute element data.");

	// 	finish timestep
		try
		{
			Eval.finish_timestep_elem(elem, time, locU);
		}
		UG_CATCH_THROW("(instationary) FinishTimestep: Cannot finish timestep.");
	}
}

} // end namespace ug


#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__ */
