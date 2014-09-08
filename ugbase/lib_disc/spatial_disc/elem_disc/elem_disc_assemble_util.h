/*
 * elem_disc_assemble_util.h
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
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Stiffness matrix
 * \param[in]		u				solution
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleStiffnessMatrix(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                        	ConstSmartPtr<TDomain> spDomain,
                        	ConstSmartPtr<DoFDistribution> dd,
        					TIterator iterBegin,
        					TIterator iterEnd,
                        	int si, bool bNonRegularGrid,
                        	typename TAlgebra::matrix_type& A,
                        	const typename TAlgebra::vector_type& u,
                        	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

	try
	{
	DataEvaluator<TDomain> Eval(STIFF,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locA;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locA.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind, true);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix Cannot prepare element.");

	// 	Assemble JA
		locA = 0.0;
		try
		{
			Eval.add_jac_A_elem(locA, locU, elem, vCornerCoords);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot compute Jacobian (A).");

	//	send local to global matrix
		try{
			spAssTuner->add_local_mat_to_global(A,locA,dd);
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot add local matrix.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot finish element loop.");

	}
	UG_CATCH_THROW("AssembleStiffnessMatrix': Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleStiffnessMatrix(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                        	ConstSmartPtr<TDomain> spDomain,
                        	ConstSmartPtr<DoFDistribution> dd,
                        	int si, bool bNonRegularGrid,
                        	typename TAlgebra::matrix_type& A,
                        	const typename TAlgebra::vector_type& u,
                        	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleStiffnessMatrix<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, A, u, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleStiffnessMatrix<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, A, u, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble Mass Matrix
////////////////////////////////////////////////////////////////////////////////
/**
 * This function adds to the Mass matrix the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	M				Mass matrix
 * \param[in]		u				solution
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleMassMatrix(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                   	ConstSmartPtr<TDomain> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& M,
					const typename TAlgebra::vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(MASS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locM;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locM.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind, true);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot prepare element.");

	// 	Assemble JM
		locM = 0.0;
		try
		{
			Eval.add_jac_M_elem(locM, locU, elem, vCornerCoords);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot compute Jacobian (M).");

	// send local to global matrix
		try{
			spAssTuner->add_local_mat_to_global(M, locM, dd);
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot add local matrix.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("AssembleMassMatrix: Cannot finish element loop.");

	}
	UG_CATCH_THROW("AssembleMassMatrix: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleMassMatrix(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                   	ConstSmartPtr<TDomain> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& M,
					const typename TAlgebra::vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleMassMatrix<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, M, u, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleMassMatrix<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, M, u, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Jacobian
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the jacobian the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		u				solution
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleJacobian(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                 	ConstSmartPtr<TDomain> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locJ;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locJ.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind, true);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot prepare element.");

	//	reset local algebra
		locJ = 0.0;

	// 	Assemble JA
		try
		{
			Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot compute Jacobian (A).");

	// send local to global matrix
		try{
			spAssTuner->add_local_mat_to_global(J, locJ, dd);
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot add local matrix.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleJacobian(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                 	ConstSmartPtr<TDomain> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleJacobian<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, J, u, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleJacobian<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, J, u, spAssTuner);
	}
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
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		vSol			current and previous solutions
 * \param[in]		s_a0			scaling factor for stiffness part
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleJacobian(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                 	ConstSmartPtr<TDomain> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
					number s_a0,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	get current time and vector
	const typename TAlgebra::vector_type& u = *vSol->solution(0);

//	create local time series
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(MASS | STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid,
					   &locTimeSeries);
	Eval.set_time_point(0);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

//	local algebra
	LocalIndices ind; LocalVector locU; LocalMatrix locJ;

	EL_PROFILE_BEGIN(Elem_AssembleJacobian);
// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locJ.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//	read local values of time series
		if(Eval.time_series_needed())
			locTimeSeries.read_values(vSol, ind);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind, true);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot prepare element.");

	//	reset local algebra
		locJ = 0.0;

		//EL_PROFILE_BEGIN(Elem_add_JA);
		// 	Assemble JA
		try
		{
			Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_INSTATIONARY);
			locJ *= s_a0;

			Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_STATIONARY);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (A).");
		//EL_PROFILE_END();

	// 	Assemble JM
		try
		{
			Eval.add_jac_M_elem(locJ, locU, elem, vCornerCoords, PT_INSTATIONARY);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (M).");

	// send local to global matrix
		try{
			spAssTuner->add_local_mat_to_global(J, locJ, dd);
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot add local matrix.");

	}
	EL_PROFILE_END();

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleJacobian(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                 	ConstSmartPtr<TDomain> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
					number s_a0,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);
	
		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleJacobian<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, J, vSol, s_a0, spAssTuner);
	}
	else
	{
		AssembleJacobian<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, J, vSol, s_a0, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Defect
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the defect the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		u				solution
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleDefect(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
               	const typename TAlgebra::vector_type& u,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least one element exists, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locD, tmpLocD;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locD.resize(ind); tmpLocD.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//  ANALOG to 'domain_disc_elem()' -  modifies the solution, used
	//	for computing the defect
		if( spAssTuner->modify_solution_enabled() )
		{
			LocalVector& modLocU = locU;
			try{
				spAssTuner->modify_LocalSol(modLocU, locU, dd);
			} UG_CATCH_THROW("Cannot modify local solution.");

		// recopy modified LocalVector:
			locU = modLocU;
		}

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot prepare element.");

	//	reset local algebra
		locD = 0.0;

	// 	Assemble A
		try
		{
			Eval.add_def_A_elem(locD, locU, elem, vCornerCoords);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute Defect (A).");

	// 	Assemble rhs
		try
		{
			tmpLocD = 0.0;
			Eval.add_rhs_elem(tmpLocD, elem, vCornerCoords);
			locD.scale_append(-1, tmpLocD);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute Rhs.");

	// 	send local to global defect
		try{
			spAssTuner->add_local_vec_to_global(d, locD, dd);
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot add local vector.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(stationary) AssembleDefect: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(stationary) AssembleDefect: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleDefect(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
               	const typename TAlgebra::vector_type& u,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleDefect<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, d, u, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleDefect<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, d, u, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Defect
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the defect the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleDefect(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff,
				ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	check time scheme
	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleDefect: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleDefect: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

//	create local time series
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(MASS | STIFF | RHS | EXPL,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid,
	                   &locTimeSeries, &vScaleMass, &vScaleStiff);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locD, tmpLocD;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

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
			Eval.set_time_point(t);

		// 	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, vCornerCoords, ind, false);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot prepare element.");

		// 	Assemble M
			try
			{
				tmpLocD = 0.0;
				Eval.add_def_M_elem(tmpLocD, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locD.scale_append(vScaleMass[t], tmpLocD);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Defect (M).");

		// 	Assemble A
			try
			{
				tmpLocD = 0.0;
				Eval.add_def_A_elem(tmpLocD, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locD.scale_append(vScaleStiff[t], tmpLocD);

				if(t == 0)
					Eval.add_def_A_elem(locD, locU, elem, vCornerCoords, PT_STATIONARY);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Defect (A).");


			// Assemble defect for explicit reaction_rate, reaction and source
			if( t == 1 ) // only valid at lowest timediscretization order
			{

			   tmpLocD = 0.0;
			   try
			   {
			     Eval.add_def_A_expl_elem(tmpLocD, locU, elem, vCornerCoords, PT_INSTATIONARY);
			   }
			   UG_CATCH_THROW("(instationary) AssembleDefect explizit: Cannot compute Defect (A).");

//			   UG_ASSERT(vScaleStiff.size() == 2, "Only one step method supported.");
			   const number dt = vSol->time(0)-vSol->time(1);
			   locD.scale_append(dt, tmpLocD);
			}


		// 	Assemble rhs
			try
			{
				tmpLocD = 0.0;
				Eval.add_rhs_elem(tmpLocD, elem, vCornerCoords, PT_INSTATIONARY);
				locD.scale_append( -vScaleStiff[t], tmpLocD);

				if(t==0){
					tmpLocD = 0.0;
					Eval.add_rhs_elem(tmpLocD, elem, vCornerCoords, PT_STATIONARY);
					locD.scale_append( -1.0, tmpLocD);
				}
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Rhs.");
		}

	// 	send local to global defect
		try{
			spAssTuner->add_local_vec_to_global(d, locD, dd);
		}
		UG_CATCH_THROW("(instationary) AssembleDefect: Cannot add local vector.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleDefect: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(instationary) AssembleDefect: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleDefect(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff,
				ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);
		
		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleDefect<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleDefect<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Linear
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Matrix and to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleLinear(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locRhs; LocalMatrix locA;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locRhs.resize(ind); locA.resize(ind);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locRhs, elem, vCornerCoords, ind, true);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot prepare element.");

	//	reset local algebra
		locA = 0.0;
		locRhs = 0.0;

	// 	Assemble JA
		try
		{
			Eval.add_jac_A_elem(locA, locRhs, elem, vCornerCoords);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute Jacobian (A).");

	// 	Assemble rhs
		try
		{
			Eval.add_rhs_elem(locRhs, elem, vCornerCoords);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute Rhs.");

	// 	send local to global matrix & rhs
		try{
			spAssTuner->add_local_mat_to_global(A, locA, dd);
			spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot add local vector/matrix.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(stationary) AssembleLinear: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(stationary) AssembleLinear: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleLinear(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleLinear<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, A, rhs, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleLinear<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, A, rhs, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Linear
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Matrix and to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		mark			BoolMarker
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleLinear(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               	const std::vector<number>& vScaleMass,
               	const std::vector<number>& vScaleStiff,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	check time scheme
	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleLinear: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleLinear: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

//	create local time solution
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(MASS | STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid,
					   &locTimeSeries, &vScaleMass, &vScaleStiff);

//	prepare loop
	Eval.prepare_elem_loop(id, si);

//	local algebra
	LocalIndices ind; LocalVector locRhs, tmpLocRhs; LocalMatrix locA, tmpLocA;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locRhs.resize(ind); tmpLocRhs.resize(ind);
		locA.resize(ind); tmpLocA.resize(ind);

	//	read local values of time series
		locTimeSeries.read_values(vSol, ind);
		Eval.set_time_point(0);

	//	reset element contribution
		locA = 0.0; locRhs = 0.0;

	/////////////////////
	//	current time step

	//	get local solution at time point
		LocalVector& locU = locTimeSeries.solution(0);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind, true);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

	// 	Assemble JM
		try
		{
			tmpLocA = 0.0;
			Eval.add_jac_M_elem(tmpLocA, locU, elem, vCornerCoords, PT_INSTATIONARY);
			locA.scale_append(vScaleMass[0], tmpLocA);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (M).");

	// 	Assemble JA
		try
		{
			tmpLocA = 0.0;
			Eval.add_jac_A_elem(tmpLocA, locU, elem, vCornerCoords, PT_INSTATIONARY);
			locA.scale_append(vScaleStiff[0], tmpLocA);

			Eval.add_jac_A_elem(locA, locU, elem, vCornerCoords, PT_STATIONARY);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

	// 	Assemble rhs
		try
		{
			tmpLocRhs = 0.0;
			Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
			locRhs.scale_append(vScaleStiff[0], tmpLocRhs);

			Eval.add_rhs_elem(locRhs, elem, vCornerCoords, PT_STATIONARY);
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");

	///////////////////
	//	old time steps

	//	loop all old time points
		for(size_t t = 1; t < vScaleStiff.size(); ++t)
		{
		//	get local solution at time point
			LocalVector& locU = locTimeSeries.solution(t);
			Eval.set_time_point(t);

		// 	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, vCornerCoords, ind, false);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

		// 	Assemble dM
			try
			{
				tmpLocRhs = 0.0;
				Eval.add_def_M_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locRhs.scale_append(-vScaleMass[t], tmpLocRhs);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (M).");

		// 	Assemble dA
			try
			{
				tmpLocRhs = 0.0;
				Eval.add_def_A_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

		// 	Assemble rhs
			try
			{
				tmpLocRhs = 0.0;
				Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
				locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");
		}

		// 	send local to global matrix & rhs
			try{
				spAssTuner->add_local_mat_to_global(A, locA, dd);
				spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot add local vector/matrix.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleLinear: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(instationary) AssembleLinear: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleLinear(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               	const std::vector<number>& vScaleMass,
               	const std::vector<number>& vScaleStiff,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleLinear<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, A, rhs, vSol, vScaleMass, vScaleStiff, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleLinear<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, A, rhs, vSol, vScaleMass, vScaleStiff, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble Rhs
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		u				solution
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleRhs(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
            	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
               	const typename TAlgebra::vector_type& u,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare loop
	Eval.prepare_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locRhs;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind); locRhs.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind);
		}
		UG_CATCH_THROW("AssembleRhs: Cannot prepare element.");

	//	reset local algebra
		locRhs = 0.0;

	// 	Assemble rhs
		try
		{
			Eval.add_rhs_elem(locRhs, elem, vCornerCoords);
		}
		UG_CATCH_THROW("AssembleRhs: Cannot compute Rhs.");

	// 	send local to global rhs
		try{
			spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
		}
		UG_CATCH_THROW("AssembleRhs: Cannot add local vector.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("AssembleRhs: Cannot finish element loop.");

	}
	UG_CATCH_THROW("AssembleRhs: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleRhs(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
            	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
               	const typename TAlgebra::vector_type& u,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleRhs<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, rhs, u, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleRhs<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, rhs, u, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Rhs
////////////////////////////////////////////////////////////////////////////////

/**
 * This function adds to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleRhs(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
            	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               	const std::vector<number>& vScaleMass,
               	const std::vector<number>& vScaleStiff,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	check time scheme
	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleRhs: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleRhs: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

//	get current time
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(MASS | STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid,
					   &locTimeSeries, &vScaleMass, &vScaleStiff);

//	prepare loop
	Eval.prepare_elem_loop(id, si);

//	local algebra
	LocalIndices ind; LocalVector locRhs, tmpLocRhs;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locRhs.resize(ind); tmpLocRhs.resize(ind);

	//	read local values of time series
		locTimeSeries.read_values(vSol, ind);
		Eval.set_time_point(0);

	//	reset element contribution
		locRhs = 0.0;

	/////////////////////
	//	current time step

	//	get local solution at time point
		LocalVector& locU = locTimeSeries.solution(0);

	// 	prepare element
		try
		{
			Eval.prepare_elem(locU, elem, vCornerCoords, ind, false);
		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element.");

	// 	Assemble rhs
		try
		{
			tmpLocRhs = 0.0;
			Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
			locRhs.scale_append(vScaleStiff[0], tmpLocRhs);

			Eval.add_rhs_elem(locRhs, elem, vCornerCoords, PT_STATIONARY);
		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Rhs.");

	///////////////////
	//	old time steps

	//	loop all old time points
		for(size_t t = 1; t < vScaleStiff.size(); ++t)
		{
		//	get local solution at time point
			LocalVector& locU = locTimeSeries.solution(t);
			Eval.set_time_point(t);

		// 	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, vCornerCoords, ind, false);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element.");

		// 	Assemble dM
			try
			{
				tmpLocRhs = 0.0;
				Eval.add_def_M_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locRhs.scale_append(-vScaleMass[t], tmpLocRhs);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Jacobian (M).");

		// 	Assemble dA
			try
			{
				tmpLocRhs = 0.0;
				Eval.add_def_A_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Jacobian (A).");

		// 	Assemble rhs
			try
			{
				tmpLocRhs = 0.0;
				Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
				locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Rhs.");
		}

		// 	send local to global rhs
			try{
				spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot add local vector.");
	}

// 	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleRhs: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(instationary) AssembleRhs: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleRhs(	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
            	ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               	const std::vector<number>& vScaleMass,
               	const std::vector<number>& vScaleStiff,
               	ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		AssembleRhs<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		AssembleRhs<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Prepare Timestep (instationary)
////////////////////////////////////////////////////////////////////////////////

/**
 * This function calls the function "prepare_timestep_elem" of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vSol			current and previous solutions
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
PrepareTimestep(const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
               	int si, bool bNonRegularGrid,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
                ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	get current time and vector
	const typename TAlgebra::vector_type& u = *vSol->solution(0);

//	create local time series
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

	try
	{
	DataEvaluator<TDomain> Eval(NONE,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid,
					   &locTimeSeries);
	Eval.set_time_point(0);

//	prepare element loop
	Eval.prepare_elem_loop(id, si);

//	local algebra
	LocalIndices ind; LocalVector locU;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//	read local values of time series
		if(Eval.time_series_needed())
			locTimeSeries.read_values(vSol, ind);

	// 	prepare timestep
		try
		{
			Eval.prepare_timestep_elem(vSol->time(0), locU, elem, vCornerCoords);
		}
		UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot prepare timestep.");
	}

	}
	UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
PrepareTimestep(const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                ConstSmartPtr<TDomain> spDomain,
               	ConstSmartPtr<DoFDistribution> dd,
               	int si, bool bNonRegularGrid,
                ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
                ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		PrepareTimestep<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, vSol, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		PrepareTimestep<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, vSol, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Finish Timestep (instationary)
////////////////////////////////////////////////////////////////////////////////

/**
 * This function calls the function "finish_timestep_elem" of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vSol			current and previous solutions
 * \param[in]		spAssTuner		assemble adapter
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
FinishTimestep(const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               ConstSmartPtr<TDomain> spDomain,
               ConstSmartPtr<DoFDistribution> dd,
               TIterator iterBegin,
			   TIterator iterEnd,
               int si, bool bNonRegularGrid,
               ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
// 	check if at least on element exist, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	get current time and vector
	const typename TAlgebra::vector_type& u = *vSol->solution(0);

//	create local time series
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(NONE,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid,
					   &locTimeSeries);
	Eval.set_time_point(0);

//	prepare loop
	Eval.prepare_elem_loop(id, si);

//	local algebra
	LocalIndices ind; LocalVector locU;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	//	read local values of time series
		if(Eval.time_series_needed())
			locTimeSeries.read_values(vSol, ind);

	// 	finish timestep
		try
		{
			Eval.finish_timestep_elem(locTimeSeries.time(0), locU, elem, vCornerCoords);
		}
		UG_CATCH_THROW("(instationary) FinishTimestep: Cannot finish timestep.");
	}

	}
	UG_CATCH_THROW("(instationary) FinishTimestep: Cannot create Data Evaluator");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
FinishTimestep(const std::vector<IElemDisc<TDomain>*>& vElemDisc,
               ConstSmartPtr<TDomain> spDomain,
               ConstSmartPtr<DoFDistribution> dd,
               int si, bool bNonRegularGrid,
               ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol,
               ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
	//	check if only some elements are selected
	if(spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		FinishTimestep<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, vSol, spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		FinishTimestep<TElem,TDomain,TAlgebra>
			(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, vSol, spAssTuner);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Error estimator (stationary)
////////////////////////////////////////////////////////////////////////////////

/**
 * This function assembles the error estimator associated with all the
 * element discs in the internal data structure.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		u				solution
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleErrorEstimator
(
	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
	ConstSmartPtr<TDomain> spDomain,
	ConstSmartPtr<DoFDistribution> dd,
	TIterator iterBegin,
	TIterator iterEnd,
	int si,
	bool bNonRegularGrid,
	const typename TAlgebra::vector_type& u
)
{
// 	check if at least one element exists, else return
	if(iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	prepare for given elem discs
	try
	{
	DataEvaluator<TDomain> Eval(STIFF | RHS,
	                   vElemDisc, dd->function_pattern(), bNonRegularGrid);

//	prepare element loop
	Eval.prepare_err_est_elem_loop(id, si);

// 	local indices and local algebra
	LocalIndices ind; LocalVector locU;

// 	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	// 	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	// 	adapt local algebra
		locU.resize(ind);

	// 	read local values of u
		GetLocalVector(locU, u);

	// 	prepare element
		try
		{
			Eval.prepare_err_est_elem(locU, elem, vCornerCoords, ind, false);
		}
		UG_CATCH_THROW("(stationary) AssembleRhs: Cannot prepare element.");

	// 	assemble stiffness part
		try
		{
			Eval.compute_err_est_A_elem(locU, elem, vCornerCoords, ind);
		}
		UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

	// 	assemble rhs part
		try
		{
			Eval.compute_err_est_rhs_elem(elem, vCornerCoords, ind);
		}
		UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

	}

// 	finish element loop
	try
	{
		Eval.finish_err_est_elem_loop();
	}
	UG_CATCH_THROW("AssembleErrorEstimator: Cannot finish element loop.");

	}
	UG_CATCH_THROW("AssembleErrorEstimator: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleErrorEstimator
(
	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
	ConstSmartPtr<TDomain> spDomain,
	ConstSmartPtr<DoFDistribution> dd,
	int si,
	bool bNonRegularGrid,
	const typename TAlgebra::vector_type& u
)
{
	//	general case: assembling over all elements in subset si
	AssembleErrorEstimator<TElem,TDomain,TAlgebra>
		(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si),
			si, bNonRegularGrid, u);
}

////////////////////////////////////////////////////////////////////////////////
// Error estimator (instationary)
////////////////////////////////////////////////////////////////////////////////

/**
 * This function assembles the error estimator associated with all the
 * element discs in the internal data structure.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		iterBegin		element iterator
 * \param[in]		iterEnd			element iterator
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		vSol				solution
 */
template <typename TElem, typename TDomain, typename TAlgebra, typename TIterator>
void
AssembleErrorEstimator
(
	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
	ConstSmartPtr<TDomain> spDomain,
	ConstSmartPtr<DoFDistribution> dd,
	TIterator iterBegin,
	TIterator iterEnd,
	int si,
	bool bNonRegularGrid,
	std::vector<number> vScaleMass,
	std::vector<number> vScaleStiff,
	ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol
)
{
// 	check if at least one element exists, else return
	if (iterBegin == iterEnd) return;

//	reference object id
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	storage for corner coordinates
	MathVector<TDomain::dim> vCornerCoords[TElem::NUM_VERTICES];

//	check time scheme
	if(vScaleMass.size() != vScaleStiff.size())
		UG_THROW("(instationary) AssembleErrorEstimator: s_a and s_m must have same size.");

	if(vSol->size() < vScaleStiff.size())
		UG_THROW("(instationary) AssembleErrorEstimator: Time stepping scheme needs at "
				"least "<<vScaleStiff.size()<<" time steps, but only "<<
				vSol->size() << " passed.");

//	create local time series
	LocalVectorTimeSeries locTimeSeries;
	locTimeSeries.read_times(vSol);

//	prepare for given elem discs
	try
	{
		DataEvaluator<TDomain> Eval(MASS | STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries, &vScaleMass, &vScaleStiff);

	//	prepare element loop
		Eval.prepare_err_est_elem_loop(id, si);

	// 	local indices and local algebra
		LocalIndices ind;

	// 	loop over all elements
		for (TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		// 	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	read local values of time series
			locTimeSeries.read_values(vSol, ind);

		//	loop all time points and assemble them
			for (std::size_t t = 0; t < vScaleStiff.size(); ++t)
			{
			//	get local solution at timepoint
				LocalVector& locU = locTimeSeries.solution(t);
				Eval.set_time_point(t);

			// 	prepare element
				try
				{
					Eval.prepare_err_est_elem(locU, elem, vCornerCoords, ind, false);
				}
				UG_CATCH_THROW("AssembleErrorEstimator: Cannot prepare element.");

			// 	assemble stiffness part
				try
				{
					Eval.compute_err_est_A_elem(locU, elem, vCornerCoords, ind, vScaleMass[t], vScaleStiff[t]);
				}
				UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

			// 	assemble mass part
				try
				{
					Eval.compute_err_est_M_elem(locU, elem, vCornerCoords, ind, vScaleMass[t], vScaleStiff[t]);
				}
				UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

			// 	assemble rhs part
				try
				{
					Eval.compute_err_est_rhs_elem(elem, vCornerCoords, ind, vScaleMass[t], vScaleStiff[t]);
				}
				UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");
			}
		}

		// 	finish element loop
		try
		{
			Eval.finish_err_est_elem_loop();
		}
		UG_CATCH_THROW("AssembleErrorEstimator: Cannot finish element loop.");

	}
	UG_CATCH_THROW("AssembleErrorEstimator: Cannot create Data Evaluator.");
}

template <typename TElem, typename TDomain, typename TAlgebra>
void
AssembleErrorEstimator
(
	const std::vector<IElemDisc<TDomain>*>& vElemDisc,
	ConstSmartPtr<TDomain> spDomain,
	ConstSmartPtr<DoFDistribution> dd,
	int si,
	bool bNonRegularGrid,
	std::vector<number> vScaleMass,
	std::vector<number> vScaleStiff,
	ConstSmartPtr<VectorTimeSeries<typename TAlgebra::vector_type> > vSol
)
{
	//	general case: assembling over all elements in subset si
	AssembleErrorEstimator<TElem,TDomain,TAlgebra>
		(vElemDisc, spDomain, dd, dd->template begin<TElem>(si), dd->template end<TElem>(si),
			si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
}

} // end namespace ug


#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__ */
