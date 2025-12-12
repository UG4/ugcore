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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"

// intern headers
//#include "../../reference_element/reference_element.h"
#include "./elem_disc_interface.h"
//#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/user_data/data_evaluator.h"
#include "bridge/util_algebra_dependent.h"

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

extern DebugID DID_ELEM_DISC_ASSEMBLE_UTIL;

/// Global assembler based on the straightforward application of the local discretizations
/**
 * Template class of the global assembler that applies the local (element)
 * discretizations to the elements in given subsets and adds the local data
 * to the global ones.
 *
 * \tparam TDomain		domain type
 * \tparam TAlgebra		algebra type
 */
template <typename TDomain, typename TAlgebra>
class StdGlobAssembler
{
	///	Domain type
	using domain_type = TDomain;
	
	///	Algebra type
	using algebra_type = TAlgebra;
	
	///	Vector type in the algebra
	using vector_type = typename algebra_type::vector_type;
	
	///	Matrix type in the algebra
	using matrix_type = typename algebra_type::matrix_type;
	
////////////////////////////////////////////////////////////////////////////////
// Assemble Stiffness Matrix
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Stiffness matrix. (This version
	 * processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in,out]	A				Stiffness matrix
	 * \param[in]		u				solution
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleStiffnessMatrix(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
								ConstSmartPtr<domain_type> spDomain,
								ConstSmartPtr<DoFDistribution> dd,
								TIterator iterBegin,
								TIterator iterEnd,
								int si, bool bNonRegularGrid,
								matrix_type& A,
								const vector_type& u,
								ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

		try
		{
		DataEvaluator<domain_type> Eval(STIFF,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locU; LocalMatrix locA;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind); locA.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
			}
			UG_CATCH_THROW("AssembleStiffnessMatrix Cannot prepare element.");

		//	Assemble JA
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

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot finish element loop.");

		}
		UG_CATCH_THROW("AssembleStiffnessMatrix': Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble Mass Matrix
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Mass matrix. (This version
	 * processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in,out]	M				Mass matrix
	 * \param[in]		u				solution
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleMassMatrix( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<domain_type> spDomain,
						ConstSmartPtr<DoFDistribution> dd,
						TIterator iterBegin,
						TIterator iterEnd,
						int si, bool bNonRegularGrid,
						matrix_type& M,
						const vector_type& u,
						ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(MASS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locU; LocalMatrix locM;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind); locM.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
			}
			UG_CATCH_THROW("AssembleMassMatrix: Cannot prepare element.");

		//	Assemble JM
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

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot finish element loop.");

		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Jacobian
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Jacobian in the stationary case. (This version
	 * processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in,out]	J				jacobian
	 * \param[in]		u				solution
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleJacobian(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<domain_type> spDomain,
						ConstSmartPtr<DoFDistribution> dd,
						TIterator iterBegin,
						TIterator iterEnd,
						int si, bool bNonRegularGrid,
						matrix_type& J,
						const vector_type& u,
						ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locU; LocalMatrix locJ;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind); locJ.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
			}
			UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot prepare element.");

		//	reset local algebra
			locJ = 0.0;

		//	Assemble JA
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

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Jacobian
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Jacobian in the time-dependent case.
	 * Note, that s_m0 == 1
	 * (This version processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
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
	template <typename TElem, typename TIterator>
	static void
	AssembleJacobian(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<domain_type> spDomain,
						ConstSmartPtr<DoFDistribution> dd,
						TIterator iterBegin,
						TIterator iterEnd,
						int si, bool bNonRegularGrid,
						matrix_type& J,
						ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
						number s_a0,
						ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	get current time and vector
		const vector_type& u = *vSol->solution(0);

	//	create local time series
		LocalVectorTimeSeries locTimeSeries;
		locTimeSeries.read_times(vSol);

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(MASS | STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries);
		Eval.set_time_point(0);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local algebra
		LocalIndices ind; LocalVector locU; LocalMatrix locJ;

		EL_PROFILE_BEGIN(Elem_AssembleJacobian);
	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind); locJ.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	read local values of time series
			if(Eval.time_series_needed())
				locTimeSeries.read_values(vSol, ind);

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
			}
			UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot prepare element.");

		//	reset local algebra
			locJ = 0.0;

			//EL_PROFILE_BEGIN(Elem_add_JA);
			//	Assemble JA
			try
			{
				Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locJ *= s_a0;

				Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_STATIONARY);
			}
			UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (A).");
			//EL_PROFILE_END();

		//	Assemble JM
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

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Defect
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Defect in the stationary case. (This version
	 * processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in,out]	d				defect
	 * \param[in]		u				solution
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleDefect( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& d,
					const vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if at least one element exists, else return
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locU, locD, tmpLocD;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind); locD.resize(ind); tmpLocD.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind);
			}
			UG_CATCH_THROW("(stationary) AssembleDefect: Cannot prepare element.");

		//	ANALOG to 'domain_disc_elem()' -  modifies the solution, used
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

		//	reset local algebra
			locD = 0.0;

		//	Assemble A
			try
			{
				Eval.add_def_A_elem(locD, locU, elem, vCornerCoords);
			}
			UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute Defect (A).");

		//	Assemble rhs
			try
			{
				tmpLocD = 0.0;
				Eval.add_rhs_elem(tmpLocD, elem, vCornerCoords);
				locD.scale_append(-1, tmpLocD);

			}
			UG_CATCH_THROW("(stationary) AssembleDefect: Cannot compute Rhs.");

		//	send local to global defect
			try{
				spAssTuner->add_local_vec_to_global(d, locD, dd);
			}
			UG_CATCH_THROW("(stationary) AssembleDefect: Cannot add local vector.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(stationary) AssembleDefect: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Defect
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Defect in the instationary case. (This version
	 * processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
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
	template <typename TElem, typename TIterator>
	static void
	AssembleDefect( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& d,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					const std::vector<number>& vScaleMass,
					const std::vector<number>& vScaleStiff,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

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
		DataEvaluator<domain_type> Eval(MASS | STIFF | RHS | EXPL,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries, &vScaleMass, &vScaleStiff);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locD, tmpLocD;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locD.resize(ind); tmpLocD.resize(ind);

		//	read local values of time series
			locTimeSeries.read_values(vSol, ind);

		//	reset contribution of this element
			locD = 0.0;

		//	loop all time points and assemble them
			for(size_t t = 0; t < vScaleStiff.size(); ++t)
			{
				number scale_stiff = vScaleStiff[t];
				
			//	get local solution at timepoint
				LocalVector& locU = locTimeSeries.solution(t);
				Eval.set_time_point(t);

			//	prepare element
				try
				{
					Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, false);
				}
				UG_CATCH_THROW("(instationary) AssembleDefect: Cannot prepare element.");

			//	Assemble M
				try
				{
					tmpLocD = 0.0;
					Eval.add_def_M_elem(tmpLocD, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locD.scale_append(vScaleMass[t], tmpLocD);
				}
				UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Defect (M).");

			//	Assemble A
				try
				{
					if(scale_stiff != 0.0)
					{
						tmpLocD = 0.0;
						Eval.add_def_A_elem(tmpLocD, locU, elem, vCornerCoords, PT_INSTATIONARY);
						locD.scale_append(scale_stiff, tmpLocD);
					}

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

//					UG_ASSERT(vScaleStiff.size() == 2, "Only one step method supported.");
					const number dt = vSol->time(0)-vSol->time(1);
					locD.scale_append(dt, tmpLocD);
				}


			//	Assemble rhs
				try
				{
					if(scale_stiff != 0.0)
					{
						tmpLocD = 0.0;
						Eval.add_rhs_elem(tmpLocD, elem, vCornerCoords, PT_INSTATIONARY);
						locD.scale_append( - scale_stiff, tmpLocD);
					}

					if(t == 0)
					{
						tmpLocD = 0.0;
						Eval.add_rhs_elem(tmpLocD, elem, vCornerCoords, PT_STATIONARY);
						locD.scale_append( -1.0, tmpLocD);
					}
				}
				UG_CATCH_THROW("(instationary) AssembleDefect: Cannot compute Rhs.");
			}

		//	send local to global defect
			try{
				spAssTuner->add_local_vec_to_global(d, locD, dd);
			}
			UG_CATCH_THROW("(instationary) AssembleDefect: Cannot add local vector.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(instationary) AssembleDefect: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(instationary) AssembleDefect: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (stationary) Linear problem
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Matrix and the global Right-Hand Side
	 * of the Linear problem in the stationary case. (This version processes
	 * elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in,out]	A				Matrix
	 * \param[in,out]	rhs				Right-hand side
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleLinear( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					matrix_type& A,
					vector_type& rhs,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

		UG_DLOG(DID_ELEM_DISC_ASSEMBLE_UTIL, 2, ">>OCT_DISC_DEBUG: " << "elem_disc_assemble_util.h: " << "AssembleLinear(): DataEvaluator(): " << id << std::endl);

	//	prepare loop
		Eval.prepare_elem_loop(id, si);

		UG_DLOG(DID_ELEM_DISC_ASSEMBLE_UTIL, 2, ">>OCT_DISC_DEBUG: " << "elem_disc_assemble_util.h: " << "AssembleLinear(): prepare_elem_loop(): " << id << std::endl);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locRhs; LocalMatrix locA;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locRhs.resize(ind); locA.resize(ind);

		//	prepare element
			try
			{
				UG_DLOG(DID_ELEM_DISC_ASSEMBLE_UTIL, 2, ">>OCT_DISC_DEBUG: " << "elem_disc_assemble_util.h: " << "AssembleLinear(): prepare_elem(): " << id << std::endl);
				for(int i = 0; i < 8; ++i)
					UG_DLOG(DID_ELEM_DISC_ASSEMBLE_UTIL, 2, ">>OCT_DISC_DEBUG: " << "elem_disc_assemble_util.h: " << "AssembleLinear(): prepare_elem(): " << "vCornerCoords " << i << ": " << vCornerCoords[i] << std::endl);

				Eval.prepare_elem(locRhs, elem, id, vCornerCoords, ind, true);
			}
			UG_CATCH_THROW("(stationary) AssembleLinear: Cannot prepare element.");

		//	reset local algebra
			locA = 0.0;
			locRhs = 0.0;

		//	Assemble JA
			try
			{
				Eval.add_jac_A_elem(locA, locRhs, elem, vCornerCoords);
			}
			UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute Jacobian (A).");

		//	Assemble rhs
			try
			{
				Eval.add_rhs_elem(locRhs, elem, vCornerCoords);
			}
			UG_CATCH_THROW("(stationary) AssembleLinear: Cannot compute Rhs.");

		//	send local to global matrix & rhs
			try{
				spAssTuner->add_local_mat_to_global(A, locA, dd);
				spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
			}
			UG_CATCH_THROW("(stationary) AssembleLinear: Cannot add local vector/matrix.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(stationary) AssembleLinear: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Linear problem
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Matrix and the global Right-Hand Side
	 * of the Linear problem in the time-dependent case. (This version processes
	 * elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
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
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleLinear( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					matrix_type& A,
					vector_type& rhs,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					const std::vector<number>& vScaleMass,
					const std::vector<number>& vScaleStiff,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

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
		DataEvaluator<domain_type> Eval(MASS | STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries, &vScaleMass, &vScaleStiff);

	//	prepare loop
		Eval.prepare_elem_loop(id, si);

	//	local algebra
		LocalIndices ind; LocalVector locRhs, tmpLocRhs; LocalMatrix locA, tmpLocA;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
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

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

		if (!spAssTuner->matrix_is_const())
		{
		//	Assemble JM
			try
			{
				tmpLocA = 0.0;
				Eval.add_jac_M_elem(tmpLocA, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locA.scale_append(vScaleMass[0], tmpLocA);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (M).");

		//	Assemble JA
			try
			{
				if (vScaleStiff[0] != 0.0)
				{
					tmpLocA = 0.0;
					Eval.add_jac_A_elem(tmpLocA, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locA.scale_append(vScaleStiff[0], tmpLocA);
				}
				Eval.add_jac_A_elem(locA, locU, elem, vCornerCoords, PT_STATIONARY);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");
		}
		//	Assemble rhs
			try
			{
				if (vScaleStiff[0] != 0.0)
				{
					tmpLocRhs = 0.0;
					Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(vScaleStiff[0], tmpLocRhs);
				}
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
				number scaleStiff = vScaleStiff[t];

			//	prepare element
				try
				{
					Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, false);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

			//	Assemble dM
				try
				{
					tmpLocRhs = 0.0;
					Eval.add_def_M_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(-vScaleMass[t], tmpLocRhs);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (M).");

			//	Assemble dA
				try
				{
					if (scaleStiff != 0.0)
					{
						tmpLocRhs = 0.0;
						Eval.add_def_A_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
						locRhs.scale_append(-scaleStiff, tmpLocRhs);
					}
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

			//	Assemble rhs
				try
				{
					if (scaleStiff != 0.0)
					{
						tmpLocRhs = 0.0;
						Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
						locRhs.scale_append(scaleStiff, tmpLocRhs);
					}
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");
			}

		//	send local to global matrix & rhs
			try{
				if (!spAssTuner->matrix_is_const())
					spAssTuner->add_local_mat_to_global(A, locA, dd);
				spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot add local vector/matrix.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(instationary) AssembleLinear: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble Rhs (of a stationary problem)
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Right-Hand Side. (This version processes
	 * elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in,out]	rhs				Right-hand side
	 * \param[in]		u				solution
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleRhs(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& rhs,
					const vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

	//	prepare loop
		Eval.prepare_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locU, locRhs;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind); locRhs.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind);
			}
			UG_CATCH_THROW("AssembleRhs: Cannot prepare element.");

		//	reset local algebra
			locRhs = 0.0;

		//	Assemble rhs
			try
			{
				Eval.add_rhs_elem(locRhs, elem, vCornerCoords);
			}
			UG_CATCH_THROW("AssembleRhs: Cannot compute Rhs.");

		//	send local to global rhs
			try{
				spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
			}
			UG_CATCH_THROW("AssembleRhs: Cannot add local vector.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("AssembleRhs: Cannot finish element loop.");

		}
		UG_CATCH_THROW("AssembleRhs: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Assemble (instationary) Rhs
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function adds the contributions of all passed element discretizations
	 * on one given subset to the global Right-Hand Side in the time-dependent case.
	 * (This version processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
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
	template <typename TElem, typename TIterator>
	static void
	AssembleRhs(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& rhs,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					const std::vector<number>& vScaleMass,
					const std::vector<number>& vScaleStiff,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

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
		DataEvaluator<domain_type> Eval(MASS | STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries, &vScaleMass, &vScaleStiff);

	//	prepare loop
		Eval.prepare_elem_loop(id, si);

	//	local algebra
		LocalIndices ind; LocalVector locRhs, tmpLocRhs;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
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

		//	prepare element
			try
			{
				Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, false);
			}
			UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element.");

		//	Assemble rhs
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

			//	prepare element
				try
				{
					Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, false);
				}
				UG_CATCH_THROW("(instationary) AssembleRhs: Cannot prepare element.");

			//	Assemble dM
				try
				{
					tmpLocRhs = 0.0;
					Eval.add_def_M_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(-vScaleMass[t], tmpLocRhs);
				}
				UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Jacobian (M).");

			//	Assemble dA
				try
				{
					tmpLocRhs = 0.0;
					Eval.add_def_A_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);
				}
				UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Jacobian (A).");

			//	Assemble rhs
				try
				{
					tmpLocRhs = 0.0;
					Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
				}
				UG_CATCH_THROW("(instationary) AssembleRhs: Cannot compute Rhs.");
			}

			//	send local to global rhs
				try{
					spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
				}
				UG_CATCH_THROW("(instationary) AssembleRhs: Cannot add local vector.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot finish element loop.");

		}
		UG_CATCH_THROW("(instationary) AssembleRhs: Cannot create Data Evaluator.");
	}

// //////////////////////////////////////////////////////////////////////////////
// Prepare Timestep (for instationary problems)
// /////////////////////////////////////////////////////////////////////////////
public:
	/**
	 * This function prepares the global discretization for a time-stepping scheme
	 * by calling the "prepare_timestep" methods of all passed element
	 * discretizations.
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		vSol			current and previous solutions
	 * \param[in]		spAssTuner		assemble adapter
	 */
	static void
	PrepareTimestep
	(
		const std::vector<IElemDisc<domain_type>*>& vElemDisc,
		ConstSmartPtr<DoFDistribution> dd,
		bool bNonRegularGrid,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		number future_time,
		ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner
	)
	{
	//	get current time and vector
		const vector_type& u = *vSol->solution(0);

	//	create local time series
		LocalVectorTimeSeries locTimeSeries;
		locTimeSeries.read_times(vSol);

		try
		{
			DataEvaluator<domain_type> Eval(NONE,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries);
			Eval.set_time_point(0);

		//	prepare timestep
			try
			{
				VectorProxy<vector_type> vp(u);
				size_t algebra_index = bridge::AlgebraTypeIDProvider::instance().id<algebra_type>();
				Eval.prepare_timestep(future_time, vSol->time(0), &vp, algebra_index);
			}
			UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot prepare timestep.");

		}
		UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Prepare Timestep Elem (for instationary problems)
////////////////////////////////////////////////////////////////////////////////
public:
	/**
	 * This function prepares the global discretization for a time-stepping scheme
	 * by calling the "prepare_timestep_elem" methods of all passed element
	 * discretizations on one given subset.
	 * (This version processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		vSol			current and previous solutions
	 * \param[in]		spAssTuner		assemble adapter
	 */
	template <typename TElem, typename TIterator>
	static void
	PrepareTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	get current time and vector
		const vector_type& u = *vSol->solution(0);

	//	create local time series
		LocalVectorTimeSeries locTimeSeries;
		locTimeSeries.read_times(vSol);

		try
		{
		DataEvaluator<domain_type> Eval(NONE,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries);
		Eval.set_time_point(0);

	//	prepare element loop
		Eval.prepare_elem_loop(id, si);

	//	local algebra
		LocalIndices ind; LocalVector locU;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	read local values of time series
			if(Eval.time_series_needed())
				locTimeSeries.read_values(vSol, ind);

		//	prepare timestep
			try
			{
				Eval.prepare_timestep_elem(vSol->time(0), locU, elem, vCornerCoords);
			}
			UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot prepare timestep.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot finish element loop.");
		
		}
		UG_CATCH_THROW("(instationary) PrepareTimestep: Cannot create Data Evaluator.");
	}

// //////////////////////////////////////////////////////////////////////////////
// Finish Timestep (for instationary problems)
// /////////////////////////////////////////////////////////////////////////////
public:
	/**
	 * This function finishes the global discretization for a time-stepping scheme
	 * by calling the "finish_timestep" methods of all passed element
	 * discretizations.
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		vSol			current and previous solutions
	 * \param[in]		spAssTuner		assemble adapter
	 */
	static void
	FinishTimestep
	(
		const std::vector<IElemDisc<domain_type>*>& vElemDisc,
		ConstSmartPtr<DoFDistribution> dd,
		bool bNonRegularGrid,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner
	)
	{
	//	get current time and vector
		const vector_type& u = *vSol->solution(0);

	//	create local time series
		LocalVectorTimeSeries locTimeSeries;
		locTimeSeries.read_times(vSol);

		try
		{
			DataEvaluator<domain_type> Eval(NONE,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries);
			Eval.set_time_point(0);

		//	finish time step
			try
			{
				VectorProxy<vector_type> vp(u);
				size_t algebra_index = bridge::AlgebraTypeIDProvider::instance().id<algebra_type>();
				Eval.finish_timestep(vSol->time(0), &vp, algebra_index);
			}
			UG_CATCH_THROW("(instationary) FinishTimestep: Cannot prepare time step.");

		}
		UG_CATCH_THROW("(instationary) FinishTimestep: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Finish Timestep Elem (for instationary problems)
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function finalizes the global discretization in a time-stepping scheme
	 * by calling the "finish_timestep_elem" methods of all passed element
	 * discretizations on one given subset.
	 * (This version processes elements in a given interval.)
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
	template <typename TElem, typename TIterator>
	static void
	FinishTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				   ConstSmartPtr<domain_type> spDomain,
				   ConstSmartPtr<DoFDistribution> dd,
				   TIterator iterBegin,
				   TIterator iterEnd,
				   int si, bool bNonRegularGrid,
				   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				   ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	get current time and vector
		const vector_type& u = *vSol->solution(0);

	//	create local time series
		LocalVectorTimeSeries locTimeSeries;
		locTimeSeries.read_times(vSol);

	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(NONE,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   &locTimeSeries);
		Eval.set_time_point(0);

	//	prepare loop
		Eval.prepare_elem_loop(id, si);

	//	local algebra
		LocalIndices ind; LocalVector locU;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	check if elem is skipped from assembling
			if(!spAssTuner->element_used(elem)) continue;

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	read local values of time series
			if(Eval.time_series_needed())
				locTimeSeries.read_values(vSol, ind);

		//	finish timestep
			try
			{
				Eval.finish_timestep_elem(locTimeSeries.time(0), locU, elem, vCornerCoords);
			}
			UG_CATCH_THROW("(instationary) FinishTimestepElem: Cannot finish timestep.");
		}

	//	finish element loop
		try
		{
			Eval.finish_elem_loop();
		}
		UG_CATCH_THROW("AssembleMassMatrix: Cannot finish element loop.");
		
		}
		UG_CATCH_THROW("(instationary) FinishTimestepElem: Cannot create Data Evaluator");
	}

////////////////////////////////////////////////////////////////////////////////
// Init. all exports (an optional operation, to use the exports for plotting etc.)
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function finalizes the global discretization in a time-stepping scheme
	 * by calling the "finish_timestep_elem" methods of all passed element
	 * discretizations on one given subset.
	 * (This version processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		bAsTimeDependent flag to simulate the instationary case for the initialization
	 */
	template <typename TElem, typename TIterator>
	static void
	InitAllExports(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				   ConstSmartPtr<DoFDistribution> dd,
				   TIterator iterBegin,
				   TIterator iterEnd,
				   int si, bool bNonRegularGrid, bool bAsTimeDependent)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	dummy local time series (only to simulate the instationary case for the initialization)
		LocalVectorTimeSeries locTimeSeries;
		
	//	prepare for given elem discs
		try
		{
		DataEvaluator<domain_type> Eval(MASS | STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid,
						   bAsTimeDependent? &locTimeSeries : nullptr);
		Eval.set_time_point(0);

	//	prepare loop
		Eval.prepare_elem_loop(id, si);


	//	finish element loop
		Eval.finish_elem_loop();
		
		}
		UG_CATCH_THROW("InitAllExports: Data Evaluator failed");
	}

////////////////////////////////////////////////////////////////////////////////
// Error estimator (stationary)
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function assembles the error estimator associated with all the
	 * element discretizations in the internal data structure for one given subset.
	 * (This version processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		u				solution
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleErrorEstimator
	(
		const std::vector<IElemError<domain_type>*>& vElemDisc,
		ConstSmartPtr<domain_type> spDomain,
		ConstSmartPtr<DoFDistribution> dd,
		TIterator iterBegin,
		TIterator iterEnd,
		int si,
		bool bNonRegularGrid,
		const vector_type& u
	)
	{
	//	check if there are any elements at all, otherwise return immediately
		if(iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

	//	prepare for given elem discs
		try
		{
		ErrorEvaluator<domain_type> Eval(STIFF | RHS,
						   vElemDisc, dd->function_pattern(), bNonRegularGrid);

	//	prepare element loop
		Eval.prepare_err_est_elem_loop(id, si);

	//	local indices and local algebra
		LocalIndices ind; LocalVector locU;

	//	Loop over all elements
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get Element
			TElem* elem = *iter;

		//	get corner coordinates
			FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

		//	get global indices
			dd->indices(elem, ind, Eval.use_hanging());

		//	adapt local algebra
			locU.resize(ind);

		//	read local values of u
			GetLocalVector(locU, u);

		//	prepare element
			try
			{
				Eval.prepare_err_est_elem(locU, elem, vCornerCoords, ind, false);
			}
			UG_CATCH_THROW("(stationary) AssembleRhs: Cannot prepare element.");

		//	assemble stiffness part
			try
			{
				Eval.compute_err_est_A_elem(locU, elem, vCornerCoords, ind);
			}
			UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

		//	assemble rhs part
			try
			{
				Eval.compute_err_est_rhs_elem(elem, vCornerCoords, ind);
			}
			UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

		}

	//	finish element loop
		try
		{
			Eval.finish_err_est_elem_loop();
		}
		UG_CATCH_THROW("AssembleErrorEstimator: Cannot finish element loop.");

		}
		UG_CATCH_THROW("AssembleErrorEstimator: Cannot create Data Evaluator.");
	}

////////////////////////////////////////////////////////////////////////////////
// Error estimator (for time-dependent problems)
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function assembles the error estimator associated with all the
	 * element discretizations in the internal data structure for one given subset.
	 * (This version processes elements in a given interval.)
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		spDomain		domain
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		iterBegin		element iterator
	 * \param[in]		iterEnd			element iterator
	 * \param[in]		si				subset index
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		vScaleMass		scaling factors for mass part
	 * \param[in]		vScaleStiff		scaling factors for stiffness part
	 * \param[in]		vSol				solution
	 */
	template <typename TElem, typename TIterator>
	static void
	AssembleErrorEstimator
	(
		const std::vector<IElemError<domain_type>*>& vElemDisc,
		ConstSmartPtr<domain_type> spDomain,
		ConstSmartPtr<DoFDistribution> dd,
		TIterator iterBegin,
		TIterator iterEnd,
		int si,
		bool bNonRegularGrid,
		const std::vector<number>& vScaleMass,
		const std::vector<number>& vScaleStiff,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol
	)
	{
	//	check if there are any elements at all, otherwise return immediately
		if (iterBegin == iterEnd) return;

	//	reference object id
		static constexpr ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	//	storage for corner coordinates
		MathVector<domain_type::dim> vCornerCoords[TElem::NUM_VERTICES];

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
			ErrorEvaluator<domain_type> Eval(MASS | STIFF | RHS,
							   vElemDisc, dd->function_pattern(), bNonRegularGrid,
							   &locTimeSeries, &vScaleMass, &vScaleStiff);

		//	prepare element loop
			Eval.prepare_err_est_elem_loop(id, si);

		//	local indices and local algebra
			LocalIndices ind;

		//	loop over all elements
			for (TIterator iter = iterBegin; iter != iterEnd; ++iter)
			{
			//	get Element
				TElem* elem = *iter;

			//	get corner coordinates
				FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

			//	get global indices
				dd->indices(elem, ind, Eval.use_hanging());

			//	read local values of time series
				locTimeSeries.read_values(vSol, ind);

			//	loop all time points and assemble them
				for (std::size_t t = 0; t < vScaleStiff.size(); ++t)
				{
				//	get local solution at timepoint
					LocalVector& locU = locTimeSeries.solution(t);
					Eval.set_time_point(t);

				//	prepare element
					try
					{
						Eval.prepare_err_est_elem(locU, elem, vCornerCoords, ind, false);
					}
					UG_CATCH_THROW("AssembleErrorEstimator: Cannot prepare element.");

				//	assemble stiffness part
					try
					{
						Eval.compute_err_est_A_elem(locU, elem, vCornerCoords, ind, vScaleMass[t], vScaleStiff[t]);
					}
					UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

				//	assemble mass part
					try
					{
						Eval.compute_err_est_M_elem(locU, elem, vCornerCoords, ind, vScaleMass[t], vScaleStiff[t]);
					}
					UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");

				//	assemble rhs part
					try
					{
						Eval.compute_err_est_rhs_elem(elem, vCornerCoords, ind, vScaleMass[t], vScaleStiff[t]);
					}
					UG_CATCH_THROW("AssembleErrorEstimator: Cannot assemble the error estimator for stiffness part.");
				}
			}

			//	finish element loop
			try
			{
				Eval.finish_err_est_elem_loop();
			}
			UG_CATCH_THROW("AssembleErrorEstimator: Cannot finish element loop.");

		}
		UG_CATCH_THROW("AssembleErrorEstimator: Cannot create Data Evaluator.");
	}

}; // class StdGlobAssembler

} // end namespace ug


#endif