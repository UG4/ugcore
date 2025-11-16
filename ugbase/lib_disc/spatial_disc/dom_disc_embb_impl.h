/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/*
 * Implementation of functions from embass.h
 */

namespace ug {

/*-------- class LSGFGlobAssembler --------*/

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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, 0);
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locU.resize(ind); locA.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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
				spAssTuner->add_local_mat_to_global(A, locA, dd);
			}
			UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot add local matrix.");
		}
		else
		{
		//	assemble the matrix for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
			//	read local values of u
				GetLocalVector(locU, u);
				
			//	extrapolate the inner values to the outer ones
				m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

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
				
			//	eliminate the outer connections in the local matrix
				m_extrapol.template eliminate_extrapolated<TElem> (locA, base_co);

			//	send local to global matrix
				try{
					spAssTuner->add_local_mat_to_global(A, locA, dd);
				}
				UG_CATCH_THROW("AssembleStiffnessMatrix: Cannot add local matrix.");
			}
		}
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;
		
	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, 0);
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());
		
	//	adapt local algebra
		locU.resize(ind); locM.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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
		else // an interface element
		{
		//	assemble the matrix for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
			//	read local values of u
				GetLocalVector(locU, u);

			//	extrapolate the inner values to the outer ones
				m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

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

			//	eliminate the outer connections in the local matrix
				m_extrapol.template eliminate_extrapolated<TElem> (locM, base_co);

			// send local to global matrix
				try{
					spAssTuner->add_local_mat_to_global(M, locM, dd);
				}
				UG_CATCH_THROW("AssembleMassMatrix: Cannot add local matrix.");
			}
		}
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, 0);
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locU.resize(ind); locJ.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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
		else // an interface element
		{
		//	assemble the jacobian for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
			//	read local values of u
				GetLocalVector(locU, u);

			//	extrapolate the inner values to the outer ones
				m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

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

			//	eliminate the outer connections in the local matrix
				m_extrapol.template eliminate_extrapolated<TElem> (locJ, base_co);

			// send local to global matrix
				try{
					spAssTuner->add_local_mat_to_global(J, locJ, dd);
				}
				UG_CATCH_THROW("(stationary) AssembleJacobian: Cannot add local matrix.");
			}
		}
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, vSol->time(0));
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locU.resize(ind); locJ.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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

		//	Assemble JA
			try
			{
				Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locJ *= s_a0;

				Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_STATIONARY);
			}
			UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (A).");
			
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
		else // an interface element
		{
		//	assemble the jacobian for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
			//	read local values of u
				GetLocalVector(locU, u);

			//	extrapolate the inner values to the outer ones
				m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

			//	read local values of time series
				if(Eval.time_series_needed())
				{
					locTimeSeries.read_values(vSol, ind);
					for (size_t i = 0; i < locTimeSeries.size (); i++)
						m_extrapol.template extrapolate_sol_by_lsf<TElem>
							(locTimeSeries.solution (i), base_co);
				}

			//	prepare element
				try
				{
					Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
				}
				UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot prepare element.");

			//	reset local algebra
				locJ = 0.0;

				//	Assemble JA
				try
				{
					Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locJ *= s_a0;

					Eval.add_jac_A_elem(locJ, locU, elem, vCornerCoords, PT_STATIONARY);
				}
				UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (A).");

			//	Assemble JM
				try
				{
					Eval.add_jac_M_elem(locJ, locU, elem, vCornerCoords, PT_INSTATIONARY);
				}
				UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot compute Jacobian (M).");

			//	eliminate the outer connections in the local matrix
				m_extrapol.template eliminate_extrapolated<TElem> (locJ, base_co);

			// send local to global matrix
				try{
					spAssTuner->add_local_mat_to_global(J, locJ, dd);
				}
				UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot add local matrix.");
			}
		}
	}

//	finish element loop
	try
	{
		Eval.finish_elem_loop();
	}
	UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot finish element loop.");

	}
	UG_CATCH_THROW("(instationary) AssembleJacobian: Cannot create Data Evaluator.");
}

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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, 0);
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locU.resize(ind); locD.resize(ind); tmpLocD.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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
		else // an interface element
		{
		//	assemble the defect for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
			//	read local values of u
				GetLocalVector(locU, u);

			//	extrapolate the inner values to the outer ones
				m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

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

			//	eliminate the outer connections in the local matrix
				m_extrapol.template clear_outer_vectors<TElem> (locD, base_co);

			//	send local to global defect
				try{
					spAssTuner->add_local_vec_to_global(d, locD, dd);
				}
				UG_CATCH_THROW("(stationary) AssembleDefect: Cannot add local vector.");
			}
		}
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, vSol->time(0));
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locD.resize(ind); tmpLocD.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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
		else // an interface element
		{
		//	assemble the defect for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
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

				//	extrapolate the inner values to the outer ones
					m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

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

//						UG_ASSERT(vScaleStiff.size() == 2, "Only one step method supported.");
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

			//	eliminate the outer connections in the local matrix
				m_extrapol.template clear_outer_vectors<TElem>(locD, base_co);

			//	send local to global defect
				try{
					spAssTuner->add_local_vec_to_global(d, locD, dd);
				}
				UG_CATCH_THROW("(instationary) AssembleDefect: Cannot add local vector.");
			}
		}
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	prepare loop
	Eval.prepare_elem_loop(id, si);

//	local indices and local algebra
	LocalIndices ind; LocalVector locRhs; LocalMatrix locA;

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, 0);
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locRhs.resize(ind); locA.resize(ind);

	//	prepare element
		try
		{
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
		
	//	For interface elements, correct the matrix and the right-hand side
		if (elem_status == 0)
			m_extrapol.template eliminate_extrapolated<TElem> (locA, locRhs);

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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, vSol->time(0));
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

	//	get global indices
		dd->indices(elem, ind, Eval.use_hanging());

	//	adapt local algebra
		locRhs.resize(ind); tmpLocRhs.resize(ind);
		locA.resize(ind); tmpLocA.resize(ind);

	//	assemble depending on whether this is an inner or an interface element
		if (elem_status > 0) // an inner element, usual assembling
		{
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
				tmpLocA = 0.0;
				Eval.add_jac_A_elem(tmpLocA, locU, elem, vCornerCoords, PT_INSTATIONARY);
				locA.scale_append(vScaleStiff[0], tmpLocA);

				Eval.add_jac_A_elem(locA, locU, elem, vCornerCoords, PT_STATIONARY);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

		//	Assemble rhs
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
					tmpLocRhs = 0.0;
					Eval.add_def_A_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

			//	Assemble rhs
				try
				{
					tmpLocRhs = 0.0;
					Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
					locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");
			}

		//	send local to global matrix & rhs
			try{
				spAssTuner->add_local_mat_to_global(A, locA, dd);
				spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
			}
			UG_CATCH_THROW("(instationary) AssembleLinear: Cannot add local vector/matrix.");
		}
		else // an interface element
		{
		//	assemble the system for every corner separately
			for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
			if (m_extrapol.corner_inside (base_co)) /* ... only for the corners 'inside' */
			{
			//	read local values of time series
				locTimeSeries.read_values(vSol, ind);
				Eval.set_time_point(0);

			//	reset element contribution
				locA = 0.0; locRhs = 0.0;

			/////////////////////
			//	current time step

			//	get local solution at time point
				LocalVector& locU = locTimeSeries.solution(0);

			//	extrapolate the inner values to the outer ones
				m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

			//	prepare element
				try
				{
					Eval.prepare_elem(locU, elem, id, vCornerCoords, ind, true);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot prepare element.");

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
					tmpLocA = 0.0;
					Eval.add_jac_A_elem(tmpLocA, locU, elem, vCornerCoords, PT_INSTATIONARY);
					locA.scale_append(vScaleStiff[0], tmpLocA);

					Eval.add_jac_A_elem(locA, locU, elem, vCornerCoords, PT_STATIONARY);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

			//	Assemble rhs
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

				//	extrapolate the inner values to the outer ones
					m_extrapol.template extrapolate_sol_by_lsf<TElem> (locU, base_co);

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
						tmpLocRhs = 0.0;
						Eval.add_def_A_elem(tmpLocRhs, locU, elem, vCornerCoords, PT_INSTATIONARY);
						locRhs.scale_append(-vScaleStiff[t], tmpLocRhs);
					}
					UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Jacobian (A).");

				//	Assemble rhs
					try
					{
						tmpLocRhs = 0.0;
						Eval.add_rhs_elem(tmpLocRhs, elem, vCornerCoords, PT_INSTATIONARY);
						locRhs.scale_append(vScaleStiff[t], tmpLocRhs);
					}
					UG_CATCH_THROW("(instationary) AssembleLinear: Cannot compute Rhs.");
				}

			//	eliminate the outer connections in the local matrix
				m_extrapol.template eliminate_extrapolated<TElem> (locA, locRhs, base_co);

			//	send local to global matrix & rhs
				try{
					spAssTuner->add_local_mat_to_global(A, locA, dd);
					spAssTuner->add_local_vec_to_global(rhs, locRhs, dd);
				}
				UG_CATCH_THROW("(instationary) AssembleLinear: Cannot add local vector/matrix.");
			}
		}
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

template <typename TDomain, typename TAlgebra, typename TExtrapolation>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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


/**
 * This function prepares the global discretization for a time-stepping scheme
 * by calling the "prepare_timestep_elem" methods of all passed element
 * discretizations on one given subset.
 * (This version processes elements in a given interval.)
 *
 * REMARK: This version does not extrapolate anything. It merely skips outer elements.
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
PrepareTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<domain_type> spDomain,
				ConstSmartPtr<DoFDistribution> dd,
				TIterator iterBegin,
				TIterator iterEnd,
				int si, bool bNonRegularGrid,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, vSol->time(0));
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

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
		UG_CATCH_THROW("(instationary) PrepareTimestepElem: Cannot prepare timestep.");
	}

	}
	UG_CATCH_THROW("(instationary) PrepareTimestepElem: Cannot create Data Evaluator.");
}

template <typename TDomain, typename TAlgebra, typename TExtrapolation>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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

/**
 * This function finalizes the global discretization in a time-stepping scheme
 * by calling the "finish_timestep_elem" methods of all passed element
 * discretizations on one given subset.
 * (This version processes elements in a given interval.)
 *
 * REMARK: This version does not extrapolate anything. It merely skips outer elements.
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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
FinishTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
			   ConstSmartPtr<domain_type> spDomain,
			   ConstSmartPtr<DoFDistribution> dd,
			   TIterator iterBegin,
			   TIterator iterEnd,
			   int si, bool bNonRegularGrid,
			   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			   ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
{
//	check the dof distribution for the ghost-fluid method
	m_extrapol.check_dd (dd);

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

//	grid level of all the elements
	int g_level = dd->grid_level().level ();
	if (g_level == GridLevel::TOP)
		g_level = dd->multi_grid()->top_level ();
	
//	Loop over all elements
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(!spAssTuner->element_used(elem)) continue;

	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *spDomain);

	//	check whether we are inside
		int elem_status = m_extrapol.check_elem_lsf
			(TElem::NUM_VERTICES, elem, si, g_level, Eval.use_hanging(), vCornerCoords, vSol->time(0));
		if (elem_status < 0) continue; // this is an outer element, do not assemble it at all
		if (m_bAssembleOnlyCut && elem_status > 0) continue; // exclude elements that are not cut

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

	}
	UG_CATCH_THROW("(instationary) FinishTimestepElem: Cannot create Data Evaluator");
}

////////////////////////////////////////////////////////////////////////////////
// Init. all exports (an optional operation, to use the exports for plotting etc.)
////////////////////////////////////////////////////////////////////////////////

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
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
template <typename TElem, typename TIterator>
void LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation>::
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

} // end namespace ug

/* End of File */
