/*
 * assemble_elem_disc.h
 *
 *  Created on: 08.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// intern headers
#include "../../reference_element/reference_element.h"
#include "./elem_disc_interface.h"
#include "lib_discretization/common/function_group.h"
#include "lib_discretization/spatial_discretization/ip_data/data_evaluator.h"

//#define PROFILE_ELEM_LOOP
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

//////////////////////////////////
// Assemble Stiffness Matrix
//////////////////////////////////

/**
 * This function adds to the Stiffness matrix the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Stiffness matrix
 * \param[in]		u				solution
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleStiffnessMatrix(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
                        	const IDoFDistribution<TDoFDistribution>& dofDistr,
                        	int si, bool bNonRegularGrid,
                        	typename TAlgebra::matrix_type& A,
                        	const typename TAlgebra::vector_type& u,
                        	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> locA;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time-independent
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		vElemDisc[i]->set_time_dependent(false);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleStiffnessMatrix': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locA.set_indices(ind, ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locA.set_indices(ind, ind);
		}

	// 	read local values of u
		locU.read_values(u);

	// 	reset local matrix and rhs
		locA.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleStiffnessMatrix': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, true))
		{
			UG_LOG("ERROR in 'AssembleStiffnessMatrix': Cannot compute elem data.\n");
			return false;
		}

	// 	Assemble JA
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_JA(locA, locU))
			{
				UG_LOG("ERROR in 'AssembleStiffnessMatrix': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	//	Compute element data
		Eval.add_coupl_JA(locA, ind);

	// 	send local to global matrix
		A.add(locA);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleStiffnessMatrix: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}


//////////////////////////////////
// Assemble Mass Matrix
//////////////////////////////////
/**
 * This function adds to the Mass matrix the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	M				Mass matrix
 * \param[in]		u				solution
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleMassMatrix(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
					const IDoFDistribution<TDoFDistribution>& dofDistr,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& M,
					const typename TAlgebra::vector_type& u,
                	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> locM;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time-independent
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		vElemDisc[i]->set_time_dependent(false);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleMassMatrix': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locM.set_indices(ind, ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN))
	{
		UG_LOG("ERROR in 'AssembleMassMatrix': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locM.set_indices(ind, ind);
		}

	// 	read local values of u
		locU.read_values(u);

	// 	reset local matrix and rhs
		locM.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleMassMatrix': Cannot prepare element.\n");
				return false;
			}
		}

	// 	Assemble JM
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JM
			if(!vElemDisc[i]->assemble_JM(locM, locU))
			{
				UG_LOG("ERROR in 'AssembleMassMatrix': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	// 	send local to global matrix
		M.add(locM);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleMassMatrix: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}

//////////////////////////////////
// Assemble Jacobian
//////////////////////////////////

/**
 * This function adds to the jacobian the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		u				solution
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleJacobian(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
					const IDoFDistribution<TDoFDistribution>& dofDistr,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u,
                	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> locJ;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare indices.\n");
		return false;
	}

//	set time-independent
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		vElemDisc[i]->set_time_dependent(false);

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locJ.set_indices(ind, ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locJ.set_indices(ind, ind);
		}

	// 	read local values of u
		locU.read_values(u);

	// 	reset local matrix and rhs
		locJ.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, true))
		{
			UG_LOG("ERROR in 'AssembleJacobian': Cannot compute elem data.\n");
			return false;
		}

	//	Compute element data
		if(!Eval.compute_lin_defect_JA(locU, ind))
		{
			UG_LOG("ERROR in 'AssembleJacobian': Cannot compute lin_defect_JA.\n");
			return false;
		}

	// 	Assemble JA
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_JA(locJ, locU))
			{
				UG_LOG("ERROR in 'AssembleJacobian': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	//	Compute element data
		Eval.add_coupl_JA(locJ, ind);


	// 	send local to global matrix
		J.add(locJ);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}

/**
 * This function adds to the jacobian the entries of one subset for all passed
 * element discretizations.
 * Note, that it is assumed	to have s_m0 == 1
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		u				solution
 * \param[in]		s_a0			scaling factor for stiffness part
 * \param[in]		time			current time
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleJacobian(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
					const IDoFDistribution<TDoFDistribution>& dofDistr,
					int si, bool bNonRegularGrid,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u, number time,
	                const SolutionTimeSeries<typename TAlgebra::vector_type>& solList,
					number s_a0,
                	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> locJ;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time dependent part
	LocalVectorTimeSeries<typename TAlgebra::vector_type> locTimeSeries(solList);
	bool bNeedLocTimeSeries = false;
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		bNeedLocTimeSeries |= vElemDisc[i]->set_time_dependent(true, time, &locTimeSeries);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locJ.set_indices(ind, ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN, time))
	{
		UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locJ.set_indices(ind, ind);
		}

	// 	read local values of u
		locU.read_values(u);

	//	read local values of time series
		if(bNeedLocTimeSeries) locTimeSeries.read_values(ind);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleJacobian': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, true))
		{
			UG_LOG("ERROR in 'AssembleJacobian': Cannot compute elem data.\n");
			return false;
		}

	//	Compute element data
		if(!Eval.compute_lin_defect_JA(locU, ind))
		{
			UG_LOG("ERROR in 'AssembleJacobian': Cannot compute lin_defect_JA.\n");
			return false;
		}

	//	Compute element data
		if(!Eval.compute_lin_defect_JM(locU, ind))
		{
			UG_LOG("ERROR in 'AssembleJacobian': Cannot compute lin_defect_JM.\n");
			return false;
		}

	// 	reset local matrix
		locJ.set(0.0);

	// 	Assemble JA
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_JA(locJ, locU))
			{
				UG_LOG("ERROR in 'AssembleJacobian': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	//	Compute element data
		Eval.add_coupl_JA(locJ, ind);

	//	scale stiffness part
		locJ *= s_a0;

	// 	Assemble JM
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JM
			if(!vElemDisc[i]->assemble_JM(locJ, locU))
			{
				UG_LOG("ERROR in 'AssembleJacobian': Cannot assemble local "
						"Mass Matrix.\n");
				return false;
			}
		}

	//	Compute element data
		Eval.add_coupl_JM(locJ, ind);

	// 	send local to global matrix
		J.add(locJ);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}

//////////////////////////////////
// Assemble Defect
//////////////////////////////////

/**
 * This function adds to the defect the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		u				solution
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleDefect(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
               	const IDoFDistribution<TDoFDistribution>& dofDistr,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
               	const typename TAlgebra::vector_type& u,
            	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleDefect': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalVector<typename TAlgebra::vector_type::value_type> locD;
	LocalVector<typename TAlgebra::vector_type::value_type> locRhs;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time-independent
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		vElemDisc[i]->set_time_dependent(false);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleDefect': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locD.set_indices(ind);
		locRhs.set_indices(ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_DEFECT))
	{
		UG_LOG("ERROR in 'AssembleDefect': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locD.set_indices(ind);
			locRhs.set_indices(ind);
		}

	// 	read local values of u
		locU.read_values(u);

	// 	reset local defect and rhs
		locD.set(0.0);
		locRhs.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, false))
		{
			UG_LOG("ERROR in 'AssembleDefect': Cannot compute elem data.\n");
			return false;
		}

	// 	Assemble defect A
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_A(locD, locU))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}


	// 	Assemble rhs
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_f(locRhs))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locD -= locRhs;

	// 	send local to global matrix
		d.add(locD);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleDefect: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}


/**
 * This function adds to the defect the entries of one subset for all passed
 * element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		u				solution
 * \param[in]		s_m				scaling factor for mass part
 * \param[in]		s_a				scaling factor for stiffness part
 * \param[in]		time			current time
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleDefect(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
               	const IDoFDistribution<TDoFDistribution>& dofDistr,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& d,
               	const typename TAlgebra::vector_type& u, number time,
                const SolutionTimeSeries<typename TAlgebra::vector_type>& solList,
               	number s_m, number s_a,
            	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleDefect': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalVector<typename TAlgebra::vector_type::value_type> locD;
	LocalVector<typename TAlgebra::vector_type::value_type> tmpLocD;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleDefect': Cannot prepare indices.\n");
		return false;
	}

//	set time dependent part
	LocalVectorTimeSeries<typename TAlgebra::vector_type> locTimeSeries(solList);
	bool bNeedLocTimeSeries = false;
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		bNeedLocTimeSeries |= vElemDisc[i]->set_time_dependent(true, time, &locTimeSeries);

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locD.set_indices(ind);
		tmpLocD.set_indices(ind);

	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_DEFECT, time))
	{
		UG_LOG("ERROR in 'AssembleDefect': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locD.set_indices(ind);
			tmpLocD.set_indices(ind);
		}

	// 	read local values of u
		locU.read_values(u);

	//	read local values of time series
		if(bNeedLocTimeSeries) locTimeSeries.read_values(ind);

	// 	reset local matrix and rhs
		locD.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, false))
		{
			UG_LOG("ERROR in 'AssembleDefect': Cannot compute elem data.\n");
			return false;
		}

	// 	Assemble defect M
		tmpLocD.set(0.0);
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_M(tmpLocD, locU))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locD += tmpLocD * s_m;

	// 	Assemble defect A
		tmpLocD.set(0.0);
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_A(tmpLocD, locU))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locD += tmpLocD * s_a;

	// 	Assemble defect rhs
		tmpLocD.set(0.0);
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_f(tmpLocD))
			{
				UG_LOG("ERROR in 'AssembleDefect': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locD -= tmpLocD * s_a;

	// 	send local to global matrix
		d.add(locD);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleDefect: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}


//////////////////////////////////
// Assemble Linear
//////////////////////////////////

/**
 * This function adds to the Matrix and to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		u				solution
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleLinear(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
               	const IDoFDistribution<TDoFDistribution>& dofDistr,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
               	const typename TAlgebra::vector_type& u,
            	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalVector<typename TAlgebra::vector_type::value_type> locRhs;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> locA;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time-independent
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		vElemDisc[i]->set_time_dependent(false);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locRhs.set_indices(ind);
		locA.set_indices(ind,ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locRhs.set_indices(ind);
			locA.set_indices(ind,ind);
		}

	// 	read local values of u
		locU.read_values(u);

	// 	reset local matrix and rhs
		locRhs.set(0.0);
		locA.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, false))
		{
			UG_LOG("ERROR in 'AssembleLinear': Cannot compute elem data.\n");
			return false;
		}

	// 	Assemble JA
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_JA(locA, locU))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	// 	Assemble rhs
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble rhs
			if(!vElemDisc[i]->assemble_f(locRhs))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	// 	send local to global matrix
		A.add(locA);

	// 	send local to global rhs
		rhs.add(locRhs);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleLinear: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}

/**
 * This function adds to the Matrix and to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		u				solution
 * \param[in]		s_m				scaling factor for mass part
 * \param[in]		s_a				scaling factor for stiffness part
 * \param[in]		time			current time
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleLinear(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
               	const IDoFDistribution<TDoFDistribution>& dofDistr,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::matrix_type& A,
               	typename TAlgebra::vector_type& rhs,
               	const typename TAlgebra::vector_type& u, number time,
                const SolutionTimeSeries<typename TAlgebra::vector_type>& solList,
               	number s_m, number s_a,
            	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalVector<typename TAlgebra::vector_type::value_type> locRhs;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> locA;
	LocalMatrix<typename TAlgebra::matrix_type::value_type> tmpLocA;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time dependent part
	LocalVectorTimeSeries<typename TAlgebra::vector_type> locTimeSeries(solList);
	bool bNeedLocTimeSeries = false;
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		bNeedLocTimeSeries |= vElemDisc[i]->set_time_dependent(true, time, &locTimeSeries);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locRhs.set_indices(ind);
		locA.set_indices(ind,ind);
		tmpLocA.set_indices(ind,ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN, time))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locRhs.set_indices(ind);
			locA.set_indices(ind,ind);
			tmpLocA.set_indices(ind,ind);
		}

	// 	read local values of u
		locU.read_values(u);

	//	read local values of time series
		if(bNeedLocTimeSeries) locTimeSeries.read_values(ind);

	// 	reset local matrix and rhs
		locRhs.set(0.0);
		locA.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, false))
		{
			UG_LOG("ERROR in 'AssembleLinear': Cannot compute elem data.\n");
			return false;
		}

	// 	Assemble JM
		tmpLocA.set(0.0);
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JM
			if(!vElemDisc[i]->assemble_JM(tmpLocA, locU))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locA += tmpLocA * s_m;

	// 	Assemble JA
		tmpLocA.set(0.0);
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble JA
			if(!vElemDisc[i]->assemble_JA(tmpLocA, locU))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locA += tmpLocA * s_a;

	// 	Assemble rhs
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble rhs
			if(!vElemDisc[i]->assemble_f(locRhs))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}
		locRhs *= s_a;

	// 	send local to global matrix
		A.add(locA);

	// 	send local to global rhs
		rhs.add(locRhs);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleLinear: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}


/**
 * This function adds to the Rhs the entries of one subset
 * for all passed element discretizations.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dofDistr		DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		u				solution
 */
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleRhs(	const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
               	const IDoFDistribution<TDoFDistribution>& dofDistr,
               	int si, bool bNonRegularGrid,
               	typename TAlgebra::vector_type& rhs,
               	const typename TAlgebra::vector_type& u,
            	ISelector* sel = NULL)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	static const ReferenceObjectID refID =
				reference_element_type::REFERENCE_OBJECT_ID;

// 	check if at least on element exist, else return
	if(dofDistr.template num<TElem>(si) == 0) return true;

//	get element iterator
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = dofDistr.template begin<TElem>(si);
	iterEnd = dofDistr.template end<TElem>(si);

//	create data evaluator
	DataEvaluator<TAlgebra> Eval;

//	prepare for given elem discs
	if(!Eval.set_elem_discs(vElemDisc))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot init Evaluation of "
				"Elem Discs. Aborting\n");
		return false;
	}

// 	flag, iff use hanging nodes as well
	bool useHanging = false;

//	check if hanging nodes are needed
	if(bNonRegularGrid)
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use non-regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(true))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support non-regular"
				       " grids, but this is requested.\n");
				return false;
			}

		//	check if hanging dofs are really used
			useHanging |= vElemDisc[i]->use_hanging();
		}
	}
	else
	{
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//  let disc use regular grid assemblings
			if(!vElemDisc[i]->treat_non_regular_grid(false))
			{
				UG_LOG("ERROR: Elem Disc " << i << " does not support regular"
					   " grids, but this is requested.\n");
				return false;
			}
		}
	}

// 	local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::value_type> locU;
	LocalVector<typename TAlgebra::vector_type::value_type> locRhs;

// 	set common functions
	ind.set_function_group(Eval.fct_group());

//	set time-independent
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		vElemDisc[i]->set_time_dependent(false);

// 	prepare local indices for elem type
	if(!dofDistr.template prepare_indices<TElem>(si, ind, useHanging))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot prepare indices.\n");
		return false;
	}

// 	adjust local algebra
	if(!useHanging)
	{
		locU.set_indices(ind);
		locRhs.set_indices(ind);
	}

//	prepare element discs
	if(!Eval.prepare_elem_loop(refID, ind, IEDN_JACOBIAN))
	{
		UG_LOG("ERROR in 'AssembleLinear': Cannot prepare Elem Discs.\n");
		return false;
	}

// 	Loop over all elements
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get Element
		TElem* elem = *iter;

	//	check if elem is skipped from assembling
		if(sel)
			if(!sel->is_selected(elem))
				continue;

	// 	get global indices
		ind.access_all();
		dofDistr.update_indices(elem, ind, useHanging);

	// 	adapt local algebra
		if(useHanging)
		{
			locU.set_indices(ind);
			locRhs.set_indices(ind);
		}

	// 	read local values of u
		locU.read_values(u);

	// 	reset local matrix and rhs
		locRhs.set(0.0);

	// 	prepare element
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	prepare for elem disc
			if(!vElemDisc[i]->prepare_element(elem, locU, ind))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot prepare element.\n");
				return false;
			}
		}

	//	Compute element data
		if(!Eval.compute_elem_data(locU, ind, false))
		{
			UG_LOG("ERROR in 'AssembleLinear': Cannot compute elem data.\n");
			return false;
		}

	// 	Assemble rhs
		for(size_t i = 0; i < vElemDisc.size(); ++i)
		{
		//	access disc functions
			ind.access_by_map(Eval.map(i));

		//	assemble rhs
			if(!vElemDisc[i]->assemble_f(locRhs))
			{
				UG_LOG("ERROR in 'AssembleLinear': Cannot assemble local "
						"Stiffness Matrix.\n");
				return false;
			}
		}

	// 	send local to global rhs
		rhs.add(locRhs);
	}

// finish element loop
	for(size_t i = 0; i < vElemDisc.size(); ++i)
		if(!vElemDisc[i]->finish_element_loop())
		{
			UG_LOG("ERROR in AssembleLinear: Cannot finish element loop.\n");
			return false;
		}

//	we're done
	return true;
}


} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__ */
