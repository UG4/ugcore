/*
 * Copyright (c) 2011-2017:  G-CSC, Goethe University Frankfurt
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

#include <sstream>

#include "data_evaluator.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

DebugID DID_DATA_EVALUATOR("DATA_EVALUATOR");

///////////////////////////////////////////////////////////////////////////////
// prepare / finish
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void DataEvaluator<TDomain>::
prepare_elem_loop(const ReferenceObjectID id, int si)
{
// 	prepare loop (elem disc set local ip series here)
	try{
		for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_prep_elem_loop(id, si);
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_elem_loop: "
						"Cannot prepare element loop.");

//	extract data imports and userdatas
	try{
		extract_imports_and_userdata(si, m_discPart);
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_elem_loop: "
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
	UG_CATCH_THROW("DataEvaluatorBase::prepare_elem_loop: Import correctly implemented.");

//	prepare and check dependent data
	try{
		for(size_t i = 0; i < m_vDependentData.size(); ++i){
			m_vDependentData[i]->check_setup();
		}
	}
	UG_CATCH_THROW("DataEvaluatorBase::prepare_elem_loop: Dependent UserData "
				   " (e.g. Linker or Export) is not ready for evaluation.");

//	evaluate constant data
	for(size_t i = 0; i < m_vConstData.size(); ++i)
		m_vConstData[i]->compute((LocalVector*)nullptr, nullptr, nullptr, false);
}

template <typename TDomain>
void DataEvaluator<TDomain>::finish_elem_loop()
{
//	finish each elem disc
	try{
		for(size_t d = 0; d < m_vElemDisc[PT_ALL].size(); ++d)
		{
			IElemDisc<TDomain>* disc = m_vElemDisc[PT_ALL][d];
			
			disc->do_fsh_elem_loop();
			
			/* TODO:
			 * In prepare_elem_loop, the elemdiscs initialize the local ip's independently
			 * on if they are used. For ex., the ip's used only for the mass matrix are
			 * initialized, too, even if only the stiffness part is assembled. These ip's
			 * are not cleared below as they do not get into the lists, and this creates
			 * issues with the linkers that share subordinated userdata objects. For that,
			 * we clear here all the assigned ip's.
			 *
			 * Should it be done here on in do_fsh_elem_loop?
			 */
			for (size_t i = 0; i < disc->num_imports(); ++i)
			{
				IDataImport<dim>& imp = disc->get_import(i);
				if(imp.data_given())
					imp.data()->clear();
			}
		}
	}
	UG_CATCH_THROW("DataEvaluatorBase::fsh_elem_loop: Cannot finish element loop");

//	clear positions at user data
	/* TODO:
	 * Could it be done in a more elegant way? For ex., why clearing here all the ip
	 * series and not only the ones assigned to the particular userdata objects?
	 */
	clear_positions_in_user_data();
}


///////////////////////////////////////////////////////////////////////////////
// Assemble routines
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void DataEvaluator<TDomain>::prepare_timestep(number future_time, const number time, VectorProxyBase* u, size_t algebra_id)
{
	try
	{
		for (size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_prep_timestep(future_time, time, u, algebra_id);
	}
	UG_CATCH_THROW("DataEvaluatorBase::prep_timestep: Cannot prepare time step.");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
prepare_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	try{
		for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_prep_timestep_elem(time, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::prep_timestep_elem: Cannot prepare timestep.");
}


template <typename TDomain>
void DataEvaluator<TDomain>::
prepare_elem(LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[],
             const LocalIndices& ind,
             bool bDeriv)
{
// 	prepare element
	try{
		for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i){
			UG_DLOG(DID_DATA_EVALUATOR, 2, ">>OCT_DISC_DEBUG: " << "data_evaluator.cpp: " << "DataEvaluatorBase.prepare_elem(): m_vElemDisc[PT_ALL][i]->do_prep_elem() " << roid << std::endl);
			m_vElemDisc[PT_ALL][i]->do_prep_elem(u, elem, roid, vCornerCoords);
		}
	}
	UG_CATCH_THROW("DataEvaluatorBase::prep_elem: Cannot prepare element.");

//	adjust lin defect array of imports and derivative array of exports
//	INFO: This is place here, since the 'prepare_elem' method of an element
//			disc may change the number of integration points, even if the type
//			of the element (e.g. triangle, quad) stays the same. This is the
//			case for, e.g., the NeumannBoundary element disc.
	if(bDeriv)
	{
		for(size_t i = 0; i < m_vImport[PT_ALL][MASS].size(); ++i)
			m_vImport[PT_ALL][MASS][i]->update_dof_sizes(ind);
		for(size_t i = 0; i < m_vImport[PT_ALL][STIFF].size(); ++i)
			m_vImport[PT_ALL][STIFF][i]->update_dof_sizes(ind);
		for(size_t i = 0; i < m_vImport[PT_ALL][RHS].size(); ++i)
			m_vImport[PT_ALL][RHS][i]->update_dof_sizes(ind);

		for(size_t i = 0; i < m_vDependentData.size(); ++i)
			m_vDependentData[i]->update_dof_sizes(ind);
	}

//	evaluate position data
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->compute(&u, elem, vCornerCoords, false);

// 	process dependent data:
//	We can not simply compute exports first, then Linker, because an export
//	itself could depend on other data if implemented somehow in the IElemDisc
//	(e.g. using data from some DataImport). Thus, we have to loop the sorted
//	vector of all dependent data (that is correctly sorted the way that always
//	needed data has previously computed).

//	compute the data
	try{
		if (! time_series_needed ()) { // assemble for the given LocalVector
			for(size_t i = 0; i < m_vDependentData.size(); ++i){
				u.access_by_map(m_vDependentData[i]->map());
				m_vDependentData[i]->compute(&u, elem, vCornerCoords, bDeriv);
			}
		}
		else { // assemble for LocalVectorTimeSeries
			for(size_t i = 0; i < m_vDependentData.size(); ++i){
				u.access_by_map(m_vDependentData[i]->map());
				m_vDependentData[i]->compute(m_pLocTimeSeries, elem, vCornerCoords, bDeriv);
			}
		}
	}
	UG_CATCH_THROW("DataEvaluatorBase::prep_elem: Cannot compute data for Export or Linker.");
}

template <typename TDomain>
void DataEvaluator<TDomain>::finish_timestep(const number time, VectorProxyBase* u, size_t algebra_id)
{
	try
	{
		for (size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_fsh_timestep(time, u, algebra_id);
	}
	UG_CATCH_THROW("DataEvaluatorBase::finish_timestep: Cannot prepare time step.");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
finish_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	try{
		for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
			m_vElemDisc[PT_ALL][i]->do_fsh_timestep_elem(time, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::fsh_timestep_elem: Cannot finish timestep.");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
add_jac_A_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type)
{
	UG_ASSERT(m_discPart & STIFF, "Using add_jac_A_elem, but not STIFF requested.");

	// compute elem-owned contribution
	try{
		for(size_t i = 0; i < m_vElemDisc[type].size(); ++i)
			m_vElemDisc[type][i]->do_add_jac_A_elem(J, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_jac_A_elem: Cannot assemble Jacobian (A)");

	//	compute linearized defect
	try{
		for(size_t i = 0; i < m_vImport[type][STIFF].size(); ++i)
			m_vImport[type][STIFF][i]->compute_lin_defect(u);

		for(size_t i = 0; i < m_vImport[type][RHS].size(); ++i)
			m_vImport[type][RHS][i]->compute_lin_defect(u);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_jac_A_elem: Cannot compute"
			" linearized defect for Import.");

	//	add off diagonal coupling
	try{
		//	loop all imports located in the stiffness part
		for(size_t i = 0; i < m_vImport[type][STIFF].size(); ++i)
			m_vImport[type][STIFF][i]->add_jacobian(J, 1.0);

		//	loop all imports located in the rhs part
		for(size_t i = 0; i < m_vImport[type][RHS].size(); ++i)
			m_vImport[type][RHS][i]->add_jacobian(J, -1.0);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_jac_A_elem: Cannot add couplings.");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
add_jac_M_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type)
{
	UG_ASSERT(m_discPart & MASS, "Using add_jac_M_elem, but not MASS requested.");

	// compute elem-owned contribution
	try{
		for(size_t i = 0; i < m_vElemDisc[type].size(); ++i)
			m_vElemDisc[type][i]->do_add_jac_M_elem(J, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_jac_M_elem: Cannot assemble Jacobian (M)");

	//	compute linearized defect
	try{
		for(size_t i = 0; i < m_vImport[type][MASS].size(); ++i)
			m_vImport[type][MASS][i]->compute_lin_defect(u);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_coupl_JM: Cannot compute"
			" linearized defect for Import.");

	//	add off diagonal coupling
	try{
		//	loop all imports located in the mass part
		for(size_t i = 0; i < m_vImport[type][MASS].size(); ++i)
			m_vImport[type][MASS][i]->add_jacobian(J, 1.0);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_coupl_JM: Cannot add couplings.");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
add_def_A_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type)
{
	UG_ASSERT(m_discPart & STIFF, "Using add_def_A_elem, but not STIFF requested.");

	try{
		for(size_t i = 0; i < m_vElemDisc[type].size(); ++i)
			m_vElemDisc[type][i]->do_add_def_A_elem(d, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_def_A_elem: Cannot assemble Defect (A)");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
add_def_A_expl_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type)
{
	UG_ASSERT(m_discPart & EXPL, "Using add_def_A_elem, but not EXPL requested.");

	try{
		for(size_t i = 0; i < m_vElemDisc[type].size(); ++i)
			m_vElemDisc[type][i]->do_add_def_A_expl_elem(d, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_def_A_expl_elem: Cannot assemble Defect (A)");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
add_def_M_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type)
{
	UG_ASSERT(m_discPart & MASS, "Using add_def_M_elem, but not MASS requested.");

	try{
		for(size_t i = 0; i < m_vElemDisc[type].size(); ++i)
			m_vElemDisc[type][i]->do_add_def_M_elem(d, u, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_def_M_elem: Cannot assemble Defect (M)");
}

template <typename TDomain>
void DataEvaluator<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type)
{
	UG_ASSERT(m_discPart & RHS, "Using add_rhs_elem, but not RHS requested.");

	try{
		for(size_t i = 0; i < m_vElemDisc[type].size(); ++i)
			m_vElemDisc[type][i]->do_add_rhs_elem(rhs, elem, vCornerCoords);
	}
	UG_CATCH_THROW("DataEvaluatorBase::add_rhs_elem: Cannot assemble rhs");
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class DataEvaluatorBase<Domain1d, IElemDisc<Domain1d> >;
template class ErrorEvaluator<Domain1d>;
template class DataEvaluator<Domain1d>;
#endif
#ifdef UG_DIM_2
template class DataEvaluatorBase<Domain2d, IElemDisc<Domain2d> >;
template class ErrorEvaluator<Domain2d>;
template class DataEvaluator<Domain2d>;
#endif
#ifdef UG_DIM_3
template class DataEvaluatorBase<Domain3d, IElemDisc<Domain3d> >;
template class ErrorEvaluator<Domain3d>;
template class DataEvaluator<Domain3d>;
#endif

} // end namespace ug

