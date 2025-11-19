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


#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__

#include "elem_disc_interface.h"
#include "common/util/string_util.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
//	setting of elem ass functions
///////////////////////////////////////////////////////////////////////////////

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_prep_timestep_fct(size_t algebra_id, TAssFunc func)
{
	m_vPrepareTimestepFct[algebra_id] = static_cast<PrepareTimestepFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_prep_timestep_fct(size_t algebra_id)
{
	m_vPrepareTimestepFct[algebra_id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_prep_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareTimestepElemFct[id] = static_cast<PrepareTimestepElemFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_prep_timestep_elem_fct(ReferenceObjectID id)
{
	m_vPrepareTimestepElemFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_prep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareElemFct[id] = static_cast<PrepareElemFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_prep_elem_fct(ReferenceObjectID id)
{
	m_vPrepareElemFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_prep_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareElemLoopFct[id] = static_cast<PrepareElemLoopFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_prep_elem_loop_fct(ReferenceObjectID id)
{
	m_vPrepareElemLoopFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_fsh_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishElemLoopFct[id] = static_cast<FinishElemLoopFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_fsh_elem_loop_fct(ReferenceObjectID id)
{
	m_vFinishElemLoopFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_add_jac_A_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemJAFct[id] = static_cast<ElemJAFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_add_jac_A_elem_fct(ReferenceObjectID id)
{
	m_vElemJAFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_add_jac_M_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemJMFct[id] = static_cast<ElemJMFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_add_jac_M_elem_fct(ReferenceObjectID id)
{
	m_vElemJMFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_add_def_A_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdAFct[id] = static_cast<ElemdAFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_add_def_A_elem_fct(ReferenceObjectID id)
{
	m_vElemdAFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_add_def_A_expl_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdAExplFct[id] = static_cast<ElemdAFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_add_def_A_expl_elem_fct(ReferenceObjectID id)
{
	m_vElemdAExplFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_add_def_M_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdMFct[id] = static_cast<ElemdMFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_add_def_M_elem_fct(ReferenceObjectID id)
{
	m_vElemdMFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_add_rhs_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemRHSFct[id] = static_cast<ElemRHSFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_add_rhs_elem_fct(ReferenceObjectID id)
{
	m_vElemRHSFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_fsh_timestep_fct(size_t algebra_id, TAssFunc func)
{
	m_vFinishTimestepFct[algebra_id] = static_cast<FinishTimestepFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_fsh_timestep_fct(size_t algebra_id)
{
	m_vFinishTimestepFct[algebra_id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemAssembleFuncs<TLeaf, TDomain>::set_fsh_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishTimestepElemFct[id] = static_cast<FinishTimestepElemFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::remove_fsh_timestep_elem_fct(ReferenceObjectID id)
{
	m_vFinishTimestepElemFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_prep_err_est_elem_loop(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareErrEstElemLoopFct[id] = static_cast<PrepareErrEstElemLoopFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::remove_prep_err_est_elem_loop(ReferenceObjectID id)
{
	m_vPrepareErrEstElemLoopFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_prep_err_est_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareErrEstElemFct[id] = static_cast<PrepareErrEstElemFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::remove_prep_err_est_elem(ReferenceObjectID id)
{
	m_vPrepareErrEstElemFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_compute_err_est_A_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vElemComputeErrEstAFct[id] = static_cast<ElemComputeErrEstAFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::remove_compute_err_est_A_elem(ReferenceObjectID id)
{
	m_vElemComputeErrEstAFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_compute_err_est_M_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vElemComputeErrEstMFct[id] = static_cast<ElemComputeErrEstMFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::remove_compute_err_est_M_elem(ReferenceObjectID id)
{
	m_vElemComputeErrEstMFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_compute_err_est_rhs_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vElemComputeErrEstRhsFct[id] = static_cast<ElemComputeErrEstRhsFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::remove_compute_err_est_rhs_elem(ReferenceObjectID id)
{
	m_vElemComputeErrEstRhsFct[id] = nullptr;
};

template <typename TLeaf, typename TDomain>
template<typename TAssFunc>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_fsh_err_est_elem_loop(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishErrEstElemLoopFct[id] = static_cast<FinishErrEstElemLoopFct>(func);
};
template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::remove_fsh_err_est_elem_loop(ReferenceObjectID id)
{
	m_vFinishErrEstElemLoopFct[id] = nullptr;
};

}

#endif