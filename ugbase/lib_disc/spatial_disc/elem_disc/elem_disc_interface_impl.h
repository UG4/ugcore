/*
 * elem_disc_interface_impl.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */


#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__

#include "elem_disc_interface.h"
#include "common/util/string_util.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
//	setting of elem ass functions
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_timestep_fct(size_t algebra_id, TAssFunc func)
{
	m_vPrepareTimestepFct[algebra_id] = static_cast<PrepareTimestepFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_prep_timestep_fct(size_t algebra_id)
{
	m_vPrepareTimestepFct[algebra_id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareTimestepElemFct[id] = static_cast<PrepareTimestepElemFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_prep_timestep_elem_fct(ReferenceObjectID id)
{
	m_vPrepareTimestepElemFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareElemFct[id] = static_cast<PrepareElemFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_prep_elem_fct(ReferenceObjectID id)
{
	m_vPrepareElemFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareElemLoopFct[id] = static_cast<PrepareElemLoopFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_prep_elem_loop_fct(ReferenceObjectID id)
{
	m_vPrepareElemLoopFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_fsh_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishElemLoopFct[id] = static_cast<FinishElemLoopFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_fsh_elem_loop_fct(ReferenceObjectID id)
{
	m_vFinishElemLoopFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_jac_A_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemJAFct[id] = static_cast<ElemJAFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_add_jac_A_elem_fct(ReferenceObjectID id)
{
	m_vElemJAFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_jac_M_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemJMFct[id] = static_cast<ElemJMFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_add_jac_M_elem_fct(ReferenceObjectID id)
{
	m_vElemJMFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_def_A_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdAFct[id] = static_cast<ElemdAFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_add_def_A_elem_fct(ReferenceObjectID id)
{
	m_vElemdAFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_def_A_expl_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdAExplFct[id] = static_cast<ElemdAFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_add_def_A_expl_elem_fct(ReferenceObjectID id)
{
	m_vElemdAExplFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_def_M_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdMFct[id] = static_cast<ElemdMFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_add_def_M_elem_fct(ReferenceObjectID id)
{
	m_vElemdMFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_rhs_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemRHSFct[id] = static_cast<ElemRHSFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_add_rhs_elem_fct(ReferenceObjectID id)
{
	m_vElemRHSFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_fsh_timestep_fct(size_t algebra_id, TAssFunc func)
{
	m_vFinishTimestepFct[algebra_id] = static_cast<FinishTimestepFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_fsh_timestep_fct(size_t algebra_id)
{
	m_vFinishTimestepFct[algebra_id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_fsh_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishTimestepElemFct[id] = static_cast<FinishTimestepElemFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_fsh_timestep_elem_fct(ReferenceObjectID id)
{
	m_vFinishTimestepElemFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_err_est_elem_loop(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareErrEstElemLoopFct[id] = static_cast<PrepareErrEstElemLoopFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_prep_err_est_elem_loop(ReferenceObjectID id)
{
	m_vPrepareErrEstElemLoopFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_err_est_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareErrEstElemFct[id] = static_cast<PrepareErrEstElemFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_prep_err_est_elem(ReferenceObjectID id)
{
	m_vPrepareErrEstElemFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_compute_err_est_A_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vElemComputeErrEstAFct[id] = static_cast<ElemComputeErrEstAFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_compute_err_est_A_elem(ReferenceObjectID id)
{
	m_vElemComputeErrEstAFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_compute_err_est_M_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vElemComputeErrEstMFct[id] = static_cast<ElemComputeErrEstMFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_compute_err_est_M_elem(ReferenceObjectID id)
{
	m_vElemComputeErrEstMFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_compute_err_est_rhs_elem(ReferenceObjectID id, TAssFunc func)
{
	m_vElemComputeErrEstRhsFct[id] = static_cast<ElemComputeErrEstRhsFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_compute_err_est_rhs_elem(ReferenceObjectID id)
{
	m_vElemComputeErrEstRhsFct[id] = NULL;
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_fsh_err_est_elem_loop(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishErrEstElemLoopFct[id] = static_cast<FinishErrEstElemLoopFct>(func);
};
template <typename TDomain>
void IElemDisc<TDomain>::remove_fsh_err_est_elem_loop(ReferenceObjectID id)
{
	m_vFinishErrEstElemLoopFct[id] = NULL;
};

}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
