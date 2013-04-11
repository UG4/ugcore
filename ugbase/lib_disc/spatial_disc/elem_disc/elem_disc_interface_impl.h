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
//	setting of fast elem ass functions
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vPrepareTimestepElemFct[id] = static_cast<PrepareTimestepElemFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vPrepareElemFct[id] = static_cast<PrepareElemFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_prep_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vPrepareElemLoopFct[id] = static_cast<PrepareElemLoopFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_fsh_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vFinishElemLoopFct[id] = static_cast<FinishElemLoopFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_jac_A_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemJAFct[id] = static_cast<ElemJAFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_jac_M_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemJMFct[id] = static_cast<ElemJMFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_def_A_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemdAFct[id] = static_cast<ElemdAFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_def_A_expl_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem_explicit-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem_explicit(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem_explicit-function.");
	m_vElemdAExplFct[id] = static_cast<ElemdAFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_def_M_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemdMFct[id] = static_cast<ElemdMFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_add_rhs_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemRHSFct[id] = static_cast<ElemRHSFct>(func);
};

template <typename TDomain>
template<typename TAssFunc>
void IElemDisc<TDomain>::set_fsh_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_add_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_add_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vFinishTimestepElemFct[id] = static_cast<FinishTimestepElemFct>(func);
};

}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
