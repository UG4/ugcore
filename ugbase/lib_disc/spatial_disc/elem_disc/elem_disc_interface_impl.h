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

template <typename TDomain>
template <typename TElem>
inline void IElemDisc<TDomain>::fast_prep_timestep_elem(TElem* elem, const LocalVector& u)
{
	if (this->m_vPrepareTimestepElemFct[m_id] != NULL)
	{
	//	cast the method pointer back to the original type
		typedef void (IElemDisc<TDomain>::*Func)(TElem*, const LocalVector&);
		Func pFunc = reinterpret_cast<Func>(m_vPrepareTimestepElemFct[m_id]);
		(this->*(pFunc))(elem, u);
	}
}

template <typename TDomain>
template <typename TElem>
inline void IElemDisc<TDomain>::fast_prep_elem(TElem* elem, const LocalVector& u)
{
	UG_ASSERT(m_vPrepareElemFct[m_id]!=NULL, "Fast-Assemble Method missing.");
//	cast the method pointer back to the original type
	typedef void (IElemDisc<TDomain>::*Func)(TElem*, const LocalVector&);
	Func pFunc = reinterpret_cast<Func>(m_vPrepareElemFct[m_id]);
	(this->*(pFunc))(elem, u);
}

template <typename TDomain>
template <typename TElem>
inline void IElemDisc<TDomain>::fast_fsh_timestep_elem(TElem* elem, const number time, const LocalVector& u)
{
	if (this->m_vFinishTimestepElemFct[m_id] != NULL)
	{
	//	cast the method pointer back to the original type
		typedef void (IElemDisc<TDomain>::*Func)(TElem*, const LocalVector&);
		Func pFunc = reinterpret_cast<Func>(m_vFinishTimestepElemFct[m_id]);
		(this->*(pFunc))(elem, u);
	}
}

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
	m_vPrepareTimestepElemFct[id] = reinterpret_cast<PrepareTimestepElemFct>(func);
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
	m_vPrepareElemFct[id] = reinterpret_cast<PrepareElemFct>(func);
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

/////////////////////////////////////////////////////////////////////////////

//explicit reaction, reaction_rate and source

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

/////////////////////////////////////////////////////////////////////////////

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
	m_vFinishTimestepElemFct[id] = reinterpret_cast<FinishTimestepElemFct>(func);
};

}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
