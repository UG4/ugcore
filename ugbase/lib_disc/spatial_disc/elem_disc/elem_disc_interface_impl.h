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

template <typename TElem>
inline void IElemDisc::fast_prepare_timestep_elem(TElem* elem, const LocalVector& u)
{
//	cast the method pointer back to the original type
	typedef void (IElemDisc::*Func)(TElem*, const LocalVector&);
	Func pFunc = reinterpret_cast<Func>(m_vPrepareTimestepElemFct[m_id]);
	(this->*(pFunc))(elem, u);
}

template <typename TElem>
inline void IElemDisc::fast_prepare_elem(TElem* elem, const LocalVector& u)
{
//	cast the method pointer back to the original type
	typedef void (IElemDisc::*Func)(TElem*, const LocalVector&);
	Func pFunc = reinterpret_cast<Func>(m_vPrepareElemFct[m_id]);
	(this->*(pFunc))(elem, u);
}

template <typename TElem>
inline void IElemDisc::fast_finish_timestep_elem(TElem* elem, const number time, const LocalVector& u)
{
//	cast the method pointer back to the original type
	typedef void (IElemDisc::*Func)(TElem*, const LocalVector&);
	Func pFunc = reinterpret_cast<Func>(m_vFinishTimestepElemFct[m_id]);
	(this->*(pFunc))(elem, u);
}

///////////////////////////////////////////////////////////////////////////////
//	setting of fast elem ass functions
///////////////////////////////////////////////////////////////////////////////

template<typename TAssFunc>
void IElemDisc::set_prep_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vPrepareTimestepElemFct[id] = reinterpret_cast<PrepareTimestepElemFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_prep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vPrepareElemFct[id] = reinterpret_cast<PrepareElemFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_prep_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vPrepareElemLoopFct[id] = static_cast<PrepareElemLoopFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_fsh_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vFinishElemLoopFct[id] = static_cast<FinishElemLoopFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_ass_JA_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemJAFct[id] = static_cast<ElemJAFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_ass_JM_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemJMFct[id] = static_cast<ElemJMFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_ass_dA_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemdAFct[id] = static_cast<ElemdAFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_ass_dM_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemdMFct[id] = static_cast<ElemdMFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_ass_rhs_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vElemRHSFct[id] = static_cast<ElemRHSFct>(func);
};

template<typename TAssFunc>
void IElemDisc::set_fsh_timestep_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	if(!fast_ass_elem_enabled())
		UG_THROW("IElemDisc: trying to register fast_ass_elem-function, but"
						" IElemDisc has not been set to fast assembling. Please"
						" use 'enable_fast_ass_elem(true)' in your IElemDisc "
						" prior to the setting of any fast_ass_elem-function.");
	m_vFinishTimestepElemFct[id] = reinterpret_cast<FinishTimestepElemFct>(func);
};

}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
