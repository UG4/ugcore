/*
 * elem_disc_interface_impl.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */


#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__

#include "elem_disc_interface.h"
#include "common/util/string_util.h"

namespace ug{

template <typename TElem>
inline bool IElemDisc::prepare_elem(TElem* elem, const local_vector_type& u)
{
//	cast the method pointer back to the original type
	typedef bool (IElemDisc::*Func)(TElem*, const local_vector_type&);
	Func pFunc = reinterpret_cast<Func>(m_vPrepareElemFct[m_id]);
	return (this->*(pFunc))(elem, u);
}


template<typename TAssFunc>
void IElemDisc::reg_prepare_elem_fct(ReferenceObjectID id, TAssFunc func)
{
//	we cast the method pointer to a different type
	m_vPrepareElemFct[id] = reinterpret_cast<PrepareElemFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_prepare_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vPrepareElemLoopFct[id] = static_cast<PrepareElemLoopFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_finish_elem_loop_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vFinishElemLoopFct[id] = static_cast<FinishElemLoopFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_ass_JA_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemJAFct[id] = static_cast<ElemJAFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_ass_JM_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemJMFct[id] = static_cast<ElemJMFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_ass_dA_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdAFct[id] = static_cast<ElemdAFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_ass_dM_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemdMFct[id] = static_cast<ElemdMFct>(func);
};

template<typename TAssFunc>
void IElemDisc::reg_ass_rhs_elem_fct(ReferenceObjectID id, TAssFunc func)
{
	m_vElemRHSFct[id] = static_cast<ElemRHSFct>(func);
};

}

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
