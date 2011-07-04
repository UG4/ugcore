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

template<typename TAssFunc>
void IElemDisc::reg_prepare_elem_loop_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElemLoopFct.size())
		m_vPrepareElemLoopFct.resize(id+1, 0);

	m_vPrepareElemLoopFct[id] = (PrepareElemLoopFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_prepare_elem_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElemFct.size())
		m_vPrepareElemFct.resize(id+1, 0);

	m_vPrepareElemFct[id] = (PrepareElemFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_finish_elem_loop_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vFinishElemLoopFct.size())
		m_vFinishElemLoopFct.resize(id+1, 0);

	m_vFinishElemLoopFct[id] = (FinishElemLoopFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_JA_elem_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vElemJAFct.size())
		m_vElemJAFct.resize(id+1, 0);

	m_vElemJAFct[id] = (ElemJAFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_JM_elem_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vElemJMFct.size())
		m_vElemJMFct.resize(id+1, 0);

	m_vElemJMFct[id] = (ElemJMFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_dA_elem_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vElemdAFct.size())
		m_vElemdAFct.resize(id+1, 0);

	m_vElemdAFct[id] = (ElemdAFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_dM_elem_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vElemdMFct.size())
		m_vElemdMFct.resize(id+1, 0);

	m_vElemdMFct[id] = (ElemdMFct)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_rhs_elem_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vElemRHSFct.size())
		m_vElemRHSFct.resize(id+1, 0);

	m_vElemRHSFct[id] = (ElemRHSFct)func;
};

}

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
