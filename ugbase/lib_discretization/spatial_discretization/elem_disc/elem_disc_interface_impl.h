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
void IElemDisc::reg_prepare_vol_loop_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElemLoopFunc.size())
		m_vPrepareElemLoopFunc.resize(id+1, 0);

	m_vPrepareElemLoopFunc[id] = (PrepareElementLoopFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_prepare_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElemFunc.size())
		m_vPrepareElemFunc.resize(id+1, 0);

	m_vPrepareElemFunc[id] = (PrepareElementFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_finish_vol_loop_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vFinishElemLoopFunc.size())
		m_vFinishElemLoopFunc.resize(id+1, 0);

	m_vFinishElemLoopFunc[id] = (FinishElementLoopFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_JA_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssJAFunc.size())
		m_vAssJAFunc.resize(id+1, 0);

	m_vAssJAFunc[id] = (AssembleJAFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_JM_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssJMFunc.size())
		m_vAssJMFunc.resize(id+1, 0);

	m_vAssJMFunc[id] = (AssembleJMFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_dA_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssAFunc.size())
		m_vAssAFunc.resize(id+1, 0);

	m_vAssAFunc[id] = (AssembleAFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_dM_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssMFunc.size())
		m_vAssMFunc.resize(id+1, 0);

	m_vAssMFunc[id] = (AssembleMFunc)func;
};

template<typename TAssFunc>
void IElemDisc::reg_ass_rhs_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssFFunc.size())
		m_vAssFFunc.resize(id+1, 0);

	m_vAssFFunc[id] = (AssembleFFunc)func;
};

}

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
