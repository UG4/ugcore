/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef CUT_ELEMENT_HANDLER_BASE_IMPL_H_
#define CUT_ELEMENT_HANDLER_BASE_IMPL_H_


namespace ug{

    
template <int TWorldDim>
CutElementHandlerBase<TWorldDim>::
CutElementHandlerBase(SmartPtr<MultiGrid> mg, const char* fctNames,
                             SmartPtr<ParticleProviderSphere<dim> > interfaceProvider)
    : m_spMG(mg.operator->()), m_fctNames(fctNames),
 	  m_spInterfaceProvider(interfaceProvider)
{
    
        m_vThresholdOnLevel.resize(mg->num_levels(), 0.0);
        
        //	get position attachment
        m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
        
        // 	let position accessor access Vertex Coordinates
        if(!m_spMG->has_attachment<Vertex>(m_aPos))
            m_spMG->attach_to<Vertex>(m_aPos);
        m_aaPos.access(*m_spMG, m_aPos);
        
        // initialize data
        m_spCutMarker = make_sp(new BoolMarker);
        m_spOutsideMarker = make_sp(new BoolMarker);
        m_spFlatTopVrtMarker = make_sp(new BoolMarker);
        m_spNearInterfaceVrtMarker = make_sp(new BoolMarker);
        
        
        m_spCutMarker->assign_grid(*m_spMG);
        m_spOutsideMarker->assign_grid(*m_spMG);
        m_spFlatTopVrtMarker->assign_grid(*m_spMG);
        m_spNearInterfaceVrtMarker->assign_grid(*m_spMG);
        
        m_spCutMarker->enable_mark_inheritance(false);
        m_spOutsideMarker->enable_mark_inheritance(false);
        m_spFlatTopVrtMarker->enable_mark_inheritance(false);
        m_spNearInterfaceVrtMarker->enable_mark_inheritance(false);
        
        UG_LOG("init bool_marker done!\n");
        
    }
    

};
