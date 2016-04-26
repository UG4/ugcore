/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

namespace ug{

template<class TAPosition>
inline RefinementProjectionHandler<TAPosition>::
RefinementProjectionHandler(ISubsetHandler* sh, TAPosition aPos)
{
	UG_ASSERT(sh, "Invalid subset handler specified (null pointer)");
	UG_ASSERT(sh->grid(), "Specified subset handler has to operate on a grid!");
	m_sh = sh;
	m_aaPos.access(*m_sh->grid(), aPos);
	m_defaultCallback = SmartPtr<IRefinementCallback>(
					new RefinementCallbackLinear<TAPosition>(*m_sh->grid(), aPos));
}

template<class TAPosition>
inline RefinementProjectionHandler<TAPosition>::
RefinementProjectionHandler(SmartPtr<ISubsetHandler> sh, TAPosition aPos)
{
	UG_ASSERT(sh.valid(), "Invalid subset handler specified (through invalid SmartPtr)");
	UG_ASSERT(sh->grid(), "Specified subset handler has to operate on a grid!");
	m_spSH = sh;
	m_sh = m_spSH.get();
	m_aaPos.access(*m_sh->grid(), aPos);
	m_defaultCallback = SmartPtr<IRefinementCallback>(
					new RefinementCallbackLinear<TAPosition>(*m_sh->grid(), aPos));
}

template<class TAPosition>
inline RefinementProjectionHandler<TAPosition>::
~RefinementProjectionHandler()
{
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
set_default_callback(SmartPtr<IRefinementCallback> callback)
{
	m_defaultCallback = callback;
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
set_callback(int subsetIndex, SmartPtr<IRefinementCallback> callback)
{
	UG_COND_THROW(subsetIndex < 0,
			"Invalid subset specified during refinement projector registration");

	if(subsetIndex >= (int)m_callbacks.size())
		m_callbacks.resize(subsetIndex + 1);

	m_callbacks[subsetIndex] = callback;
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
set_callback(std::string subsetName, SmartPtr<IRefinementCallback> callback)
{
	set_callback(m_sh->get_subset_index(subsetName.c_str()), callback);
}

template<class TAPosition>
template <class TParent>
inline void RefinementProjectionHandler<TAPosition>::
handle_new_vertex(Vertex* vrt, TParent* parent)
{
	UG_ASSERT(parent, "A vertex can only be created from a parent element!");

	int si = m_sh->get_subset_index(parent);
	if((si == -1) || (si >= (int)m_callbacks.size()) || (!m_callbacks[si].valid()))
		m_defaultCallback->new_vertex(vrt, parent);
	else
		m_callbacks[si]->new_vertex(vrt, parent);
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	handle_new_vertex(vrt, parent);
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	handle_new_vertex(vrt, parent);
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
	handle_new_vertex(vrt, parent);
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
	handle_new_vertex(vrt, parent);
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
flat_grid_vertex_encountered(Vertex* vrt)
{
	int si = m_sh->get_subset_index(vrt);
	if((si == -1) || (si >= (int)m_callbacks.size()) || (!m_callbacks[si].valid()))
		m_defaultCallback->flat_grid_vertex_encountered(vrt);
	else
		m_callbacks[si]->flat_grid_vertex_encountered(vrt);
}

template<class TAPosition>
inline int RefinementProjectionHandler<TAPosition>::
current_pos(number* coordsOut, Vertex* vrt, int maxCoords)
{
	return IRefinementCallback::current_pos_helper(coordsOut, vrt, maxCoords, m_aaPos);
}

}//	end of namespace
