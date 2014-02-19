//	created by Sebastian Reiter
//	s.b.reiter@gmail.com
//	april 2013

namespace ug{

template<class TAPosition>
inline RefinementProjectionHandler<TAPosition>::
RefinementProjectionHandler(SmartPtr<ISubsetHandler> sh, TAPosition aPos)
{
	m_sh = sh;
	UG_ASSERT(m_sh.valid(), "Invalid smart pointer specified");
	UG_ASSERT(m_sh->grid(), "Specified subset handler has to operate on a grid!");
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
	if(subsetIndex < 0)
		return;

	if(subsetIndex >= (int)m_callbacks.size())
		m_callbacks.resize(subsetIndex + 1);

	m_callbacks[subsetIndex] = callback;
}

template<class TAPosition>
inline void RefinementProjectionHandler<TAPosition>::
set_callback(std::string subsetName, SmartPtr<IRefinementCallback> callback)
{
	UG_ASSERT(m_sh.valid(), "Invalid smart pointer specified");
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
new_vertex(Vertex* vrt, EdgeBase* parent)
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
