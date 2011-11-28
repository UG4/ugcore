// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.11.2011 (m,d,y)

#ifndef __H__UG__surface_view_impl_TMP__
#define __H__UG__surface_view_impl_TMP__

namespace ug{
namespace tmp{

template <class TElem>
bool SurfaceViewElementIterator<TElem>::
equal(SurfaceViewElementIterator<TElem> const& other) const
{
	return m_elemIter == other.m_elemIter;
}

template <class TElem>
void SurfaceViewElementIterator<TElem>::
increment()
{
	++m_elemIter;
	while(m_lvl < m_topLvl){
		MGSubsetHandler& sh = m_surfView->subset_handler();
		if(m_elemIter == sh.end<TElem>(m_si, m_lvl)){
		//	we reached the end of the current level.
		//	ascend one level
			++m_lvl;
			m_elemIter = sh.begin<TElem>(m_si, m_lvl);
		}
		else{
			while(m_surfView->is_shadow(*m_elemIter)){
				++m_elemIter;
				if(m_elemIter == sh.end<TElem>(m_si, m_lvl)){
				//	we reached the end of the current level.
				//	ascend one level and start again
					++m_lvl;
					m_elemIter = sh.begin<TElem>(m_si, m_lvl);
					break;
				}
			}
		}
	}

	return *this;
}

template <class TElem>
TElem& SurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}

}//	end of namespace
}//	end of namespace

#endif
