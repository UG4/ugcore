// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)
 
#include "surface_view.h"
#include "common/assert.h"
#include "lib_grid/parallelization/util/compol_boolmarker.h"
#ifdef UG_PARALLEL
	#include "pcl/pcl_interface_communicator.h"
	#include "lib_grid/parallelization/copy_policy.h"
#endif

namespace ug{

///	adds marking at extracting side
template <class TLayout>
class ComPol_GatherSurfaceStates : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;

	///	Construct the communication policy with a ug::BoolMarker.
		ComPol_GatherSurfaceStates(MultiGrid& mg,
						MultiElementAttachmentAccessor<AByte>& aaElemSurfState)
			 :	m_mg(mg), m_aaESS(aaElemSurfState)
		{}

		virtual ~ComPol_GatherSurfaceStates()	{}

		virtual int get_required_buffer_size(Interface& interface)
		{
			return interface.size() * sizeof(byte);
		}

	///	write surface state for each entry
		virtual bool collect(ug::BinaryBuffer& buff, Interface& intfc)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = intfc.begin(); iter != intfc.end(); ++iter)
			{
				Element elem = intfc.get_element(iter);
				buff.write((char*)&m_aaESS[elem], sizeof(byte));
			}
			return true;
		}

	///	reads marks from the given stream
		virtual bool extract(ug::BinaryBuffer& buff, Interface& intfc)
		{
			for(InterfaceIter iter = intfc.begin(); iter != intfc.end(); ++iter)
			{
				Element elem = intfc.get_element(iter);
				byte nv;
				buff.read((char*)&nv, sizeof(byte));
				m_aaESS[elem] |= nv;
			}
			return true;
		}

	protected:
		MultiGrid&								m_mg;
		MultiElementAttachmentAccessor<AByte>	m_aaESS;
};

////////////////////////////////////////////////////////////////////////////////
//	Create Surface View
////////////////////////////////////////////////////////////////////////////////
template <class TElem>
bool SurfaceView::
is_local_surface_view_element(TElem* elem)
{
	#ifdef UG_PARALLEL
		return !(m_pMG->has_children(elem)
				|| m_distGridMgr->contains_status(elem, ES_V_MASTER));
	#else
		return !m_pMG->has_children(elem);
	#endif
}

void SurfaceView::
refresh_surface_states()
{
	if(m_pMG->num_volumes() > 0)
		refresh_surface_states<Volume>();
	else if(m_pMG->num_faces() > 0)
		refresh_surface_states<Face>();
	else if(m_pMG->num_edges() > 0)
		refresh_surface_states<EdgeBase>();
	else
		refresh_surface_states<VertexBase>();
}

template <class TElem>
void SurfaceView::
refresh_surface_states()
{
//	some typedefs
	typedef typename geometry_traits<TElem>::iterator ElemIter;

	MultiGrid& mg = *m_pMG;

//	reset surface states of all elements. Initially, we'll set all states to hidden
	SetAttachmentValues(m_aaElemSurfState, mg.begin<VertexBase>(), mg.end<VertexBase>(), ESS_NONE);
	SetAttachmentValues(m_aaElemSurfState, mg.begin<EdgeBase>(), mg.end<EdgeBase>(), ESS_NONE);
	SetAttachmentValues(m_aaElemSurfState, mg.begin<Face>(), mg.end<Face>(), ESS_NONE);
	SetAttachmentValues(m_aaElemSurfState, mg.begin<Volume>(), mg.end<Volume>(), ESS_NONE);

//	iterate through all levels of the mgsh
	for(size_t level = 0; level < mg.num_levels(); ++level){

	//	iterate through all elements on that level
		for(ElemIter iter = mg.begin<TElem>(level);
			iter != mg.end<TElem>(level); ++iter)
		{
			TElem* elem = *iter;

		//	assign local surface state. In a serial environment, this is sufficient.
		//	In a parallel environment, we possibly mark vertical masters as surface
		//	elements, even though they have children on other processes. This will
		//	be fixed in a post-processing step.
			if(is_local_surface_view_element(elem)){
				set_surface_state(elem, ESS_SURFACE);
				mark_sides_as_surface_or_shadow<TElem, typename TElem::side>(elem);
				if(GeometricObject* p = m_pMG->get_parent(elem))
					set_surface_state(p, surface_state(p) | ESS_HIDDEN);
			}
		}
	}
/*
//	in a parallel environment, we'll mark all vertical masters as unknown
	mark_vmasters_as_unknown<VertexBase>();
	mark_vmasters_as_unknown<EdgeBase>();
	mark_vmasters_as_unknown<Face>();
	mark_vmasters_as_unknown<Volume>();
*/
//	we have to make sure that all copies have the same surface states on all processes
	adjust_parallel_surface_states<VertexBase>();
	adjust_parallel_surface_states<EdgeBase>();
	adjust_parallel_surface_states<Face>();
	adjust_parallel_surface_states<Volume>();
}

template <class TElem, class TSide>
void SurfaceView::
mark_sides_as_surface_or_shadow(TElem* elem)
{
	if(!TElem::HAS_SIDES)
		return;

	typename Grid::traits<TSide>::secure_container	sides;

	m_pMG->associated_elements(sides, elem);
	for(size_t i = 0; i < sides.size(); ++i){
		TSide* s = sides[i];
		if(surface_state(s) == ESS_NONE){
			set_surface_state(s, ESS_SURFACE);
			if(GeometricObject* p = m_pMG->get_parent(s))
				set_surface_state(p, surface_state(p) | ESS_HIDDEN);
		}
	}

	if(TSide::HAS_SIDES)
		mark_sides_as_surface_or_shadow<TElem, typename TSide::side>(elem);
}

template <class TElem>
void SurfaceView::
mark_vmasters_as_unknown()
{
/*
	#ifdef UG_PARALLEL
		typedef typename GridLayoutMap::Types<TElem>::Layout	Layout;
		typedef typename Layout::iterator					IntfcIter;
		typedef typename Layout::Interface					Intfc;
		typedef typename Intfc::iterator					ElemIter;

		if(m_distGridMgr->grid_layout_map().has_layout<TElem>(INT_V_MASTER)){
			Layout& layout = m_distGridMgr->grid_layout_map().get_layout<TElem>(INT_V_MASTER);

			for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
				for(IntfcIter iiter = layout.begin(lvl); iiter != layout.end(lvl); ++iiter)
				{
					Intfc& intfc = layout.interface(iiter);
					for(ElemIter eiter = intfc.begin(); eiter != intfc.end(); ++eiter)
					{
						set_surface_state(intfc.get_element(eiter), ESS_UNKNOWN);
					}
				}
			}
		}
	#endif
*/
}

template <class TElem>
void SurfaceView::
adjust_parallel_surface_states()
{
	#ifdef UG_PARALLEL
		typedef typename GridLayoutMap::Types<TElem>::Layout	Layout;

		GridLayoutMap& glm = m_distGridMgr->grid_layout_map();
		ComPol_GatherSurfaceStates<Layout>	cpAdjust(*m_pMG, m_aaElemSurfState);
		pcl::InterfaceCommunicator<Layout> com;

		com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, cpAdjust);

		com.communicate();

		CopyPolicy<Layout, AByte> cpCopyStates(*m_pMG, m_aElemSurfState);
		com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, cpCopyStates);
		com.communicate();
	#endif
}


////////////////////////////////////////////////////////////////////////////////
// SurfaceView
////////////////////////////////////////////////////////////////////////////////

SurfaceView::SurfaceView(SmartPtr<MGSubsetHandler> spMGSH,
                         bool adaptiveMG) :
	m_spMGSH(spMGSH),
	m_adaptiveMG(adaptiveMG),
	m_pMG(m_spMGSH->multi_grid()),
	m_distGridMgr(m_spMGSH->multi_grid()->distributed_grid_manager())
{
	UG_ASSERT(m_pMG, "A MultiGrid has to be assigned to the given subset handler");

	m_pMG->attach_to_all_dv(m_aElemSurfState, 0);
	m_aaElemSurfState.access(*m_pMG, m_aElemSurfState);

	refresh_surface_states();
}

SurfaceView::~SurfaceView()
{
	m_pMG->detach_from_all(m_aElemSurfState);
}



SurfaceLevelView::
SurfaceLevelView(SmartPtr<SurfaceView> spSV, int topLvl) :
	m_spSV(spSV),
	m_topLvl(topLvl)
{
}

}// end of namespace
