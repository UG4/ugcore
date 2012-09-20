// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)
 
#include "surface_view.h"
#include "common/assert.h"
#include "lib_grid/parallelization/util/compol_boolmarker.h"

namespace ug{


////////////////////////////////////////////////////////////////////////
template <class TElem, class TSideElem>
void MarkAssociated(BoolMarker& boolMarker)
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	
	typename Grid::traits<TSideElem>::secure_container sides;
	Grid& grid = *boolMarker.grid();

	iterator iterEnd = grid.end<TElem>();
	for(iterator iter = grid.begin<TElem>(); iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;
		if(!boolMarker.is_marked(elem)) continue;

		grid.associated_elements(sides, elem);
		
		for(size_t i = 0; i < sides.size(); ++i)
			boolMarker.mark(sides[i]);
	}
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem>
void MarkAssociatedLowerDimElems(BoolMarker& boolMarker,
                                 const Volume&)
{
//	we have to find all associated elements of lower dimension.
	MarkAssociated<TElem, Face>(boolMarker);
	MarkAssociated<TElem, EdgeBase>(boolMarker);
	MarkAssociated<TElem, VertexBase>(boolMarker);
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem>
void MarkAssociatedLowerDimElems(BoolMarker& boolMarker,
                                 const Face&)
{
//	we have to find all associated elements of lower dimension.
	MarkAssociated<TElem, EdgeBase>(boolMarker);
	MarkAssociated<TElem, VertexBase>(boolMarker);
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem>
void MarkAssociatedLowerDimElems(BoolMarker& boolMarker,
                                 const EdgeBase&)
{
//	we have to find all associated elements of lower dimension.
	MarkAssociated<TElem, VertexBase>(boolMarker);
}

template <class TElem>
void MarkAssociatedLowerDimElems(BoolMarker& boolMarker)
{
	MarkAssociatedLowerDimElems<TElem>(boolMarker, TElem());
}

template <class TElem>
void MarkAssociatedSides(BoolMarker& boolMarker)
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	typedef typename TElem::lower_dim_base_object Side;

	typename Grid::traits<Side>::secure_container sides;
	Grid& grid = *boolMarker.grid();

	iterator iterEnd = grid.end<TElem>();
	for(iterator iter = grid.begin<TElem>(); iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;
		if(!boolMarker.is_marked(elem)) continue;

		grid.associated_elements(sides, elem);

		for(size_t i = 0; i < sides.size(); ++i)
			boolMarker.mark(sides[i]);
	}
}

////////////////////////////////////////////////////////////////////////////////
//	Create Surface View
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
void SetSurfaceViewMarks(BoolMarker& boolMarker,
                         MultiGrid* pMG
#ifdef UG_PARALLEL
                         ,DistributedGridManager& distGridMgr
#endif
                         )
{
//	some typedefs
	typedef typename geometry_traits<TElem>::iterator ElemIter;

//	iterate through all levels of the mgsh
	for(size_t level = 0; level < pMG->num_levels(); ++level){

	//	iterate through all elements in the subset on that level
		for(ElemIter iter = pMG->begin<TElem>(level);
			iter != pMG->end<TElem>(level); ++iter)
		{
			TElem* elem = *iter;

		//	check whether the element has children
#ifdef UG_PARALLEL
			if(distGridMgr.is_ghost(elem)) continue;
#endif
			if(pMG->has_children(elem)) continue;

			boolMarker.mark(elem);
		}
	}
}

template <class TElem>
void RemoveSurfaceViewMarks(BoolMarker& boolMarker,
                            MultiGrid* pMG
#ifdef UG_PARALLEL
                            ,DistributedGridManager& distGridMgr
#endif
                            )
{
//	some typedefs
	typedef typename geometry_traits<TElem>::iterator ElemIter;

//	iterate through all levels of the mgsh
	for(size_t level = 0; level < pMG->num_levels(); ++level){

	//	iterate through all elements in the subset on that level
		for(ElemIter iter = pMG->begin<TElem>(level);
			iter != pMG->end<TElem>(level); ++iter)
		{
			TElem* elem = *iter;

		//	check whether the element has children
#ifdef UG_PARALLEL
			if(distGridMgr.is_ghost(elem)) continue;
#endif
			if(pMG->has_children(elem)) continue;

			boolMarker.unmark(elem);
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	Marks shadow elements. m_bMarker only holds elements which are shadows afterwards.
void MarkShadows(BoolMarker& boolMarker,
                 MultiGrid* pMG
#ifdef UG_PARALLEL
                 ,DistributedGridManager& distGridMgr
#endif
				)
{
//	Marks all elements, that do not have children and are non-ghosts
	boolMarker.clear();
#ifdef UG_PARALLEL
	SetSurfaceViewMarks<VertexBase>(boolMarker, pMG, distGridMgr);
	SetSurfaceViewMarks<EdgeBase>(boolMarker, pMG, distGridMgr);
	SetSurfaceViewMarks<Face>(boolMarker, pMG, distGridMgr);
	SetSurfaceViewMarks<Volume>(boolMarker, pMG, distGridMgr);
#else
	SetSurfaceViewMarks<VertexBase>(boolMarker, pMG);
	SetSurfaceViewMarks<EdgeBase>(boolMarker, pMG);
	SetSurfaceViewMarks<Face>(boolMarker, pMG);
	SetSurfaceViewMarks<Volume>(boolMarker, pMG);
#endif

//	assign associated elements of lower dimension to the surface view
	bool assignSidesOnly = true;
	if((pMG->num<Volume>() > 0) && !pMG->option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		assignSidesOnly = false;
	else if((pMG->num<Volume>() > 0) && !(pMG->option_is_enabled(VOLOPT_AUTOGENERATE_EDGES
										  || pMG->option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))))
		assignSidesOnly = false;
	else if((pMG->num<Face>() > 0) && !pMG->option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		assignSidesOnly = false;

	if(assignSidesOnly){
		MarkAssociatedSides<Volume>(boolMarker);
		MarkAssociatedSides<Face>(boolMarker);
		MarkAssociatedSides<EdgeBase>(boolMarker);
	}
	else
	{
		MarkAssociatedLowerDimElems<Volume>(boolMarker);
		MarkAssociatedLowerDimElems<Face>(boolMarker);
		MarkAssociatedLowerDimElems<EdgeBase>(boolMarker);
	}

#ifdef UG_PARALLEL
	RemoveSurfaceViewMarks<VertexBase>(boolMarker, pMG, distGridMgr);
	RemoveSurfaceViewMarks<EdgeBase>(boolMarker, pMG, distGridMgr);
	RemoveSurfaceViewMarks<Face>(boolMarker, pMG, distGridMgr);
	RemoveSurfaceViewMarks<Volume>(boolMarker, pMG, distGridMgr);
#else
	RemoveSurfaceViewMarks<VertexBase>(boolMarker, pMG);
	RemoveSurfaceViewMarks<EdgeBase>(boolMarker, pMG);
	RemoveSurfaceViewMarks<Face>(boolMarker, pMG);
	RemoveSurfaceViewMarks<Volume>(boolMarker, pMG);
#endif

#ifdef UG_PARALLEL
//	for a parallel subset handler we have to copy the subset-indices to
//	avoid problems at interfaces.
//	This happens e.g. in 2d, where a triangle is an element of the surface view only
//	on one process. Associated vertices on other processes wouldn't know that
//	they are surface view element, too. This has to be communicated.
	ComPol_BoolMarker_AddMarks<VertexLayout> cpSubsetVRT(boolMarker);
	ComPol_BoolMarker_AddMarks<EdgeLayout> cpSubsetEDGE(boolMarker);
	ComPol_BoolMarker_AddMarks<FaceLayout> cpSubsetFACE(boolMarker);
	ComPol_BoolMarker_AddMarks<VolumeLayout> cpSubsetVOL(boolMarker);
	pcl::InterfaceCommunicator<VertexLayout> comVRT;
	pcl::InterfaceCommunicator<EdgeLayout> comEDGE;
	pcl::InterfaceCommunicator<FaceLayout> comFACE;
	pcl::InterfaceCommunicator<VolumeLayout> comVOL;

	comVRT.send_data(distGridMgr.grid_layout_map().get_layout<VertexBase>(INT_H_SLAVE), cpSubsetVRT);
	comEDGE.send_data(distGridMgr.grid_layout_map().get_layout<EdgeBase>(INT_H_SLAVE), cpSubsetEDGE);
	comFACE.send_data(distGridMgr.grid_layout_map().get_layout<Face>(INT_H_SLAVE), cpSubsetFACE);
	comVOL.send_data(distGridMgr.grid_layout_map().get_layout<Volume>(INT_H_SLAVE), cpSubsetVOL);

	comVRT.receive_data(distGridMgr.grid_layout_map().get_layout<VertexBase>(INT_H_MASTER), cpSubsetVRT);
	comEDGE.receive_data(distGridMgr.grid_layout_map().get_layout<EdgeBase>(INT_H_MASTER), cpSubsetEDGE);
	comFACE.receive_data(distGridMgr.grid_layout_map().get_layout<Face>(INT_H_MASTER), cpSubsetFACE);
	comVOL.receive_data(distGridMgr.grid_layout_map().get_layout<Volume>(INT_H_MASTER), cpSubsetVOL);

	comVRT.communicate();
	comEDGE.communicate();
	comFACE.communicate();
	comVOL.communicate();

	comVRT.send_data(distGridMgr.grid_layout_map().get_layout<VertexBase>(INT_H_MASTER), cpSubsetVRT);
	comEDGE.send_data(distGridMgr.grid_layout_map().get_layout<EdgeBase>(INT_H_MASTER), cpSubsetEDGE);
	comFACE.send_data(distGridMgr.grid_layout_map().get_layout<Face>(INT_H_MASTER), cpSubsetFACE);
	comVOL.send_data(distGridMgr.grid_layout_map().get_layout<Volume>(INT_H_MASTER), cpSubsetVOL);

	comVRT.receive_data(distGridMgr.grid_layout_map().get_layout<VertexBase>(INT_H_SLAVE), cpSubsetVRT);
	comEDGE.receive_data(distGridMgr.grid_layout_map().get_layout<EdgeBase>(INT_H_SLAVE), cpSubsetEDGE);
	comFACE.receive_data(distGridMgr.grid_layout_map().get_layout<Face>(INT_H_SLAVE), cpSubsetFACE);
	comVOL.receive_data(distGridMgr.grid_layout_map().get_layout<Volume>(INT_H_SLAVE), cpSubsetVOL);

	comVRT.communicate();
	comEDGE.communicate();
	comFACE.communicate();
	comVOL.communicate();
#endif
}

////////////////////////////////////////////////////////////////////////////////
// SurfaceView
////////////////////////////////////////////////////////////////////////////////

SurfaceView::SurfaceView(SmartPtr<MGSubsetHandler> spMGSH,
#ifdef UG_PARALLEL
                         DistributedGridManager* pDistGridMgr,
#endif
                         bool adaptiveMG) :
	m_spMGSH(spMGSH),
	m_adaptiveMG(adaptiveMG),
	m_pMG(m_spMGSH->multi_grid())
#ifdef UG_PARALLEL
	,m_pDistGridMgr(pDistGridMgr)
#endif
{
	UG_ASSERT(m_pMG, "A MultiGrid has to be assigned to the given subset handler");

	m_Marker.assign_grid(m_pMG);

	mark_shadows();
}


SurfaceLevelView::
SurfaceLevelView(SmartPtr<SurfaceView> spSV, int topLvl) :
	m_spSV(spSV),
	m_topLvl(topLvl)
{
}

template <>
bool SurfaceView::is_shadowed(GeometricObject* obj) const
{
	switch(obj->base_object_id())
	{
		case VERTEX: return is_shadowed(static_cast<VertexBase*>(obj));
		case EDGE: return is_shadowed(static_cast<EdgeBase*>(obj));
		case FACE: return is_shadowed(static_cast<Face*>(obj));
		case VOLUME: return is_shadowed(static_cast<Volume*>(obj));
		default: UG_THROW("Base Object type not found.");
	}
}

void SurfaceView::mark_shadows()
{
#ifdef UG_PARALLEL
//	get multigrid
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(m_pDistGridMgr->get_assigned_grid());
	if(!pMG) throw(UGError("  Can't create surface-view. A Multigrid is required.\n"));

	MarkShadows(m_Marker, pMG, *m_pDistGridMgr);
#else
	MarkShadows(m_Marker, m_pMG);
#endif
}

}// end of namespace
