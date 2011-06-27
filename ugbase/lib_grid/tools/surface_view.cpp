//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d24

#include <cassert>
#include "surface_view.h"

#define NOTIFY_OBSERVERS(observerContainer, callback)	{for(ObserverContainer::iterator iter = observerContainer.begin(); iter != observerContainer.end(); iter++) (*iter)->callback;}

namespace ug
{
////////////////////////////////////////////////////////////////////////
SurfaceView::SurfaceView() : SubsetHandler()
{
//	enable_subset_attachments(true);
	enable_subset_inheritance(false);
}

SurfaceView::SurfaceView(MultiGrid& mg) : SubsetHandler(mg)
{
	m_pMG = &mg;
//	enable_subset_attachments(true);
	enable_subset_inheritance(false);
}

SurfaceView::~SurfaceView()
{
//	tell registered grid-observers that the SurfaceView is to be destroyed.
	for(ObserverContainer::iterator iter = m_gridObservers.begin();
		iter != m_gridObservers.end(); ++iter)
	{
	//	we have to pass the underlying grid here, since SurfaceViews are
	//	not directly supported by GridObservers.
		(*iter)->grid_to_be_destroyed(m_pGrid);
	}
	
//	unregister all observers
	while(!m_gridObservers.empty())
	{
		unregister_observer(m_gridObservers.back());
	}
}

void SurfaceView::grid_to_be_destroyed(Grid* grid)
{
//	notify all observers that the underlying grid is to be destroyed.
	SubsetHandler::grid_to_be_destroyed(grid);
}

////////////////////////////////////////////////////////////////////////
void SurfaceView::assign_grid(MultiGrid& mg)
{
	m_pMG = &mg;
	SubsetHandler::assign_grid(mg);
}

////////////////////////////////////////////////////////////////////////
void SurfaceView::assign_subset(VertexBase* elem, int subsetIndex)
{
	int oldInd = get_subset_index(elem);
	if(oldInd >= 0 && subsetIndex == -1){
		NOTIFY_OBSERVERS(m_vertexObservers, vertex_to_be_erased(m_pGrid, elem));
	}
		
	SubsetHandler::assign_subset(elem, subsetIndex);
	
	if(oldInd == -1 && subsetIndex >= 0){
		NOTIFY_OBSERVERS(m_vertexObservers, vertex_created(m_pGrid, elem, NULL));
	}
}

void SurfaceView::assign_subset(EdgeBase* elem, int subsetIndex)
{
	int oldInd = get_subset_index(elem);
	if(oldInd >= 0 && subsetIndex == -1){
		NOTIFY_OBSERVERS(m_edgeObservers, edge_to_be_erased(m_pGrid, elem));
	}
		
	SubsetHandler::assign_subset(elem, subsetIndex);
	
	if(oldInd == -1 && subsetIndex >= 0){
		NOTIFY_OBSERVERS(m_edgeObservers, edge_created(m_pGrid, elem, NULL));
	}
}

void SurfaceView::assign_subset(Face* elem, int subsetIndex)
{
	int oldInd = get_subset_index(elem);
	if(oldInd >= 0 && subsetIndex == -1){
		NOTIFY_OBSERVERS(m_faceObservers, face_to_be_erased(m_pGrid, elem));
	}
		
	SubsetHandler::assign_subset(elem, subsetIndex);
	
	if(oldInd == -1 && subsetIndex >= 0){
		NOTIFY_OBSERVERS(m_faceObservers, face_created(m_pGrid, elem, NULL));
	}
}

void SurfaceView::assign_subset(Volume* elem, int subsetIndex)
{
	int oldInd = get_subset_index(elem);
	if(oldInd >= 0 && subsetIndex == -1){
		NOTIFY_OBSERVERS(m_volumeObservers, volume_to_be_erased(m_pGrid, elem));
	}
		
	SubsetHandler::assign_subset(elem, subsetIndex);
	
	if(oldInd == -1 && subsetIndex >= 0){
		NOTIFY_OBSERVERS(m_volumeObservers, volume_created(m_pGrid, elem, NULL));
	}
}

////////////////////////////////////////////////////////////////////////
void SurfaceView::register_observer(GridObserver* observer, uint observerType)
{
//	check which elements have to be observed and store pointers to the observers.
//	avoid double-registration!
	if((observerType & OT_GRID_OBSERVER) == OT_GRID_OBSERVER)
	{
		ObserverContainer::iterator iter = find(m_gridObservers.begin(),
												m_gridObservers.end(), observer);
		if(iter == m_gridObservers.end())
			m_gridObservers.push_back(observer);
	}

	if((observerType & OT_VERTEX_OBSERVER) == OT_VERTEX_OBSERVER)
	{
		ObserverContainer::iterator iter = find(m_vertexObservers.begin(),
												m_vertexObservers.end(), observer);
		if(iter == m_vertexObservers.end())
			m_vertexObservers.push_back(observer);
	}

	if((observerType & OT_EDGE_OBSERVER) == OT_EDGE_OBSERVER)
	{
		ObserverContainer::iterator iter = find(m_edgeObservers.begin(),
												m_edgeObservers.end(), observer);
		if(iter == m_edgeObservers.end())
			m_edgeObservers.push_back(observer);
	}

	if((observerType & OT_FACE_OBSERVER) == OT_FACE_OBSERVER)
	{
		ObserverContainer::iterator iter = find(m_faceObservers.begin(),
												m_faceObservers.end(), observer);
		if(iter == m_faceObservers.end())
			m_faceObservers.push_back(observer);
	}

	if((observerType & OT_VOLUME_OBSERVER) == OT_VOLUME_OBSERVER)
	{
		ObserverContainer::iterator iter = find(m_volumeObservers.begin(),
												m_volumeObservers.end(), observer);
		if(iter == m_volumeObservers.end())
			m_volumeObservers.push_back(observer);
	}
}

void SurfaceView::unregister_observer(GridObserver* observer)
{
//	check where the observer has been registered and erase the corresponding entries.

	{
		ObserverContainer::iterator iter = find(m_gridObservers.begin(),
												m_gridObservers.end(), observer);
		if(iter != m_gridObservers.end())
			m_gridObservers.erase(iter);
	}

	{
		ObserverContainer::iterator iter = find(m_vertexObservers.begin(),
												m_vertexObservers.end(), observer);
		if(iter != m_vertexObservers.end())
			m_vertexObservers.erase(iter);
	}

	{
		ObserverContainer::iterator iter = find(m_edgeObservers.begin(),
												m_edgeObservers.end(), observer);
		if(iter != m_edgeObservers.end())
			m_edgeObservers.erase(iter);
	}

	{
		ObserverContainer::iterator iter = find(m_faceObservers.begin(),
												m_faceObservers.end(), observer);
		if(iter != m_faceObservers.end())
			m_faceObservers.erase(iter);
	}

	{
		ObserverContainer::iterator iter = find(m_volumeObservers.begin(),
												m_volumeObservers.end(), observer);
		if(iter != m_volumeObservers.end())
			m_volumeObservers.erase(iter);
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Check Surface Grid
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/// check if all elements are assigned correctly to the surface view
template <typename TElem>
bool CheckSurfaceViewElements(const SurfaceView& surfView,
                                     const MultiGrid& mg)
{
	typename geometry_traits<TElem>::const_iterator iter;

	std::vector<EdgeBase*> vEdge;
	std::vector<Face*> vFace;
	std::vector<Volume*> vVolume;

	for(size_t lev = 0; lev < mg.num_levels(); ++lev)
		for(iter = mg.begin<TElem>(lev); iter != mg.end<TElem>(lev); ++iter)
		{
		//	get element
			TElem* elem = *iter;

		//	get subset index
			const int si = surfView.get_subset_index(elem);

		//	check if in surface view
			bool inSurfView = (si >= 0);

			if(!mg.has_children(elem))
			{
				if(!inSurfView)
				{
					UG_LOG("ERROR in CheckSurfaceView: Element on level "
							<< lev << " without child is"
							" not contained in SurfaceView.\n");
					return false;
				}
			}
			else
			{
			//	collect all associated elements
				CollectAssociated(vEdge, *const_cast<MultiGrid*>(&mg), elem);
				CollectAssociated(vFace, *const_cast<MultiGrid*>(&mg), elem);
				CollectAssociated(vVolume, *const_cast<MultiGrid*>(&mg), elem);

				bool bIsShadow = false;

			//	search for associated surface grid edge
				for(size_t i = 0; i < vEdge.size(); ++i)
					if(!mg.has_children(vEdge[i]))
						bIsShadow = true;

			//	search for associated surface grid face
				for(size_t i = 0; i < vFace.size(); ++i)
					if(!mg.has_children(vFace[i]))
						bIsShadow = true;

			//	search for associated surface grid volume
				for(size_t i = 0; i < vVolume.size(); ++i)
					if(!mg.has_children(vVolume[i]))
						bIsShadow = true;

				if(bIsShadow && !inSurfView)
				{
					UG_LOG("ERROR in CheckSurfaceView: Shadow element on level "
							<< lev << " is not contained in SurfaceView.\n");
					return false;
				}

				if(!bIsShadow && inSurfView)
				{
					UG_LOG("ERROR in CheckSurfaceView: Element in Surface view "
							"on level "	<< lev <<
							" though has children and is no shadow.\n");
					return false;
				}
			}
		}

//	all fine
	return true;
}

/// checks if surface view is correct
bool CheckSurfaceView(const SurfaceView& surfView)
{
//	get underlying grid
	Grid* grid = surfView.get_assigned_grid();

//	check that grid is a MultiGrid
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(grid);
	if(pMG == NULL)
	{
		UG_LOG("ERROR in CheckSurfaceView: underlying grid not a Multigrid.\n");
		return false;
	}

//	check all elements
	if(!CheckSurfaceViewElements<VertexBase>(surfView, *pMG))
	{
		UG_LOG("ERROR in CheckSurfaceView: wrong VertexBase found.\n");
		return false;
	}
	if(!CheckSurfaceViewElements<EdgeBase>(surfView, *pMG))
	{
		UG_LOG("ERROR in CheckSurfaceView: wrong EdgeBase found.\n");
		return false;
	}
	if(!CheckSurfaceViewElements<Face>(surfView, *pMG))
	{
		UG_LOG("ERROR in CheckSurfaceView: wrong Face found.\n");
		return false;
	}
	if(!CheckSurfaceViewElements<Volume>(surfView, *pMG))
	{
		UG_LOG("ERROR in CheckSurfaceView: wrong Volume found.\n");
		return false;
	}

//	everything ok
	return true;
}


}//	end of namespace

