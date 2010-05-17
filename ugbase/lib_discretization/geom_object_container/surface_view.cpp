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
	enable_subset_attachments(true);
}

SurfaceView::SurfaceView(MultiGrid& mg) : SubsetHandler(mg)
{
	m_pMG = &mg;
	enable_subset_attachments(true);
}

SurfaceView::~SurfaceView()
{
//	unregister all observers
	while(!m_gridObservers.empty())
	{
		unregister_observer(m_gridObservers.back());
	}
}

////////////////////////////////////////////////////////////////////////
void SurfaceView::assign_grid(MultiGrid& mg)
{
	m_pMG = &mg;
	SubsetHandler::assign_grid(mg);
}

void SurfaceView::registered_at_grid(Grid* grid)
{
	assert(static_cast<MultiGrid*>(grid) == m_pMG && "you should assign a grid using SurfaceView::assign_grid(...)");
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
/*
//	if the observer is a grid observer, notify him about the registration
	if((observerType & OT_GRID_OBSERVER) == OT_GRID_OBSERVER)
		observer->registered_at_grid(this);
*/
}

void SurfaceView::unregister_observer(GridObserver* observer)
{
//	check where the observer has been registered and erase the corresponding entries.
	bool unregisterdFromGridObservers = false;

	{
		ObserverContainer::iterator iter = find(m_gridObservers.begin(),
												m_gridObservers.end(), observer);
		if(iter != m_gridObservers.end())
			m_gridObservers.erase(iter);

		unregisterdFromGridObservers = true;
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

//	if the observer is a grid observer, notify him about the unregistration
	if(unregisterdFromGridObservers)
		observer->unregistered_from_grid(NULL);
}

}//	end of namespace

