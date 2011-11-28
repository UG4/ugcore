// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)
 
#include "surface_view.h"

namespace ug{
namespace tmp{




#ifdef OUT_OF_ORDER
SurfaceView::SurfaceView() :
	m_pMG(NULL),
	m_aSharedEntry("MGSubsetHandler_SharedListEntry", false)
{
}

SurfaceView::SurfaceView(MultiGrid& mg) :
	m_pMG(NULL),
	m_aSharedEntry("MGSubsetHandler_SharedListEntry", false)
{
	assign_grid(mg);
}

void SurfaceView::grid_to_be_destroyed(Grid* grid)
{
	cleanup();
	ISubsetHandler::grid_to_be_destroyed(grid);
}

void SurfaceView::cleanup()
{
	//erase_subset_lists();
	//m_levels.clear();

	if(m_pMG){
		m_pMG->detach_from_vertices(m_aSharedEntry);
		m_pMG->detach_from_edges(m_aSharedEntry);
		m_pMG->detach_from_faces(m_aSharedEntry);
		m_pMG->detach_from_volumes(m_aSharedEntry);

		m_pMG = NULL;
	}

	ISubsetHandler::set_grid(NULL);
}

void SurfaceView::assign_grid(MultiGrid& mg)
{
	if(m_pMG)
		cleanup();

	m_pMG = &mg;

	ISubsetHandler::set_grid(&mg);

//	attach shared entries
	m_pGrid->attach_to_vertices(m_aSharedEntry);
	m_pGrid->attach_to_edges(m_aSharedEntry);
	m_pGrid->attach_to_faces(m_aSharedEntry);
	m_pGrid->attach_to_volumes(m_aSharedEntry);
}

void SurfaceView::erase_subset_lists()
{
}

void SurfaceView::clear_subset_lists(int index)
{
}

void SurfaceView::
change_subset_indices(int indOld, int indNew)
{
}

void SurfaceView::add_required_subset_lists(int maxIndex)
{
}

void SurfaceView::erase_subset_lists(int index)
{
}

void SurfaceView::swap_subset_lists(int ind1, int ind2)
{
}

void SurfaceView::move_subset_lists(int indexFrom, int indexTo)
{
}

void SurfaceView::
register_subset_elements_at_pipe()
{
}
#endif

}// end of namespace
}// end of namespace
