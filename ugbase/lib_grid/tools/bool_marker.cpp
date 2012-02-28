// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.01.2012 (m,d,y)
 
#include "bool_marker.h"

namespace ug{

BoolMarker::BoolMarker() :
	m_pGrid(NULL)
{
}

BoolMarker::BoolMarker(Grid& g) :
	m_pGrid(NULL)
{
	assign_grid(&g);
}

BoolMarker::~BoolMarker()
{
	assign_grid(NULL);
}

void BoolMarker::assign_grid(Grid* g)
{
	if(m_pGrid == g)
		return;

	if(m_pGrid){
		m_pGrid->detach_from_all(m_aBool);
		m_pGrid->unregister_observer(this);
		m_aaMarkVRT.invalidate();
		m_aaMarkEDGE.invalidate();
		m_aaMarkFACE.invalidate();
		m_aaMarkVOL.invalidate();
	}

	m_pGrid = g;
	if(g){
		g->register_observer(this, OT_GRID_OBSERVER);
		g->attach_to_all_dv(m_aBool, false);
		m_aaMarkVRT.access(*g, m_aBool);
		m_aaMarkEDGE.access(*g, m_aBool);
		m_aaMarkFACE.access(*g, m_aBool);
		m_aaMarkVOL.access(*g, m_aBool);
	}
}

void BoolMarker::grid_to_be_destroyed(Grid* grid)
{
	assign_grid(NULL);
}

void BoolMarker::clear()
{
	assert(m_pGrid);
	unmark(m_pGrid->begin<VertexBase>(), m_pGrid->end<VertexBase>());
	unmark(m_pGrid->begin<EdgeBase>(), m_pGrid->end<EdgeBase>());
	unmark(m_pGrid->begin<Face>(), m_pGrid->end<Face>());
	unmark(m_pGrid->begin<Volume>(), m_pGrid->end<Volume>());
}

}// end of namespace
