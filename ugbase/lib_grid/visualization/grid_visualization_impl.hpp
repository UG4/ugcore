// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// Nov 2011

#ifndef __H__UG__grid_visualization_impl__
#define __H__UG__grid_visualization_impl__

#include "grid_visualization.h"

namespace ug{

template <class TNumber, class TInt, class TAVrtPos>
GridVisualization::GridVisualization() :
	m_pGrid(NULL)
{
}

template <class TNumber, class TInt, class TAVrtPos>
GridVisualization::~GridVisualization()
{
	release_grid();
}

template <class TNumber, class TInt, class TAVrtPos>
void GridVisualization::release_grid()
{
	if(m_pGrid){
		if(m_pGrid.has_face_attachment(m_aNormal))
			m_pGrid.detach_from_faces(m_aNormal);
		m_pGrid = NULL;
	}
}

template <class TNumber, class TInt, class TAVrtPos>
void GridVisualization::set_grid(Grid& grid, TAVrtPos& aPos)
{
	release_grid();
	m_pGrid = &grid;
	grid.attach_to_faces(m_aNormal);
}

template <class TNumber, class TInt, class TAVrtPos>
void GridVisualization::update_visuals()
{
	if(!m_pGrid){
		UG_THROW("A grid is required in GridVisualization::update_geometry.");
	}
	
	Grid& g = *m_pGrid;
	
//	currently we only have one visual
	while(m_visuals.size() < 1){
		m_visuals.push_back(SPVisual(new Visual));
	}
	
//	resize the arrays in the visual
	vector<TNumber>& positions = m_visuals[0]->m_vrtPositions;
	vector<TNumber>& faceNormals = m_visuals[0]->m_faceNormals;
	vector<TInt>& tris = m_visuals[0]->m_triangles;
	TNumber color[4] = m_visuals[0]->m_color;
	
}

}//	end of namespace

#endif
