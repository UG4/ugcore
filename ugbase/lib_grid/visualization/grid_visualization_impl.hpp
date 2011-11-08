// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// Nov 2011

#ifndef __H__UG__grid_visualization_impl__
#define __H__UG__grid_visualization_impl__

#include "grid_visualization.h"

namespace ug{

/**	Calculate the normal of a face. If the position accessor is not for
 * three dimensions, then (0, 0, 1) is returned.
 * \{ */
static void CalculateNormalForVisualization(vector3& nOut, Face* f,
						Grid::AttachmentAccessor<VertexBase, APosition>& aaPos)
{
	CalculateNormal(nOut, f, aaPos);
}

static void CalculateNormalForVisualization(vector3& nOut, Face* f,
						Grid::AttachmentAccessor<VertexBase, APosition2>& aaPos)
{
	nOut = vector3(0, 0, 1);
}

static void CalculateNormalForVisualization(vector3& nOut, Face* f,
						Grid::AttachmentAccessor<VertexBase, APosition1>& aaPos)
{
	nOut = vector3(0, 0, 1);
}
/**	\} */


template <class TNumber, class TInt, class TAVrtPos>
GridVisualization<TNumber, TInt, TAVrtPos>::
GridVisualization() :
	m_pGrid(NULL)
{
}

template <class TNumber, class TInt, class TAVrtPos>
GridVisualization<TNumber, TInt, TAVrtPos>::
~GridVisualization()
{
	release_grid();
}

template <class TNumber, class TInt, class TAVrtPos>
void GridVisualization<TNumber, TInt, TAVrtPos>::
release_grid()
{
	if(m_pGrid){
		if(m_pGrid->has_vertex_attachment(m_aVrtInd))
			m_pGrid->detach_from_vertices(m_aVrtInd);
		if(m_pGrid->has_face_attachment(m_aNormal))
			m_pGrid->detach_from_faces(m_aNormal);
		m_pGrid = NULL;
	}
}

template <class TNumber, class TInt, class TAVrtPos>
void GridVisualization<TNumber, TInt, TAVrtPos>::
set_grid(Grid& grid, TAVrtPos& aPos)
{
	release_grid();
	m_pGrid = &grid;
	m_aVrtPos = aPos;

//	access attachments
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	grid.attach_to_vertices(m_aVrtInd);
	grid.attach_to_faces(m_aNormal);

	m_aaIndVRT.access(grid, m_aVrtInd);
	m_aaPosVRT.access(grid, aPos);
	m_aaNormFACE.access(grid, m_aNormal);
}

template <class TNumber, class TInt, class TAVrtPos>
void GridVisualization<TNumber, TInt, TAVrtPos>::
update_visuals()
{
	using namespace std;

	if(!m_pGrid){
		UG_THROW("A grid is required in GridVisualization::update_geometry.");
	}
	
	Grid& g = *m_pGrid;
	
//	prepare vertex data
	vector<TNumber>& positions = m_vrtPositions;
	positions.resize(g.num<VertexBase>() * 3);


//	copy position data and set vertex indices
	size_t counter = 0;
	if(TAVrtPos::ValueType::Size == 1){
		for(VertexBaseIterator iter = g.begin<VertexBase>();
			iter != g.end<VertexBase>(); ++iter, counter += 3)
		{
			m_aaIndVRT[*iter] = counter / 3;
			positions[counter] = (TNumber)m_aaPosVRT[*iter][0];
			positions[counter+1] = 0;
			positions[counter+2] = 0;
		}
	}
	else if(TAVrtPos::ValueType::Size == 2){
		for(VertexBaseIterator iter = g.begin<VertexBase>();
			iter != g.end<VertexBase>(); ++iter, counter += 3)
		{
			m_aaIndVRT[*iter] = counter / 3;
			positions[counter] = (TNumber)m_aaPosVRT[*iter][0];
			positions[counter+1] = (TNumber)m_aaPosVRT[*iter][1];
			positions[counter+2] = 0;
		}
	}
	else if(TAVrtPos::ValueType::Size == 3){
		for(VertexBaseIterator iter = g.begin<VertexBase>();
			iter != g.end<VertexBase>(); ++iter, counter += 3)
		{
			m_aaIndVRT[*iter] = counter / 3;
			positions[counter] = (TNumber)m_aaPosVRT[*iter][0];
			positions[counter+1] = (TNumber)m_aaPosVRT[*iter][1];
			positions[counter+2] = (TNumber)m_aaPosVRT[*iter][2];
		}
	}


//	currently we only have one visual
	while(m_visuals.size() < 1){
		m_visuals.push_back(SPVisual(new Visual));
	}


//	access visual and initialize values
	SPVisual vis = m_visuals[0];
	vis->m_visType = E_GRID_SUBSET;

//	resize the arrays in the visual
	vector<TNumber>& faceNormals = vis->m_faceNormals;
	vector<TInt>& tris = vis->m_triIndices;
	TNumber* color = vis->m_color;

//todo: take color from subset-info
	color[0] = 1.f;
	color[1] = 1.f;
	color[2] = 1.f;
	color[3] = 1.f;

//todo:support quadrilaterals
	faceNormals.resize(g.num<Triangle>() * 3);
	tris.resize(g.num<Triangle>() * 3);


//	fill triangle list
	counter = 0;
	for(TriangleIterator iter = g.begin<Triangle>();
		iter != g.end<Triangle>(); ++iter, counter += 3)
	{
		Triangle* t = *iter;
	//	indices
		tris[counter] = (TInt)m_aaIndVRT[t->vertex(0)];
		tris[counter+1] = (TInt)m_aaIndVRT[t->vertex(1)];
		tris[counter+2] = (TInt)m_aaIndVRT[t->vertex(2)];

	}

//todo:	quadrilateral lists


//	calculate face normals
	vector3 n;

//	if we have a three dimensional position attachment, we can calculate the normals
//	directly from that one. If not, we assume that they point upwards (in z-direction)
//todo: do this for all faces.
	counter = 0;
	for(TriangleIterator iter = g.begin<Triangle>();
		iter != g.end<Triangle>(); ++iter, counter += 3)
	{
		Triangle* t = *iter;

		CalculateNormalForVisualization(n, t, m_aaPosVRT);
		faceNormals[counter] = (TNumber)n[0];
		faceNormals[counter+1] = (TNumber)n[1];
		faceNormals[counter+2] = (TNumber)n[2];
	}
}


template <class TNumber, class TInt, class TAVrtPos>
int GridVisualization<TNumber, TInt, TAVrtPos>::
num_vertices()
{
	return (int)m_vrtPositions.size() / 3;
}

template <class TNumber, class TInt, class TAVrtPos>
const TNumber* GridVisualization<TNumber, TInt, TAVrtPos>::
vertex_positions()
{
	UG_ASSERT(!m_vrtPositions.empty(), "No vertex data available.");
	return &m_vrtPositions.front();
}


template <class TNumber, class TInt, class TAVrtPos>
int GridVisualization<TNumber, TInt, TAVrtPos>::
num_visuals()
{
	return (int)m_visuals.size();
}

template <class TNumber, class TInt, class TAVrtPos>
int GridVisualization<TNumber, TInt, TAVrtPos>::
visual_type(int visInd)
{
	UG_ASSERT(m_visuals.size() > (size_t)visInd, "Bad index!");
	return m_visuals[visInd]->m_visType;
}

template <class TNumber, class TInt, class TAVrtPos>
const TNumber* GridVisualization<TNumber, TInt, TAVrtPos>::
visual_color(int visInd)
{
	UG_ASSERT(m_visuals.size() > (size_t)visInd, "Bad index!");
	return m_visuals[visInd]->m_color;
}

template <class TNumber, class TInt, class TAVrtPos>
const TNumber* GridVisualization<TNumber, TInt, TAVrtPos>::
face_normals(int visInd)
{
	UG_ASSERT(m_visuals.size() > (size_t)visInd, "Bad index!");
	return &m_visuals[visInd]->m_faceNormals.front();
}

template <class TNumber, class TInt, class TAVrtPos>
int GridVisualization<TNumber, TInt, TAVrtPos>::
num_triangles(int visInd)
{
	UG_ASSERT(m_visuals.size() > (size_t)visInd, "Bad index!");
	return m_visuals[visInd]->m_triIndices.size() / 3;
}

template <class TNumber, class TInt, class TAVrtPos>
const TInt* GridVisualization<TNumber, TInt, TAVrtPos>::
triangle_list(int visInd)
{
	UG_ASSERT(m_visuals.size() > (size_t)visInd, "Bad index!");
	return &m_visuals[visInd]->m_triIndices.front();
}

}//	end of namespace

#endif
