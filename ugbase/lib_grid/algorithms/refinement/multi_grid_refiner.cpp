//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d09

#include <cassert>
#include "multi_grid_refiner.h"
#include "lib_grid/algorithms/algorithms.h"

using namespace std;

namespace ug
{

MultiGridRefiner::MultiGridRefiner()
{
	m_pMG = NULL;
}

MultiGridRefiner::MultiGridRefiner(MultiGrid& mg)
{
	m_pMG = NULL;
	assign_grid(mg);
}

void MultiGridRefiner::assign_grid(MultiGrid& mg)
{
	if(m_pMG)
	{
		m_pMG->unregister_observer(this);
		m_pMG = NULL;
	}

	mg.register_observer(this, OT_GRID_OBSERVER);
}

void MultiGridRefiner::registered_at_grid(Grid* grid)
{
	if(m_pMG)
		m_pMG->unregister_observer(this);

	m_pMG = dynamic_cast<MultiGrid*>(grid);
	assert(m_pMG && "MultiGridRefiner registered at incopatible grid.");

	if(!m_pMG)
		return;

	m_selMarks.assign_grid(*m_pMG);
	m_selMarks.enable_autoselection(false);
	m_selMarks.enable_selection_inheritance(false);
}

void MultiGridRefiner::unregistered_from_grid(Grid* grid)
{
	if(m_pMG)
	{
		m_pMG->unregister_observer(&m_selMarks);
		m_pMG = NULL;
	}
}

void MultiGridRefiner::clear_marks()
{
	m_selMarks.clear_selection();
}

////////////////////////////////////////////////////////////////////////
void MultiGridRefiner::refine()
{
	assert(m_pMG && "refiner not has to be assigned to a multi-grid!");
	if(!m_pMG)
		return;

//	the multi-grid
	MultiGrid& mg = *m_pMG;

//	make sure that the required options are enabled.
	if(!mg.option_is_enabled(GRIDOPT_FULL_INTERCONNECTION))
	{
		LOG("WARNING in MultiGridRefiner::refine(): auto-enabling GRIDOPT_FULL_INTERCONNECTION.\n");
		mg.enable_options(GRIDOPT_FULL_INTERCONNECTION);
	}

//	access position attachments
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	if(mg.has_vertex_attachment(aPosition))
		aaPos.access(mg, aPosition);

//	collect objects for refine
	collect_objects_for_refine();

//	cout << "num marked edges: " << m_selMarks.num<EdgeBase>() << endl;
//	cout << "num marked faces: " << m_selMarks.num<Face>() << endl;

//	we want to add new elements in a new layer.
	bool bHierarchicalInsertionWasEnabled = mg.hierarchical_insertion_enabled();
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(true);
//LOG("creating new vertices\n");
//	create new vertices from marked vertices
	for(VertexBaseIterator iter = m_selMarks.begin<VertexBase>();
		iter != m_selMarks.end<VertexBase>(); ++iter)
	{
		VertexBase* v = *iter;
		if(!mg.get_child_vertex(v))
		{
		//	create a new vertex in the next layer.
			VertexBase* nVrt = *mg.create_by_cloning(v, v);
			if(aaPos.valid())
			{
				aaPos[nVrt] = aaPos[v];
			//	change y-coord to visualise the hierarchy
				aaPos[nVrt].z;
			}
		}
	}
//LOG("creating new edges\n");
//	create new vertices and edges from marked edges
	int numNewEdges = 0;
	for(EdgeBaseIterator iter = m_selMarks.begin<EdgeBase>();
		iter != m_selMarks.end<EdgeBase>(); ++iter)
	{
	//	collect_objects_for_refine removed all edges that already were
	//	refined. No need to check that again.
		EdgeBase* e = *iter;
		Vertex* nVrt = *mg.create<Vertex>(e);
		VertexBase* substituteVrts[2];
		substituteVrts[0] = mg.get_child_vertex(e->vertex(0));
		substituteVrts[1] = mg.get_child_vertex(e->vertex(1));

		if(aaPos.valid())
		{
			VecScaleAdd(aaPos[nVrt], 0.5, aaPos[e->vertex(0)],
									0.5, aaPos[e->vertex(1)]);
			aaPos[nVrt].z;
		}

	//	split the edge
		vector<EdgeBase*> vEdges(2);
		e->refine(vEdges, nVrt, substituteVrts);
		assert((vEdges.size() == 2) && "Edge refine produced wrong number of edges.");
		mg.register_element(vEdges[0], e);
		mg.register_element(vEdges[1], e);
	}

//LOG("creating new faces\n");
//	create new vertices and faces from marked faces
	for(FaceIterator iter = m_selMarks.begin<Face>();
		iter != m_selMarks.end<Face>(); ++iter)
	{
	//TODO: improve this
		Face* f = *iter;

		vector<EdgeBase*> 	vEdges(f->num_edges());
		vector<VertexBase*> vNewEdgeVertices(f->num_edges());
		vector<VertexBase*> vNewVertexVertices(f->num_vertices());
		vector<Face*>		vFaces(f->num_vertices());// heuristic
	//	collect all associated edges.
		CollectEdges(vEdges, mg, f);
		uint numEdges = vEdges.size();

		assert(numEdges == f->num_edges() && "associated edges missing.");

	//	each vertex should have an associated vertex in the next level
		for(uint i = 0; i < vNewVertexVertices.size(); ++i)
		{
			assert(mg.get_child_vertex(f->vertex(i)) && "child vertex missing.");
			vNewVertexVertices[i] = mg.get_child_vertex(f->vertex(i));
		}

	//	each edge should have an associated vertex. sort them into vNewEdgeVertices.
		for(uint i = 0; i < vEdges.size(); ++i)
		{
			EdgeBase* e = vEdges[i];
			int edgeIndex = GetEdgeIndex(f, e);

			assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "unknown problem in CollectEdges / GetEdgeIndex.");
			assert((mg.get_child_vertex(e) != NULL) && "no new vertex on refined edge.");
			vNewEdgeVertices[edgeIndex] = mg.get_child_vertex(e);
		}

	//	we'll perform a regular refine
		VertexBase* vNewVrt = NULL;
		f->refine_regular(vFaces, &vNewVrt, vNewEdgeVertices, NULL, Vertex(), &vNewVertexVertices);

	//	if a new vertex has been created during refine, then register it at the grid.
		if(vNewVrt)
		{
			mg.register_element(vNewVrt, f);
		//	assign a new position
			if(aaPos.valid())
			{
				aaPos[vNewVrt] = CalculateCenter(f, aaPos);
				aaPos[vNewVrt].z;
			}
		}

		for(uint i = 0; i < vFaces.size(); ++i)
		{
			Face* newFace = vFaces[i];
			mg.register_element(newFace, f);
		}
	}

//	done - clean up
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(false);

	m_selMarks.clear();
}

void MultiGridRefiner::collect_objects_for_refine()
{
//	make sure, that elements that already have children won't be
//	refined again.
//	all elements that are associated with marked elements and
//	which are of lower dimension have to be marked too.
	vector<EdgeBase*> 	vEdges;
	vector<Face*> 		vFaces;

//	regard the multi-grid as a grid
	Grid& grid =*(Grid*)m_pMG;

//	LOG("num selected vrts: " << m_selMarks.num<VertexBase>() << endl);
//	LOG("num selected edges: " << m_selMarks.num<EdgeBase>() << endl);
//	LOG("num selected faces: " << m_selMarks.num<Face>() << endl);
//	LOG("num selected vols: " << m_selMarks.num<Volume>() << endl);

//	select elements that are associated with volumes
	for(VolumeIterator iter = m_selMarks.begin<Volume>();
		iter != m_selMarks.end<Volume>();)
	{
		Volume* v = *iter;
		++iter;//iterator is incremented since it would be invalidated during deselect.
		if(m_pMG->has_children(v))
		{
		//	the element has already been refined.
		//	don't refine it again.
			m_selMarks.deselect(v);
		}
		else
		{
			CollectFaces(vFaces, grid, v);
			for(uint i = 0; i < vFaces.size(); ++i)
				m_selMarks.select(vFaces[i]);
			CollectEdges(vEdges, grid, v);
			for(uint i = 0; i < vEdges.size(); ++i)
				m_selMarks.select(vEdges[i]);
			for(uint i = 0; i < v->num_vertices(); ++i)
				m_selMarks.select(v->vertex(i));
		}
	}

//	select elements that are associated with faces
	for(FaceIterator iter = m_selMarks.begin<Face>();
		iter != m_selMarks.end<Face>();)
	{
		Face* f = *iter;
		++iter;	//iterator is incremented since it would be invalidated during deselect.
		if(m_pMG->has_children(f))
		{
		//	the element has already been refined.
		//	don't refine it again.
			m_selMarks.deselect(f);
		}
		else
		{
			CollectEdges(vEdges, grid, f);
			for(uint i = 0; i < vEdges.size(); ++i)
				m_selMarks.select(vEdges[i]);
			for(uint i = 0; i < f->num_vertices(); ++i)
				m_selMarks.select(f->vertex(i));
		}
	}

//	select elements that are associated with edges
	for(EdgeBaseIterator iter = m_selMarks.begin<EdgeBase>();
		iter != m_selMarks.end<EdgeBase>();)
	{
		EdgeBase* e = *iter;
		++iter;	//iterator is incremented since it would be invalidated during deselect.
		if(m_pMG->has_children(e))
		{
		//	the element has already been refined.
		//	don't refine it again.
			m_selMarks.deselect(e);
		}
		else
		{
			for(uint i = 0; i < e->num_vertices(); ++i)
				m_selMarks.select(e->vertex(i));
		}
	}
}

}//	end of namespace
