//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#include <cassert>
#include "global_multi_grid_refiner.h"
#include "lib_grid/algorithms/algorithms.h"

using namespace std;

namespace ug
{

GlobalMultiGridRefiner::GlobalMultiGridRefiner()
{
	m_pMG = NULL;
}

GlobalMultiGridRefiner::GlobalMultiGridRefiner(MultiGrid& mg)
{
	m_pMG = NULL;
	assign_grid(mg);
}

GlobalMultiGridRefiner::~GlobalMultiGridRefiner()
{
	if(m_pMG)
		m_pMG->unregister_observer(this);
}

void GlobalMultiGridRefiner::grid_to_be_destroyed(Grid* grid)
{
	if(m_pMG)
		assign_grid(NULL);
}

void GlobalMultiGridRefiner::assign_grid(MultiGrid& mg)
{
	assign_grid(&mg);
}

void GlobalMultiGridRefiner::assign_grid(MultiGrid* mg)
{
	if(m_pMG){
		m_pMG->unregister_observer(this);
		m_pMG = NULL;
	}
	
	if(mg){
		m_pMG = mg;
		m_pMG->register_observer(this, OT_GRID_OBSERVER);
	}
}

////////////////////////////////////////////////////////////////////////
void GlobalMultiGridRefiner::refine()
{
	assert(m_pMG && "refiner has to be assigned to a multi-grid!");
	if(!m_pMG)
		return;

//	the multi-grid
	MultiGrid& mg = *m_pMG;

//	make sure that the required options are enabled.
	if(!mg.option_is_enabled(GRIDOPT_FULL_INTERCONNECTION))
	{
		LOG("WARNING in GlobalMultiGridRefiner::refine(): auto-enabling GRIDOPT_FULL_INTERCONNECTION.\n");
		mg.enable_options(GRIDOPT_FULL_INTERCONNECTION);
	}

//	access position attachments
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	if(mg.has_vertex_attachment(aPosition))
		aaPos.access(mg, aPosition);

//	notify derivates that refinement begins
	refinement_step_begins();
	
//	cout << "num marked edges: " << m_selMarks.num<EdgeBase>() << endl;
//	cout << "num marked faces: " << m_selMarks.num<Face>() << endl;

//	we want to add new elements in a new layer.
	bool bHierarchicalInsertionWasEnabled = mg.hierarchical_insertion_enabled();
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(true);


//	some buffers
	vector<VertexBase*> vVrts;
	vector<VertexBase*> vEdgeVrts;
	vector<VertexBase*> vFaceVrts;
	vector<EdgeBase*>	vEdges;
	vector<Face*>		vFaces;
	vector<Volume*>		vVols;
	
//	some repeatedly used objects
	EdgeDescriptor ed;
	FaceDescriptor fd;
	VolumeDescriptor vd;

//	the old top level
	int oldTopLevel = mg.num_levels() - 1;
	
//LOG("creating new vertices\n");
//	create new vertices from marked vertices
	for(VertexBaseIterator iter = mg.begin<VertexBase>(oldTopLevel);
		iter != mg.end<VertexBase>(oldTopLevel); ++iter)
	{
		VertexBase* v = *iter;
	//	create a new vertex in the next layer.
		VertexBase* nVrt = *mg.create_by_cloning(v, v);

		if(aaPos.valid())
		{
			aaPos[nVrt] = aaPos[v];
		//	change z-coord to visualise the hierarchy
			//aaPos[nVrt].z += 0.01;
		}
	}

//LOG("creating new edges\n");
//	create new vertices and edges from marked edges
	for(EdgeBaseIterator iter = mg.begin<EdgeBase>(oldTopLevel);
		iter != mg.end<EdgeBase>(oldTopLevel); ++iter)
	{
	//	collect_objects_for_refine removed all edges that already were
	//	refined. No need to check that again.
		EdgeBase* e = *iter;

	//	create two new edges by edge-split
		Vertex* nVrt = *mg.create<Vertex>(e);
		VertexBase* substituteVrts[2];
		substituteVrts[0] = mg.get_child_vertex(e->vertex(0));
		substituteVrts[1] = mg.get_child_vertex(e->vertex(1));

		if(aaPos.valid())
		{
			VecScaleAdd(aaPos[nVrt], 0.5, aaPos[e->vertex(0)],
									0.5, aaPos[e->vertex(1)]);
		//	change z-coord to visualise the hierarchy
			//aaPos[nVrt].z += 0.01;
		}

	//	split the edge
		e->refine(vEdges, nVrt, substituteVrts);
		assert((vEdges.size() == 2) && "Edge refine produced wrong number of edges.");
		mg.register_element(vEdges[0], e);
		mg.register_element(vEdges[1], e);
	}

//LOG("creating new faces\n");
//	create new vertices and faces from marked faces
	for(FaceIterator iter = mg.begin<Face>(oldTopLevel);
		iter != mg.end<Face>(oldTopLevel); ++iter)
	{
		Face* f = *iter;
	//	collect child-vertices
		vVrts.clear();
		for(uint j = 0; j < f->num_vertices(); ++j)
			vVrts.push_back(mg.get_child_vertex(f->vertex(j)));

	//	collect the associated edges
		vEdgeVrts.clear();
		//bool bIrregular = false;
		for(uint j = 0; j < f->num_edges(); ++j)
			vEdgeVrts.push_back(mg.get_child_vertex(mg.get_edge(f, j)));

		VertexBase* newVrt;
		if(f->refine(vFaces, &newVrt, &vEdgeVrts.front(), NULL, &vVrts.front())){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				mg.register_element(newVrt, f);
				if(aaPos.valid()){
					aaPos[newVrt] = CalculateCenter(f, aaPos);
				//	change z-coord to visualise the hierarchy
					//aaPos[newVrt].z += 0.01;
				}
			}

		//	register the new faces and assign status
			for(size_t j = 0; j < vFaces.size(); ++j)
				mg.register_element(vFaces[j], f);
		}
		else{
			LOG("  WARINING in Refine: could not refine face.\n");
		}
	}

//	create new vertices and volumes from marked volumes
	for(VolumeIterator iter = mg.begin<Volume>(oldTopLevel);
		iter != mg.end<Volume>(oldTopLevel); ++iter)
	{
		Volume* v = *iter;
	//	collect child-vertices
		vVrts.clear();
		for(uint j = 0; j < v->num_vertices(); ++j)
			vVrts.push_back(mg.get_child_vertex(v->vertex(j)));

	//	collect the associated edges
		vEdgeVrts.clear();
		//bool bIrregular = false;
		for(uint j = 0; j < v->num_edges(); ++j)
			vEdgeVrts.push_back(mg.get_child_vertex(mg.get_edge(v, j)));

	//	collect associated face-vertices
		vFaceVrts.clear();
		for(uint j = 0; j < v->num_faces(); ++j)
			vFaceVrts.push_back(mg.get_child_vertex(mg.get_face(v, j)));
		
		VertexBase* newVrt;
		if(v->refine(vVols, &newVrt, &vEdgeVrts.front(), &vFaceVrts.front(),
					NULL, Vertex(), &vVrts.front())){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				mg.register_element(newVrt, v);
				if(aaPos.valid()){
					aaPos[newVrt] = CalculateCenter(v, aaPos);
				//	change z-coord to visualise the hierarchy
					//aaPos[newVrt].z += 0.01;
				}
			}

		//	register the new faces and assign status
			for(size_t j = 0; j < vVols.size(); ++j)
				mg.register_element(vVols[j], v);
		}
		else{
			LOG("  WARINING in Refine: could not refine volume.\n");
		}
	}
	
//	done - clean up
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(false);

//	notify derivates that refinement ends
	refinement_step_ends();
}
	
}//	end of namespace
