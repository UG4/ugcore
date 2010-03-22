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
	m_copyRange = 2;
}

MultiGridRefiner::MultiGridRefiner(MultiGrid& mg)
{
	m_pMG = NULL;
	assign_grid(mg);
	m_copyRange = 2;
}

MultiGridRefiner::~MultiGridRefiner()
{
	if(m_pMG)
		m_pMG->unregister_observer(this);
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

	MultiGrid& mg = *m_pMG;

//	selection-marks
	m_selMarks.assign_grid(mg);
	m_selMarks.enable_autoselection(false);
	m_selMarks.enable_selection_inheritance(false);

//	status- and refinement-marks
	mg.attach_to_vertices_dv(m_aInt, 0);
	mg.attach_to_edges_dv(m_aInt, 0);
	mg.attach_to_faces_dv(m_aInt, 0);
	mg.attach_to_volumes_dv(m_aInt, 0);

	m_aaIntVRT.access(mg, m_aInt);
	m_aaIntEDGE.access(mg, m_aInt);
	m_aaIntFACE.access(mg, m_aInt);
	m_aaIntVOL.access(mg, m_aInt);
}

void MultiGridRefiner::unregistered_from_grid(Grid* grid)
{
	if(m_pMG)
	{
		m_pMG->detach_from_vertices(m_aInt);
		m_pMG->detach_from_edges(m_aInt);
		m_pMG->detach_from_faces(m_aInt);
		m_pMG->detach_from_volumes(m_aInt);

		m_pMG->unregister_observer(&m_selMarks);
		m_pMG = NULL;
	}
}

void MultiGridRefiner::clear_marks()
{
	m_selMarks.clear();
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


//	some buffers
	vector<VertexBase*> vVrts;
	vector<VertexBase*> vEdgeVrts;
	vector<EdgeBase*>	vEdges;
	vector<Face*>		vFaces;

//	some repeatedly used objects
	EdgeDescriptor ed;
	FaceDescriptor fd;
	VolumeDescriptor vd;

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
			//	change z-coord to visualise the hierarchy
				aaPos[nVrt].z += 0.01;
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
		int rule = get_rule(e);
		switch(rule)
		{
		case RM_COPY:
			{
			//	clone the edge.
				ed.set_vertices(mg.get_child_vertex(e->vertex(0)),
								mg.get_child_vertex(e->vertex(1)));
				EdgeBase* newEdge = *mg.create_by_cloning(e, ed, e);
				set_status(newEdge, SM_COPY);
			}break;

		default:
			{
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
					aaPos[nVrt].z += 0.01;
				}

			//	split the edge
				e->refine(vEdges, nVrt, substituteVrts);
				assert((vEdges.size() == 2) && "Edge refine produced wrong number of edges.");
				mg.register_element(vEdges[0], e);
				mg.register_element(vEdges[1], e);
				set_status(vEdges[0], SM_REGULAR);
				set_status(vEdges[1], SM_REGULAR);
			}break;
		}
	}

//LOG("creating new faces\n");
//	create new vertices and faces from marked faces
	for(FaceIterator iter = m_selMarks.begin<Face>();
		iter != m_selMarks.end<Face>(); ++iter)
	{
		Face* f = *iter;

		int rule = get_rule(f);
		switch(rule)
		{
		case RM_COPY:
			{
			//	clone the face.
				if(fd.num_vertices() != f->num_vertices())
					fd.set_num_vertices(f->num_vertices());

				for(size_t i = 0; i < fd.num_vertices(); ++i)
					fd.set_vertex(i, mg.get_child_vertex(f->vertex(i)));

				Face* nf = *mg.create_by_cloning(f, fd, f);
				set_status(nf, SM_COPY);
			}break;

		default:
			{
			//	collect child-vertices
				vVrts.clear();
				for(uint j = 0; j < f->num_vertices(); ++j){
					vVrts.push_back(mg.get_child_vertex(f->vertex(j)));
				}

			//	collect the associated edges
				vEdgeVrts.clear();
				bool bIrregular = false;
				for(uint j = 0; j < f->num_edges(); ++j){
					VertexBase* vrt = mg.get_child_vertex(mg.get_edge(f, j));
					vEdgeVrts.push_back(vrt);
					if(!vrt)
						bIrregular = true;
				}

				if(bIrregular){
					assert((get_rule(f) != RM_REGULAR) && "Bad refinement-rule set during collect_objects_for_refine!");
		//TODO:	care about anisotropy
					set_rule(f, RM_IRREGULAR);
				}

				VertexBase* newVrt;
				if(f->refine(vFaces, &newVrt, &vEdgeVrts.front(), NULL, &vVrts.front())){
				//	if a new vertex was generated, we have to register it
					if(newVrt){
						mg.register_element(newVrt, f);
						if(aaPos.valid()){
							aaPos[newVrt] = CalculateCenter(f, aaPos);
						//	change z-coord to visualise the hierarchy
							aaPos[newVrt].z += 0.01;
						}
					}

					int oldRule = get_rule(f);

				//	register the new faces and assign status
					for(size_t j = 0; j < vFaces.size(); ++j){
						mg.register_element(vFaces[j], f);
						switch(oldRule)
						{
						case RM_REGULAR:	set_status(vFaces[j], SM_REGULAR); break;
						case RM_IRREGULAR:	set_status(vFaces[j], SM_IRREGULAR); break;
						default:
							assert((oldRule == RM_REGULAR) && "rule not yet handled.");//always fails.
							break;
						}
					}
				}
				else{
					LOG("  WARINING in Refine: could not refine face.\n");
				}
			}
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

//	Some buffers and repeatedly used objects
	vector<EdgeBase*> 	vEdges;
	vector<Face*> 		vFaces;
	vector<Volume*>		vVolumes;

//	regard the multi-grid as a grid
	Grid& grid =*static_cast<Grid*>(m_pMG);

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
		else if(get_status(v) == SM_COPY){
		//	copy elements must not be refined in the current version
			LOG("  WARNING: copy-elements must not be refined in the current version!\n");
			m_selMarks.deselect(v);
		}
		else if(get_status(v) == SM_IRREGULAR){
		//	irregular elements must not be refined in the current version
			LOG("  WARNING: irregular-elements must not be refined in the current version!\n");
			m_selMarks.deselect(v);
		}
		else{
		//	mark it for regular refinement
			set_rule(v, RM_REGULAR);
		//	collect associated faces
			CollectFaces(vFaces, grid, v);
			for(uint i = 0; i < vFaces.size(); ++i)
				m_selMarks.select(vFaces[i]);
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
		else if(get_status(f) == SM_COPY){
		//	copy elements must not be refined in the current version
			LOG("  WARNING: copy-elements must not be refined in the current version!\n");
			m_selMarks.deselect(f);
		}
		else if(get_status(f) == SM_IRREGULAR){
		//	irregular elements must not be refined in the current version
			LOG("  WARNING: irregular-elements must not be refined in the current version!\n");
			m_selMarks.deselect(f);
		}
		else{
		//	mark it for regular refinement
			set_rule(f, RM_REGULAR);
		//	collect associated edges
			CollectEdges(vEdges, grid, f);
			for(uint i = 0; i < vEdges.size(); ++i)
				m_selMarks.select(vEdges[i]);
		}
	}

//	select elements that are associated with edges
	for(EdgeBaseIterator iter = m_selMarks.begin<EdgeBase>();
		iter != m_selMarks.end<EdgeBase>();)
	{
		EdgeBase* e = *iter;
		++iter;	//iterator is incremented since it would be invalidated during deselect.
		if(m_pMG->has_children(e)){
		//	the element has already been refined.
		//	don't refine it again.
			m_selMarks.deselect(e);
		}
		else if(get_status(e) == SM_COPY){
		//	copy elements must not be refined in the current version
			LOG("  WARNING: copy-elements must not be refined in the current version!\n");
			m_selMarks.deselect(e);
		}
		else if(get_status(e) == SM_IRREGULAR){
		//	irregular elements must not be refined in the current version
			LOG("  WARNING: irregular-elements must not be refined in the current version!\n");
			m_selMarks.deselect(e);
		}
		else{
		//	mark it for regular refinement
			set_rule(e, RM_REGULAR);
		//	select associated vertices
			for(uint i = 0; i < e->num_vertices(); ++i)
				m_selMarks.select(e->vertex(i));
		
		//	collect associated elements of higher dimension that
		//	will be used for the closure.
			CollectVolumes(vVolumes, grid, e);
			for(size_t i = 0; i < vVolumes.size(); ++i){
				if(!m_selMarks.is_selected(vVolumes[i])){
					Volume* v = vVolumes[i];
					m_selMarks.select(v);
				
				//	since all neighbours could possibly be marked, we have to
				//	assign rule UNKNOWN. in most cases this volume will be refined
				//	with irregular refinement later on.
					set_rule(v, RM_UNKNOWN);

				//	select associated vertices
					for(size_t j = 0; j < v->num_vertices(); ++j)
						m_selMarks.select(v->vertex(j));
				}
			}

			CollectFaces(vFaces, grid, e);
			for(size_t i = 0; i < vFaces.size(); ++i){
				if(!m_selMarks.is_selected(vFaces[i])){
					Face* f = vFaces[i];
					m_selMarks.select(f);
				
				//	since all neighbours could possibly be marked, we have to
				//	assign rule UNKNOWN. in most cases this face will be refined
				//	with irregular refinement later on.
					set_rule(f, RM_UNKNOWN);

				//	select associated vertices
					for(size_t j = 0; j < f->num_vertices(); ++j)
						m_selMarks.select(f->vertex(j));
				}
			}
		}
	}

//	collect copy-elements
//	we'll collect unselected edges, faces and volumes that are in
//	a neighbourhood to the selection
//TODO:	This could be optimized!
	if(get_copy_range() > 0){
		vector<VertexBase*> vVrts;
		vVrts.reserve(size_t((float)m_selMarks.num<VertexBase>() * 1.5f));
		
		for(VertexBaseIterator iter = m_selMarks.begin<VertexBase>();
			iter != m_selMarks.end<VertexBase>(); ++iter)
		{
			vVrts.push_back(*iter);
		}

	//	we have to make sure that we'll only collect copy-elements in
	//	the correct neighbourhood.
	//	After we processed all elements between iFirst and iEnd, the
	//	first neighbourhood is done. We may then set iFirst to iEnd and
	//	iEnd to vVrts.size(), to process the next neighbourhood.
		size_t iFirst = 0;
		size_t iEnd = vVrts.size();
	//	iterate for each neighbourhood
		for(size_t iNbr = 0; iNbr < get_copy_range(); ++iNbr)
		{
		//	iterate over candidates
			for(size_t i = iFirst; i != iEnd; ++i){
				VertexBase* vrt = vVrts[i];
			//	check all associated elements.
			//	If we can find an unselected, it is a copy-element.
			//	we then have to add all its associated unselected vertices
			//	to vVrts.
			
			//	check edges
				{
					EdgeBaseIterator iterEnd = grid.associated_edges_end(vrt);
					for(EdgeBaseIterator iter = grid.associated_edges_begin(vrt);
						iter != iterEnd; ++iter)
					{
						if(!m_selMarks.is_selected(*iter)){
							EdgeBase* e = *iter;
						//	we found a copy-element
							m_selMarks.select(e);
							set_rule(e, RM_COPY);
						//	schedule associated unselected vertices for further checks
							for(size_t j = 0; j < 2; ++j){
								if(!m_selMarks.is_selected(e->vertex(j))){
									m_selMarks.select(e->vertex(j));
									vVrts.push_back(e->vertex(j));
								}
							}
						}
					}
				}

			//	check faces
				{
					FaceIterator iterEnd = grid.associated_faces_end(vrt);
					for(FaceIterator iter = grid.associated_faces_begin(vrt);
						iter != iterEnd; ++iter)
					{
						if(!m_selMarks.is_selected(*iter)){
							Face* f = *iter;
						//	we found a copy-element
							m_selMarks.select(f);
							set_rule(f, RM_COPY);
						//	schedule associated unselected vertices for further checks
							for(size_t j = 0; j < f->num_vertices(); ++j){
								if(!m_selMarks.is_selected(f->vertex(j))){
									m_selMarks.select(f->vertex(j));
									vVrts.push_back(f->vertex(j));
								}
							}
						}
					}
				}

			//	check volumes
				{
					VolumeIterator iterEnd = grid.associated_volumes_end(vrt);
					for(VolumeIterator iter = grid.associated_volumes_begin(vrt);
						iter != iterEnd; ++iter)
					{
						if(!m_selMarks.is_selected(*iter)){
							Volume* v = *iter;
						//	we found a copy-element
							m_selMarks.select(v);
							set_rule(v, RM_COPY);
						//	schedule associated unselected vertices for further checks
							for(size_t j = 0; j < v->num_vertices(); ++j){
								if(!m_selMarks.is_selected(v->vertex(j))){
									m_selMarks.select(v->vertex(j));
									vVrts.push_back(v->vertex(j));
								}
							}
						}
					}
				}
			}

		//	the neighbourhood is done. check the next neighbourhood.
			iFirst = iEnd;
			iEnd = vVrts.size();
		}
	}
}

}//	end of namespace
