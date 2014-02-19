// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 12.09.2012 (d,m,y)

#include "lib_grid/algorithms/selection_util.h"
#include "adaptive_regular_mg_refiner.h"

using namespace std;

namespace ug{

AdaptiveRegularRefiner_MultiGrid::
AdaptiveRegularRefiner_MultiGrid(IRefinementCallback* refCallback) :
	HangingNodeRefiner_MultiGrid(refCallback)
{
}

AdaptiveRegularRefiner_MultiGrid::
AdaptiveRegularRefiner_MultiGrid(MultiGrid& mg,
								 IRefinementCallback* refCallback) :
	HangingNodeRefiner_MultiGrid(refCallback)
{
	set_grid(&mg);
}

AdaptiveRegularRefiner_MultiGrid::
~AdaptiveRegularRefiner_MultiGrid()
{
}

void AdaptiveRegularRefiner_MultiGrid::
assign_grid(MultiGrid& mg)
{
	set_grid(&mg);
}

void AdaptiveRegularRefiner_MultiGrid::
set_grid(MultiGrid* mg)
{
	HangingNodeRefiner_MultiGrid::set_grid(mg);

	m_closureElems.enable_autoselection(false);
	m_closureElems.enable_selection_inheritance(false);
	m_closureElems.enable_strict_inheritance(false);

	m_closureElems.assign_grid(mg);
}


void AdaptiveRegularRefiner_MultiGrid::
remove_closure_elements()
{
//	remove all closure elements from the current grid

	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

	EraseSelectedObjects(m_closureElems);
}


void AdaptiveRegularRefiner_MultiGrid::
create_closure_elements()
{
	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

	MultiGrid& mg = *multi_grid();

//	temporarily create refinement callback, if none is available.
//todo:	A better solution should be used for this!
	bool localRefCallbackSet = false;
	if(!m_refCallback){
		if(mg.has_vertex_attachment(aPosition)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition>(mg, aPosition);
		}
		else if(mg.has_vertex_attachment(aPosition2)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition2>(mg, aPosition2);
		}
		else if(mg.has_vertex_attachment(aPosition1)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition1>(mg, aPosition1);
		}
	}

	if(mg.num<Volume>() > 0)
		create_closure_elements_3d();
	else if(mg.num<Face>() > 0)
		create_closure_elements_2d();

//	clear the refinement-callback if necessary
	if(localRefCallbackSet){
		delete m_refCallback;
		m_refCallback = NULL;
	}
}


void AdaptiveRegularRefiner_MultiGrid::
create_closure_elements_2d()
{
//todo:	This method currently only works for one type of elements. No manifolds
//		are currently supported.

//	for each surface element (faces in 2d, volumes in 3d) adjacent to a constraining
//	element, we generate closure elements in the level above.

	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

//	remove all closure elements. This is currently required, since the selector
//	m_closureElems is currently also used to store non-closer elements, which are
//	to be refined.
	remove_closure_elements();

	MultiGrid& mg = *multi_grid();

//	iterate over all constraining edges and collect associated surface faces
//	Those then have to be refined to generate a closure.

	Grid::face_traits::secure_container	assElems;
	Grid::edge_traits::secure_container assEdges;
	Face::ConstVertexArray vrts;
	std::vector<VertexBase*> newVrtVrts;
	std::vector<VertexBase*> newEdgeVrts;
	std::vector<Face*>	newFaces;
	EdgeDescriptor ed;

//	we'll select all new elements on the fly
	m_closureElems.enable_autoselection(true);

	for(Grid::traits<ConstrainingEdge>::iterator i_edge = mg.begin<ConstrainingEdge>();
		i_edge != mg.end<ConstrainingEdge>(); ++i_edge)
	{
	//	check all associated elements of i_edge, whether one is a surface element.
	//	If so, create the closure.
		mg.associated_elements(assElems, *i_edge);

		for(size_t i_elem = 0; i_elem < assElems.size(); ++i_elem){
			Face* elem = assElems[i_elem];
			if(mg.has_children(elem))
				continue;

			newVrtVrts.clear();
			newEdgeVrts.clear();

		//	copy associated vertices and edges to the next level
			vrts = elem->vertices();
			size_t numVrts = elem->num_vertices();
			for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
				VertexBase* vrt = vrts[i_vrt];
				if(!mg.has_children(vrt)){
					newVrtVrts.push_back(*mg.create<RegularVertex>(vrt));
					if(m_refCallback)
						m_refCallback->new_vertex(newVrtVrts.back(), vrt);
				}
				else
					newVrtVrts.push_back(mg.get_child_vertex(vrt));
			}

			mg.associated_elements_sorted(assEdges, elem);
			for(size_t i = 0; i < assEdges.size(); ++i){
				EdgeBase* e = assEdges[i];
				if(!mg.has_children(e)){
					ed.set_vertices(mg.get_child_vertex(e->vertex(0)),
									mg.get_child_vertex(e->vertex(1)));
					mg.create<Edge>(ed, e);
					newEdgeVrts.push_back(NULL);
				}
				else
					newEdgeVrts.push_back(mg.get_child_vertex(e));
			}

		//	refine the element
			VertexBase* newFaceVrt;
			if(elem->refine(newFaces, &newFaceVrt, &newEdgeVrts.front(),
							NULL, &newVrtVrts.front()))
			{
				if(newFaceVrt){
					mg.register_element(newFaceVrt, elem);
					if(m_refCallback)
						m_refCallback->new_vertex(newFaceVrt, elem);
				}
				for(size_t i = 0; i < newFaces.size(); ++i)
					mg.register_element(newFaces[i], elem);
			}
		}
	}

//	stop selecting new elements
	m_closureElems.enable_autoselection(false);
}

void AdaptiveRegularRefiner_MultiGrid::
create_closure_elements_3d()
{
//todo:	This method currently only works for one type of elements. No manifolds
//		are currently supported.

//	for each surface element (faces in 2d, volumes in 3d) adjacent to a constraining
//	element, we generate closure elements in the level above.

	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

//	remove all closure elements. This is currently required, since the selector
//	m_closureElems is currently also used to store non-closer elements, which are
//	to be refined.
	remove_closure_elements();

	MultiGrid& mg = *multi_grid();

//	iterate over all constraining edges and collect associated surface faces / volumes.
//	Those then have to be refined to generate a closure.

	Grid::volume_traits::secure_container	assElems;
	Grid::edge_traits::secure_container assEdges;
	Grid::face_traits::secure_container assFaces;
	Volume::ConstVertexArray vrts;
	std::vector<VertexBase*> newVolVrtVrts;
	std::vector<VertexBase*> newVolEdgeVrts;
	std::vector<VertexBase*> newVolFaceVrts;
	std::vector<Volume*> newVols;
	EdgeDescriptor ed;
	FaceDescriptor fd;
	
//	when refining the associated faces, we need some structs, too
	Grid::edge_traits::secure_container assFaceEdges;
	std::vector<VertexBase*> newFaceVrtVrts;
	std::vector<VertexBase*> newFaceEdgeVrts;
	std::vector<Face*> newFaces;

//	we'll select all new elements on the fly
	m_closureElems.enable_autoselection(true);

	for(Grid::traits<ConstrainingEdge>::iterator i_edge = mg.begin<ConstrainingEdge>();
		i_edge != mg.end<ConstrainingEdge>(); ++i_edge)
	{
	//	check all associated elements of i_edge, whether one is a surface element.
	//	If so, select it into m_closureElems. This is only temporary, since it
	//	isn't a real closure element.
		mg.associated_elements(assElems, *i_edge);

		for(size_t i_elem = 0; i_elem < assElems.size(); ++i_elem){
			Volume* elem = assElems[i_elem];
			if(mg.has_children(elem))
				continue;

			newVolVrtVrts.clear();
			newVolEdgeVrts.clear();
			newVolFaceVrts.clear();

		//	copy associated vertices and edges to the next level
			vrts = elem->vertices();
			size_t numVrts = elem->num_vertices();
			for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
				VertexBase* vrt = vrts[i_vrt];
				if(!mg.has_children(vrt)){
					newVolVrtVrts.push_back(*mg.create<RegularVertex>(vrt));
					if(m_refCallback)
						m_refCallback->new_vertex(newVolVrtVrts.back(), vrt);
				}
				else
					newVolVrtVrts.push_back(mg.get_child_vertex(vrt));
			}

			mg.associated_elements_sorted(assEdges, elem);
			for(size_t i = 0; i < assEdges.size(); ++i){
				EdgeBase* e = assEdges[i];
				if(!mg.has_children(e)){
					ed.set_vertices(mg.get_child_vertex(e->vertex(0)),
									mg.get_child_vertex(e->vertex(1)));
					mg.create<Edge>(ed, e);
					newVolEdgeVrts.push_back(NULL);
				}
				else
					newVolEdgeVrts.push_back(mg.get_child_vertex(e));
			}

		//	we have to either refine or copy associated faces
			mg.associated_elements_sorted(assFaces, elem);
			for(size_t i_face = 0; i_face < assFaces.size(); ++i_face){
				
				Face* f = assFaces[i_face];
				if(mg.has_children(f)){
					newVolFaceVrts.push_back(mg.get_child_vertex(f));
					continue;
				}
				
			//	collect child vertices of associated edges
				mg.associated_elements_sorted(assFaceEdges, f);
				newFaceEdgeVrts.clear();
				
				bool faceRefinement = false;
				for(size_t i = 0; i < assFaceEdges.size(); ++i){
					VertexBase* child = mg.get_child_vertex(assFaceEdges[i]);
					newFaceEdgeVrts.push_back(child);
					faceRefinement |= (child != NULL);
				}				
				
				if(faceRefinement){
					newFaceVrtVrts.clear();
					for(size_t i = 0; i < f->num_vertices(); ++i)
						newFaceVrtVrts.push_back(mg.get_child_vertex(f->vertex(i)));
					
					VertexBase* newFaceVrt = NULL;
					if(f->refine(newFaces, &newFaceVrt, &newFaceEdgeVrts.front(),
								 NULL, &newFaceVrtVrts.front()))
					{
						if(newFaceVrt){
							mg.register_element(newFaceVrt, f);
							if(m_refCallback)
								m_refCallback->new_vertex(newFaceVrt, f);
						}
						for(size_t i = 0; i < newFaces.size(); ++i)
							mg.register_element(newFaces[i], f);
					}
					newVolFaceVrts.push_back(newFaceVrt);
				}	
				else{
					Face::ConstVertexArray fvrts = f->vertices();
					if(f->num_vertices() == 3)
						mg.create<Triangle>(TriangleDescriptor(
												mg.get_child_vertex(fvrts[0]),
												mg.get_child_vertex(fvrts[1]),
												mg.get_child_vertex(fvrts[2])),
											f);
					else if(f->num_vertices() == 4)
						mg.create<Quadrilateral>(QuadrilateralDescriptor(
													mg.get_child_vertex(fvrts[0]),
													mg.get_child_vertex(fvrts[1]),
													mg.get_child_vertex(fvrts[2]),
													mg.get_child_vertex(fvrts[3])),
												f);
					newVolFaceVrts.push_back(NULL);
				}
			}

		//	refine the element
		//	if we're performing tetrahedral refinement, we have to collect
		//	the corner coordinates, so that the refinement algorithm may choose
		//	the best interior diagonal.
			vector3 corners[4];
			vector3* pCorners = NULL;
			if((elem->num_vertices() == 4) && m_refCallback){
				for(size_t i = 0; i < 4; ++i){
					m_refCallback->current_pos(&corners[i].x(), vrts[i], 3);
				}
				pCorners = corners;
			}

			VertexBase* newVolVrt;
			if(elem->refine(newVols, &newVolVrt, &newVolEdgeVrts.front(),
							&newVolFaceVrts.front(), NULL, RegularVertex(),
							&newVolVrtVrts.front(), pCorners))
			{
				if(newVolVrt){
					mg.register_element(newVolVrt, elem);
					if(m_refCallback)
						m_refCallback->new_vertex(newVolVrt, elem);
				}

				for(size_t i = 0; i < newVols.size(); ++i)
					mg.register_element(newVols[i], elem);
			}
		}
	}

//	stop selecting new elements
	m_closureElems.enable_autoselection(false);
}

template <class TElem>
void AdaptiveRegularRefiner_MultiGrid::
get_parents_of_marked_closure_elements(std::vector<GridObject*>& parents,
									   Selector::status_t mark)
{
	UG_ASSERT(multi_grid(), "A multi grid has to be assigned to the refiner.");
	MultiGrid& mg = *multi_grid();

	typedef typename BaseClass::selector_t::template traits<TElem>::iterator	TIter;
	for(TIter iter = m_selMarkedElements.begin<TElem>();
		iter != m_selMarkedElements.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		if(!m_closureElems.is_selected(e))
			continue;

		if(get_mark(e) & mark){
			GridObject* parent = mg.get_parent(e);
			if(parent && !m_closureElems.is_selected(parent))
				parents.push_back(parent);
		}
	}
}

void AdaptiveRegularRefiner_MultiGrid::
perform_refinement()
{
//	todo: copy refinement marks from closure elements to their parents
	vector<GridObject*> parents;
	Selector::status_t refMark = RM_REFINE | RM_ANISOTROPIC;
	get_parents_of_marked_closure_elements<VertexBase>(parents, refMark);
	get_parents_of_marked_closure_elements<EdgeBase>(parents, refMark);
	get_parents_of_marked_closure_elements<Face>(parents, refMark);
	get_parents_of_marked_closure_elements<Volume>(parents, refMark);

	remove_closure_elements();

//	mark parents of formerly marked closure elements for refinement
	mark(parents.begin(), parents.end(), RM_REFINE);

	HangingNodeRefiner_MultiGrid::perform_refinement();
	create_closure_elements();
}

bool AdaptiveRegularRefiner_MultiGrid::
perform_coarsening()
{
//	todo: copy coarsen marks from closure elements to their parents
	remove_closure_elements();
	bool bSuccess = HangingNodeRefiner_MultiGrid::perform_coarsening();
	create_closure_elements();
	return bSuccess;
}

}// end of namespace
