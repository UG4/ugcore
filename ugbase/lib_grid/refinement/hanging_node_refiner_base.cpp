/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <vector>
#include "hanging_node_refiner_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/file_io/file_io.h"
#include "ref_mark_adjusters/std_hnode_adjuster.h"
#include "lib_grid/tools/selector_multi_grid.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
#endif

//define PROFILE_HANGING_NODE_REFINER if you want to profile
//the refinement code.
#define PROFILE_HANGING_NODE_REFINER
#ifdef PROFILE_HANGING_NODE_REFINER
	#define HNODE_PROFILE_FUNC()	PROFILE_FUNC()
	#define HNODE_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define HNODE_PROFILE_END()	PROFILE_END()
#else
	#define HNODE_PROFILE_FUNC()
	#define HNODE_PROFILE_BEGIN(name)
	#define HNODE_PROFILE_END()
#endif


using namespace std;

namespace ug{
template <class TSelector>
HangingNodeRefinerBase<TSelector>::
HangingNodeRefinerBase(SPRefinementProjector projector) :
	IRefiner(projector),
	m_pGrid(NULL),
	m_nodeDependencyOrder1(true),
	m_adjustingRefMarks(false)
	//,m_automarkHigherDimensionalObjects(false)
{
	add_ref_mark_adjuster(StdHNodeAdjuster::create());
}

template <class TSelector>
HangingNodeRefinerBase<TSelector>::
HangingNodeRefinerBase(const HangingNodeRefinerBase&)
{
	throw(UGError("no copy construction allowed."));
}

template <class TSelector>
HangingNodeRefinerBase<TSelector>::
~HangingNodeRefinerBase()
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
set_grid(typename TSelector::grid_type* grid)
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
		m_selMarkedElements.assign_grid(NULL);
		m_pGrid = NULL;
	}

	if(grid){
		m_pGrid = grid;
		grid->register_observer(this, OT_GRID_OBSERVER);
		m_selMarkedElements.assign_grid(*grid);
		m_selMarkedElements.enable_autoselection(false);
		m_selMarkedElements.enable_selection_inheritance(false);
		set_message_hub(grid->message_hub());
	}
}

template <class TSelector> void
HangingNodeRefinerBase<TSelector>::grid_to_be_destroyed(Grid* grid)
{
	m_pGrid = NULL;
}

template <class TSelector> void
HangingNodeRefinerBase<TSelector>::clear_marks()
{
	m_selMarkedElements.clear();
	m_newlyMarkedRefVrts.clear();
	m_newlyMarkedRefEdges.clear();
	m_newlyMarkedRefFaces.clear();
	m_newlyMarkedRefVols.clear();
}

template <class TSelector>
bool HangingNodeRefinerBase<TSelector>::
mark(Vertex* v, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");
	if(get_mark(v) != refMark){
	//	if ref-mark-adjustment is taking place, we'll also allow for the selection
	//	of non-surface vertices
		if(refinement_is_allowed(v) || m_adjustingRefMarks){
			m_selMarkedElements.select(v, refMark);
			if(m_adjustingRefMarks && (refMark & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_DUMMY)))
				m_newlyMarkedRefVrts.push_back(v);
			return true;
		}
	}
	return false;
}

template <class TSelector>
bool HangingNodeRefinerBase<TSelector>::mark(Edge* e, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");
	if(get_mark(e) != refMark){
		if(refinement_is_allowed(e)){
			m_selMarkedElements.select(e, refMark);
			if(m_adjustingRefMarks && (refMark & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_DUMMY)))
				m_newlyMarkedRefEdges.push_back(e);
			return true;
		}
	}
	return false;
}

template <class TSelector>
bool HangingNodeRefinerBase<TSelector>::mark(Face* f, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");

	if(get_mark(f) != refMark){
		if(refinement_is_allowed(f)){
			m_selMarkedElements.select(f, refMark);
			if(m_adjustingRefMarks && (refMark & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_DUMMY)))
				m_newlyMarkedRefFaces.push_back(f);
			return true;
		}
	}
	return false;
}

template <class TSelector>
bool HangingNodeRefinerBase<TSelector>::mark(Volume* v, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");
	if(get_mark(v) != refMark){
		if(refinement_is_allowed(v)){
			m_selMarkedElements.select(v, refMark);
			if(m_adjustingRefMarks && (refMark & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_DUMMY)))
				m_newlyMarkedRefVols.push_back(v);
			return true;
		}
	}
	return false;
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
mark_neighborhood(size_t numIterations, RefinementMark refMark, bool sideNbrsOnly)
{
	if(!m_pGrid)
		return;

	typedef typename TSelector::template traits<Vertex>::iterator	SelVrtIter;
	typedef typename TSelector::template traits<Edge>::iterator		SelEdgeIter;
	typedef typename TSelector::template traits<Face>::iterator			SelFaceIter;
	typedef typename TSelector::template traits<Volume>::iterator		SelVolIter;

	Grid::edge_traits::secure_container		edges;
	Grid::face_traits::secure_container		faces;
	Grid::volume_traits::secure_container	vols;

	EdgeDescriptor edesc;

	TSelector& sel = m_selMarkedElements;
	Grid& g = *m_pGrid;

	#ifdef UG_PARALLEL
		bool hasVolumes = pcl::OneProcTrue(g.num_volumes() > 0);
		bool hasFaces;
		if(hasVolumes)
			hasFaces = true;
		else
			hasFaces = pcl::OneProcTrue(g.num_faces() > 0);
		bool hasEdges;
		if(hasFaces)
			hasEdges = true;
		else
			hasEdges = pcl::OneProcTrue(g.num_edges() > 0);
	#else
		bool hasEdges	= g.num_edges() > 0;
		bool hasFaces	= g.num_faces() > 0;
		bool hasVolumes = g.num_volumes() > 0;
	#endif

//todo:	Speed could be improved by only considering newly marked elements after each iteration.
	for(size_t iteration = 0; iteration < numIterations; ++iteration){
		if(sideNbrsOnly){
		//	first we'll select all unselected sides which are connected to marked elements
			for(SelVolIter iter = sel.template begin<Volume>();
				iter != sel.template end<Volume>(); ++iter)
			{
				Volume* vol = *iter;
				RefinementMark s = get_mark(vol);
				g.associated_elements(faces, vol);
				for_each_in_vec(Face* f, faces){
					if(!sel.is_selected(f) && this->refinement_is_allowed(f))
						this->mark(f, s);
				}end_for;
			}

			for(SelFaceIter iter = sel.template begin<Face>();
				iter != sel.template end<Face>(); ++iter)
			{
				Face* f = *iter;
				RefinementMark s = get_mark(f);
				g.associated_elements(edges, f);
				for_each_in_vec(Edge* e, edges){
					if(!sel.is_selected(e) && this->refinement_is_allowed(e))
						this->mark(e, s);
				}end_for;
			}

			for(SelEdgeIter iter = sel.template begin<Edge>();
				iter != sel.template end<Edge>(); ++iter)
			{
				Edge* e = *iter;
				RefinementMark s = get_mark(e);
				for(size_t i = 0; i < 2; ++i){
					Vertex* vrt = e->vertex(i);
					if(!sel.is_selected(vrt) && this->refinement_is_allowed(vrt))
						mark(vrt, s);
				}
			}
		}
		else{
		//	first we'll select all vertices which are connected to marked elements
			for(SelEdgeIter iter = sel.template begin<Edge>();
				iter != sel.template end<Edge>(); ++iter)
			{
				Edge* e = *iter;
				sel.select(e->vertex(0), sel.get_selection_status(e));
				sel.select(e->vertex(1), sel.get_selection_status(e));
			}

			for(SelFaceIter iter = sel.template begin<Face>();
				iter != sel.template end<Face>(); ++iter)
			{
				Face* f = *iter;
				ISelector::status_t s = sel.get_selection_status(f);
				Face::ConstVertexArray faceVrts = f->vertices();
				for(size_t i = 0; i < f->num_vertices(); ++i)
					sel.select(faceVrts[i], s);
			}

			for(SelVolIter iter = sel.template begin<Volume>();
				iter != sel.template end<Volume>(); ++iter)
			{
				Volume* vol = *iter;
				ISelector::status_t s = sel.get_selection_status(vol);
				Volume::ConstVertexArray volVrts = vol->vertices();
				for(size_t i = 0; i < vol->num_vertices(); ++i)
					sel.select(volVrts[i], s);
			}
		}

	//	if we're in a parallel environment, we have to broadcast the current selection
		sel.broadcast_selection_states();

		if(sideNbrsOnly){
		//	mark those elements which have a marked side
			if(hasVolumes){
				vector<Edge*> markedDummyEdges;
				lg_for_each(Volume, vol, g){
					if(this->refinement_is_allowed(vol)){
						RefinementMark s = get_mark(vol);
						if(s < refMark){
							g.associated_elements(faces, vol);
							RefinementMark rm = RM_NONE;
							for_each_in_vec(Face* f, faces){
								RefinementMark sSide = get_mark(f);
								if(sSide){
									rm = refMark;
									if(rm == RM_NONE){
										switch(sSide){
											case RM_REFINE:	rm = RM_REFINE; break;
											case RM_ANISOTROPIC: rm = RM_ANISOTROPIC; break;
											case RM_CLOSURE: rm = RM_CLOSURE; break;
											case RM_COARSEN: rm = RM_COARSEN; break;
											default: break;
										}
									}
									mark(vol, rm);
									break;
								}
							}end_for;

							if(rm == RM_ANISOTROPIC){
							//	we'll also copy edge-marks to guarantee an anisotropic refinement
								bool copying = true;
								while(copying){
									copying = false;
									for_each_in_vec(Face* f, faces){
										g.associated_elements(edges, f);
										for_each_in_vec(Edge* e, edges){
											RefinementMark erm = get_mark(e);
											if(erm == RM_REFINE){
												if(f->get_opposing_side(e, edesc)){
													Edge* oe = g.get_edge(edesc);
													if(oe && get_mark(oe) < erm){
														copying = true;
														mark(oe, RM_DUMMY);
														markedDummyEdges.push_back(oe);
													}
												}
											}	
										}end_for;
									}end_for;
								}
							}
						}
					}
				}lg_end_for;

				for_each_in_vec(Edge* e, markedDummyEdges){
					mark(e, RM_REFINE);	
				}end_for;
			}
			else if(hasFaces){
				lg_for_each(Face, f, g){
					if(this->refinement_is_allowed(f)){
						RefinementMark s = get_mark(f);
						if(s < refMark){
							g.associated_elements(edges, f);
							RefinementMark rm = RM_NONE;
							for_each_in_vec(Edge* e, edges){
								RefinementMark sSide = get_mark(e);
								if(sSide){
									rm = refMark;
									if(rm == RM_NONE){
										switch(sSide){
											case RM_REFINE:	rm = RM_REFINE; break;
											case RM_ANISOTROPIC: rm = RM_ANISOTROPIC; break;
											case RM_CLOSURE: rm = RM_CLOSURE; break;
											case RM_COARSEN: rm = RM_COARSEN; break;
											default: break;
										}
									}
									mark(f, rm);
									break;
								}
							}end_for;

						//	todo: if rm == RM_ANISOTROPIC: mark opposite edges (cf hasVolumes)
						}
					}
				}lg_end_for;
			}
			else if(hasEdges){
				lg_for_each(Edge, e, g){
					if(this->refinement_is_allowed(e)){
						RefinementMark s = get_mark(e);
						if(s < refMark){
							Edge::ConstVertexArray vrts = e->vertices();
							for(size_t i = 0; i < 2; ++i){
								ISelector::status_t sSide = sel.get_selection_status(vrts[i]);
								if(sSide){
									RefinementMark rm = refMark;
									if(rm == RM_NONE){
										switch(sSide){
											case RM_REFINE:	rm = RM_REFINE; break;
											case RM_ANISOTROPIC: rm = RM_ANISOTROPIC; break;
											case RM_CLOSURE: rm = RM_CLOSURE; break;
											case RM_COARSEN: rm = RM_COARSEN; break;
											default: break;
										}
									}
									mark(e, rm);
								}
							}
						}
					}
				}lg_end_for;
			}
		}
		else{
		//	mark all associated elements of marked vertices
			for(SelVrtIter iter = sel.template begin<Vertex>();
				iter != sel.template end<Vertex>(); ++iter)
			{
				ISelector::status_t s = sel.get_selection_status(*iter);
				RefinementMark rm = refMark;
				if(rm == RM_NONE){
					switch(s){
						case RM_REFINE:	rm = RM_REFINE; break;
						case RM_ANISOTROPIC: rm = RM_ANISOTROPIC; break;
						case RM_CLOSURE: rm = RM_CLOSURE; break;
						case RM_COARSEN: rm = RM_COARSEN; break;
						default: break;
					}
				}

				g.associated_elements(edges, *iter);
				for(size_t i = 0; i < edges.size(); ++i)
					if(!this->is_marked(edges[i]))
						this->mark(edges[i], rm);

				g.associated_elements(faces, *iter);
				for(size_t i = 0; i < faces.size(); ++i)
					if(!this->is_marked(faces[i]))
						this->mark(faces[i], rm);

				g.associated_elements(vols, *iter);
				for(size_t i = 0; i < vols.size(); ++i)
					if(!this->is_marked(vols[i]))
						this->mark(vols[i], rm);
			}

		//	since we selected vertices which possibly may not be refined, we have to
		//	deselect those now.
			for(SelVrtIter iter = sel.template begin<Vertex>();
				iter != sel.template end<Vertex>();)
			{
				Vertex* vrt = *iter;
				++iter;

				if(!refinement_is_allowed(vrt))
					sel.deselect(vrt);
			}
		}
		
	}

}

template <class TSelector>
RefinementMark HangingNodeRefinerBase<TSelector>::
get_mark(Vertex* v)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(v)
							& (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_COARSEN | RM_DUMMY));
}

template <class TSelector>
RefinementMark HangingNodeRefinerBase<TSelector>::
get_mark(Edge* e)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(e)
							& (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_COARSEN | RM_DUMMY));
}

template <class TSelector>
RefinementMark HangingNodeRefinerBase<TSelector>::
get_mark(Face* f)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(f)
							& (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_COARSEN | RM_DUMMY));
}

template <class TSelector>
RefinementMark HangingNodeRefinerBase<TSelector>::
get_mark(Volume* v)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(v)
							& (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE | RM_COARSEN | RM_DUMMY));
}

template <class TSelector>
bool HangingNodeRefinerBase<TSelector>::
save_marks_to_file(const char* filename)
{
	UG_DLOG(LIB_GRID, 1, "  saving marks to file...\n");
	if(!m_pGrid){
		UG_THROW("ERROR in HangingNodeRefinerBase::save_marks_to_file: No grid assigned!");
	}

	Grid& g = *m_pGrid;
	SubsetHandler sh(g);

	AssignGridToSubset(g, sh, 4);

	selector_t& sel = get_refmark_selector();

	for(VertexIterator iter = g.vertices_begin(); iter != g.vertices_end(); ++iter){
		typename selector_t::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_NONE: break;
			case RM_REFINE: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_CLOSURE: sh.assign_subset(*iter, 2); break;
			case RM_COARSEN: sh.assign_subset(*iter, 3); break;
			default: sh.assign_subset(*iter, 4); break;
		}
	}

	for(EdgeIterator iter = g.edges_begin(); iter != g.edges_end(); ++iter){
		typename selector_t::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_NONE: break;
			case RM_REFINE: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_CLOSURE: sh.assign_subset(*iter, 2); break;
			case RM_COARSEN: sh.assign_subset(*iter, 3); break;
			default: sh.assign_subset(*iter, 4); break;
		}
	}

	for(FaceIterator iter = g.faces_begin(); iter != g.faces_end(); ++iter){
		typename selector_t::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_NONE: break;
			case RM_REFINE: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_CLOSURE: sh.assign_subset(*iter, 2); break;
			case RM_COARSEN: sh.assign_subset(*iter, 3); break;
			default: sh.assign_subset(*iter, 4); break;
		}
	}

	for(VolumeIterator iter = g.volumes_begin(); iter != g.volumes_end(); ++iter){
		typename  selector_t::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_NONE: break;
			case RM_REFINE: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_CLOSURE: sh.assign_subset(*iter, 2); break;
			case RM_COARSEN: sh.assign_subset(*iter, 3); break;
			default: sh.assign_subset(*iter, 4); break;
		}
	}

	sh.subset_info(0).name = "refine regular";
	sh.subset_info(1).name = "refine anisotropic";
	sh.subset_info(2).name = "refine closure";
	sh.subset_info(3).name = "coarsen";
	sh.subset_info(4).name = "combi-mark";
	sh.subset_info(5).name = "no marks";

	AssignSubsetColors(sh);

	if(MultiGrid* pmg = dynamic_cast<MultiGrid*>(&g))
		return SaveGridHierarchyTransformed(*pmg, sh, filename, 0.1);
	else
		return SaveGridToFile(g, sh, filename);
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
enable_node_dependency_order_1(bool bEnable)
{
	m_nodeDependencyOrder1 = bEnable;
	for(size_t i = 0; i < m_refMarkAdjusters.size(); ++i){
		m_refMarkAdjusters[i]->enable_node_dependency_order_1(bEnable);
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::perform_refinement()
{
	HNODE_PROFILE_BEGIN(perform_hnode_refinement);
	UG_DLOG(LIB_GRID, 1, "performing hanging-node-refine:\n");

	if(!m_pGrid)
		throw(UGError("ERROR in HangingNodeRefinerBase::refine(...): No grid assigned."));

	if(m_selMarkedElements.grid() != m_pGrid)
		throw(UGError("selector not initialized properly. Use HangingNodeRefinerBase::set_grid."));

	Grid& grid = *m_pGrid;

//	check grid options.
	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING in HangingNodeRefiner::refine(): grid option GRIDOPT_AUTOGENERATE_SIDES auto-enabled." << endl);
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	containers used for temporary results
	vector<Edge*> 	vEdges;
	vector<Face*>	 	vFaces;
	vector<Volume*>		vVols;

	HNODE_PROFILE_BEGIN(href_AdjustingMarks);
//	fills the queues with the elements that have to be refined.
	collect_objects_for_refine();

//	assigns hnode marks
	assign_hnode_marks();
	HNODE_PROFILE_END();

//	{
//		UG_LOG("DEBUG SAVE...\n");
//		static int refselCount = 0;
//		stringstream ss;
//		ss << "refselbefore_" << refselCount << ".ugx";
//		save_marks_to_file(ss.str().c_str());
//		++refselCount;
//	}

//	if a debug file was specified, we'll now save the marks to that file
	if(!m_adjustedMarksDebugFilename.empty())
		save_marks_to_file(m_adjustedMarksDebugFilename.c_str());

	m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_REFINEMENT_BEGINS,
										m_selMarkedElements.get_grid_objects()));

//	call pre_refine to allow derived classes to perform some actions
	HNODE_PROFILE_BEGIN(href_PreRefine);
	pre_refine();
	HNODE_PROFILE_END();

////////////////////////////////
//	ConstrainedVertices
	UG_DLOG(LIB_GRID, 1, "  constrained vertices.\n");
	HNODE_PROFILE_BEGIN(href_ConstrainedVertices);
	for(typename selector_t::template traits<ConstrainedVertex>::iterator
			iter = m_selMarkedElements.template begin<ConstrainedVertex>();
		iter != m_selMarkedElements.template end<ConstrainedVertex>();)
	{
		ConstrainedVertex* cdv = *iter;
		++iter;
		process_constrained_vertex(cdv);
	}
	HNODE_PROFILE_END();

////////////////////////////////
//	ConstrainedEdges
	UG_DLOG(LIB_GRID, 1, "  constrained edges.\n");
	HNODE_PROFILE_BEGIN(href_ConstrainedEdges);
	for(typename selector_t::template traits<ConstrainedEdge>::iterator iter =
			m_selMarkedElements.template begin<ConstrainedEdge>();
		iter != m_selMarkedElements.template end<ConstrainedEdge>();)
	{
		ConstrainedEdge* cde = *iter;
		++iter;
		if(marked_to_normal(cde))
			remove_hmark(cde, RM_REFINE);
		process_constrained_edge(cde);
	}
	HNODE_PROFILE_END();

////////////////////////////////
//	ConstrainingEdges
//	iterate through scheduled cg-edges
	UG_DLOG(LIB_GRID, 1, "  constraining edges.\n");
	HNODE_PROFILE_BEGIN(href_ConstrainingEdges);
	{
		typename selector_t::template traits<ConstrainingEdge>::iterator iter =
			m_selMarkedElements.template begin<ConstrainingEdge>();

		while(iter != m_selMarkedElements.template end<ConstrainingEdge>()){
			ConstrainingEdge* cge = *iter;
			++iter;
			if(marked_to_normal(cge))
				remove_hmark(cge, RM_REFINE);
		//	if a constraining edge was previously created through replacement of
		//	a constrained edge, we won't process it here
			if(!marked_to_constraining(cge))
				process_constraining_edge(cge);
		}
	}
	HNODE_PROFILE_END();

////////////////////////////////
//	Normal Edges
//	iterate through scheduled edges
	UG_DLOG(LIB_GRID, 1, "  normal edges.\n");
	HNODE_PROFILE_BEGIN(href_NormalEdges);
	{
		typename selector_t::template traits<RegularEdge>::iterator
			iter = m_selMarkedElements.template begin<RegularEdge>();
		while(iter != m_selMarkedElements.template end<RegularEdge>())
		{
			RegularEdge* e = *iter;
			++iter;

		//	a normal edge may have previously been created by replacing a
		//	constrained or constraining edge. Those edges won't be considered here
			if((!refinement_is_allowed(e)) || marked_to_normal(e)){
				continue;
			}

			if(marked_to_constraining(e)){
				refine_edge_with_hanging_vertex(e);
			}
			else if(marked_refine(e) || marked_copy(e)){
				refine_edge_with_normal_vertex(e);
			}
		}
	}
	HNODE_PROFILE_END();

////////////////////////////////
//	Faces
	HNODE_PROFILE_BEGIN(href_ConstrainedFaces);
	UG_DLOG(LIB_GRID, 1, "  constrained triangles.\n");
	for(typename selector_t::template traits<ConstrainedTriangle>::iterator
			iter = m_selMarkedElements.template begin<ConstrainedTriangle>();
		iter != m_selMarkedElements.template end<ConstrainedTriangle>();)
	{
		ConstrainedTriangle* cdf = *iter;
		++iter;
		if(marked_to_normal(cdf))
			remove_hmark(cdf, RM_REFINE);
		process_constrained_face(cdf);
	}

	UG_DLOG(LIB_GRID, 1, "  constrained quadrilaterals.\n");
	for(typename selector_t::template traits<ConstrainedQuadrilateral>::iterator
			iter = m_selMarkedElements.template begin<ConstrainedQuadrilateral>();
		iter != m_selMarkedElements.template end<ConstrainedQuadrilateral>();)
	{
		ConstrainedQuadrilateral* cdf = *iter;
		++iter;
		if(marked_to_normal(cdf))
			remove_hmark(cdf, RM_REFINE);
		process_constrained_face(cdf);
	}
	HNODE_PROFILE_END();

	HNODE_PROFILE_BEGIN(href_ConstrainingFaces);
	UG_DLOG(LIB_GRID, 1, "  constraining triangles.\n");
//	constraining triangles
	{
		typename selector_t::template traits<ConstrainingTriangle>::iterator
			iter = m_selMarkedElements.template begin<ConstrainingTriangle>();
		while(iter != m_selMarkedElements.template end<ConstrainingTriangle>()){
			ConstrainingTriangle* cgf = *iter;
			++iter;
			if(marked_to_normal(cgf))
				remove_hmark(cgf, RM_REFINE);
		//	if a constraining face was previously created through replacement of
		//	a constrained face, we won't process it here
			if(!marked_to_constraining(cgf))
				process_constraining_face(cgf);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  constraining quadrilaterals.\n");
//	constraining quadrilaterals
	{
		typename selector_t::template traits<ConstrainingQuadrilateral>::iterator
			iter = m_selMarkedElements.template begin<ConstrainingQuadrilateral>();
		while(iter != m_selMarkedElements.template end<ConstrainingQuadrilateral>()){
			ConstrainingQuadrilateral* cgf = *iter;
			++iter;
			if(marked_to_normal(cgf))
				remove_hmark(cgf, RM_REFINE);
		//	if a constraining face was previously created through replacement of
		//	a constrained face, we won't process it here
			if(!marked_to_constraining(cgf))
				process_constraining_face(cgf);
		}
	}
	HNODE_PROFILE_END();

	HNODE_PROFILE_BEGIN(href_NormalFaces);
	UG_DLOG(LIB_GRID, 1, "  normal triangles.\n");
//	normal triangles
	{
		typename selector_t::template traits<Triangle>::iterator
			iter = m_selMarkedElements.template begin<Triangle>();
		while(iter != m_selMarkedElements.template end<Triangle>()){
			Face* f = *iter;
			++iter;

		//	a normal face may have previously been created by replacing a
		//	constrained or constraining face. Those faces won't be considered here
			if((!refinement_is_allowed(f)) || marked_to_normal(f)){
				continue;
			}

			if(marked_to_constraining(f))
				refine_face_with_hanging_vertex(f);
			else if(marked_refine(f) || marked_copy(f))
				refine_face_with_normal_vertex(f);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal quadrilaterals.\n");
//	normal quadrilaterals
	{
		typename selector_t::template traits<Quadrilateral>::iterator
			iter = m_selMarkedElements.template begin<Quadrilateral>();
		while(iter != m_selMarkedElements.template end<Quadrilateral>()){
			Face* f = *iter;
			++iter;

		//	a normal face may have previously been created by replacing a
		//	constrained or constraining face. Those faces won't be considered here
			if((!refinement_is_allowed(f)) || marked_to_normal(f)){
				continue;
			}

			if(marked_to_constraining(f))
				refine_face_with_hanging_vertex(f);
			else if(marked_refine(f) || marked_copy(f))
				refine_face_with_normal_vertex(f);
		}
	}
	HNODE_PROFILE_END();

////////////////////////////////
//	Volumes
	UG_DLOG(LIB_GRID, 1, "  volumes.\n");
	HNODE_PROFILE_BEGIN(href_Volumes);
	{
		typename selector_t::template traits<Volume>::iterator
			iter = m_selMarkedElements.template begin<Volume>();
		while(iter != m_selMarkedElements.template end<Volume>()){
			Volume* vol = *iter;
			++iter;
			if(refinement_is_allowed(vol)){
				refine_volume_with_normal_vertex(vol);
			}
		}
	}
	HNODE_PROFILE_END();

	UG_DLOG(LIB_GRID, 1, "  refinement done.\n");

////////////////////////////////
//	call post_refine to allow derived classes to perform some actions
	HNODE_PROFILE_BEGIN(href_PostRefine);
	post_refine();
	HNODE_PROFILE_END();



////////////////////////////////
//	before calling the refinement callback, we make sure that only elements are
//	marked, which actually received new children
	for(typename selector_t::template traits<Vertex>::iterator
			iter = m_selMarkedElements.template begin<Vertex>();
			iter != m_selMarkedElements.template end<Vertex>();)
	{
		Vertex* e = *iter;
		++iter;
		if(!(marked_refine(e) || marked_copy(e)))
			m_selMarkedElements.deselect(e);
	}
	for(typename selector_t::template traits<Edge>::iterator
			iter = m_selMarkedElements.template begin<Edge>();
			iter != m_selMarkedElements.template end<Edge>();)
	{
		Edge* e = *iter;
		++iter;
		if(!(marked_refine(e) || marked_copy(e)))
			m_selMarkedElements.deselect(e);
	}
	for(typename selector_t::template traits<Face>::iterator
			iter = m_selMarkedElements.template begin<Face>();
			iter != m_selMarkedElements.template end<Face>();)
	{
		Face* e = *iter;
		++iter;
		if(!(marked_refine(e) || marked_copy(e)))
			m_selMarkedElements.deselect(e);
	}

	// {
	// 	UG_LOG("DEBUG SAVE...\n");
	// 	static int refselCount = 0;
	// 	stringstream ss;
	// 	ss << "refselafter_" << refselCount << ".ugx";
	// 	save_marks_to_file(ss.str().c_str());
	// 	++refselCount;
	// }

//	notify the grid's message hub that refinement ends
	HNODE_PROFILE_BEGIN(href_AdaptionEndsMessage);
	m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_REFINEMENT_ENDS,
										m_selMarkedElements.get_grid_objects()));
	HNODE_PROFILE_END();

////////////////////////////////
//	Clean up
//	clear the refinement-callback if necessary
	HNODE_PROFILE_BEGIN(href_CleanUp);
	clear_marks();
	HNODE_PROFILE_END();

//	debugging utilities for the periodic boundary manager
	// if(m_pGrid->periodic_boundary_manager()){
	// 	m_pGrid->periodic_boundary_manager()->print_identification<Vertex>();
	// 	m_pGrid->periodic_boundary_manager()->print_identification<Edge>();
	// 	m_pGrid->periodic_boundary_manager()->print_identification<Face>();

	// 	UG_LOG("DEBUGGING: Checking validity of PeriodicBoundaryManager:\n");
	// 	if(m_pGrid->periodic_boundary_manager())
	// 		m_pGrid->periodic_boundary_manager()->validity_check();
	// }

	UG_DLOG(LIB_GRID, 1, "  done.\n");
}

template <class TSelector>
template <class TElem>
bool HangingNodeRefinerBase<TSelector>::remove_coarsen_marks()
{
	typedef typename selector_t::template traits<TElem>::iterator	ElemIter;
	bool removedCoarsenMark = false;
	for(ElemIter iter = m_selMarkedElements.template begin<TElem>();
		iter != m_selMarkedElements.template end<TElem>();)
	{
		TElem* e = *iter;
		++iter;
		if(get_mark(e) == RM_COARSEN){
			m_selMarkedElements.deselect(e);
			removedCoarsenMark = true;
		}
	}

	return removedCoarsenMark;
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::collect_objects_for_refine()
{
	HNODE_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "hnode_ref-start: collect_objects_for_refine\n");

	m_adjustingRefMarks = true;

//	build correct selection. see HangingVertexRefiner description.
//	unmark all elements which are marked for coarsening

	bool removedCoarseMarks = remove_coarsen_marks<Vertex>();
	removedCoarseMarks |= remove_coarsen_marks<Edge>();
	removedCoarseMarks |= remove_coarsen_marks<Face>();
	removedCoarseMarks |= remove_coarsen_marks<Volume>();

	if(removedCoarseMarks){
		UG_LOG("WARNING in HangingNodeRefinerBase::collect_objects_for_refine: "
				"Removed coarsen marks.\n");
	}

	std::vector<Vertex*>	newlyMarkedVrts;
	std::vector<Edge*>		newlyMarkedEdges;
	std::vector<Face*>			newlyMarkedFaces;
	std::vector<Volume*>		newlyMarkedVols;

	newlyMarkedVrts.assign(m_selMarkedElements.template begin<Vertex>(),
						   m_selMarkedElements.template end<Vertex>());
	newlyMarkedEdges.assign(m_selMarkedElements.template begin<Edge>(),
							m_selMarkedElements.template end<Edge>());
	newlyMarkedFaces.assign(m_selMarkedElements.template begin<Face>(),
							m_selMarkedElements.template end<Face>());
	newlyMarkedVols.assign(m_selMarkedElements.template begin<Volume>(),
						   m_selMarkedElements.template end<Volume>());

	bool continueAdjustment = true;
	bool firstAdjustment = true;

	while(continueAdjustment){
		if(!firstAdjustment){
		//	we don't simply pass m_newlyMarkedXXX to the adjusters, since we want to
		//	record newly marked elements during adjustment. Those newly marked elems
		//	are then used for the next calls, etc.
		//	This is only necessary if we're not in the first adjustment iteration.
			newlyMarkedVrts.swap(m_newlyMarkedRefVrts);
			newlyMarkedEdges.swap(m_newlyMarkedRefEdges);
			newlyMarkedFaces.swap(m_newlyMarkedRefFaces);
			newlyMarkedVols.swap(m_newlyMarkedRefVols);
		}

		firstAdjustment = false;

		m_newlyMarkedRefVrts.clear();
		m_newlyMarkedRefEdges.clear();
		m_newlyMarkedRefFaces.clear();
		m_newlyMarkedRefVols.clear();

	//	call the adjusters
		for(size_t i = 0; i < m_refMarkAdjusters.size(); ++i){
			m_refMarkAdjusters[i]->ref_marks_changed(*this, newlyMarkedVrts,
														newlyMarkedEdges,
														newlyMarkedFaces,
														newlyMarkedVols);
		}

		bool continueRequired =	(!m_newlyMarkedRefVrts.empty())
							 || (!m_newlyMarkedRefEdges.empty())
							 || (!m_newlyMarkedRefFaces.empty())
							 || (!m_newlyMarkedRefVols.empty());

		continueAdjustment = continue_collect_objects_for_refine(continueRequired);
	}

	m_adjustingRefMarks = false;
	UG_DLOG(LIB_GRID, 1, "hnode_ref-stop: collect_objects_for_refine\n");
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
assign_hnode_marks()
{
	UG_DLOG(LIB_GRID, 1, "  assigning hnode marks...\n");
//	iterate over all faces and volumes. If the element is not marked, but
//	a side is marked, the side has to be marked for hnode refinement.
//	Note that we won't mark any new elements for refinement here - we only adjust the marks
//	or add conversion marks (HNRM_TO_NORMAL etc).
//	Note also that we won't remove any marks during this algorithm (neither normal
//	nor hnode marks).
	vector<Face*> faces;
	vector<Volume*> vols;

//	the grid
	UG_ASSERT(m_pGrid, "A grid is required to perform this operation!");
	Grid& grid = *m_pGrid;

	if(grid.num<Volume>() > 0){
		for(sel_face_iter iter = m_selMarkedElements.template begin<Face>();
			iter != m_selMarkedElements.template end<Face>(); ++iter)
		{
			Face* f = *iter;
			CollectAssociated(vols, grid, f);
			for(size_t i = 0; i < vols.size(); ++i){
				if(refinement_is_allowed(vols[i])
				   && (!m_selMarkedElements.is_selected(vols[i])))
				{
					add_hmark(f, HNRM_TO_CONSTRAINING);
					break;
				}
			}
		}

	//	make sure, that constrained faces, edges and vertices of selected
	//	constraining faces will be converted to normal edges, if they are not
	//	yet marked and to constraining edges, if they are marked.
		for(typename selector_t::template traits<ConstrainingTriangle>::iterator
				iter = m_selMarkedElements.template begin<ConstrainingTriangle>();
			iter != m_selMarkedElements.template end<ConstrainingTriangle>(); ++iter)
		{
			ConstrainingTriangle* cgf = *iter;
			if(marked_refine(cgf)){
				add_hmark(cgf, HNRM_TO_NORMAL);
				for(size_t i = 0; i < cgf->num_constrained_vertices(); ++i)
					add_hmark(cgf->constrained_vertex(i), HNRM_TO_NORMAL);
				for(size_t i = 0; i < cgf->num_constrained_edges(); ++i){
					Edge* cde = cgf->constrained_edge(i);
					if(marked_refine(cde))
						add_hmark(cde, HNRM_TO_CONSTRAINING);
					else
						add_hmark(cde, HNRM_TO_NORMAL);
				}
				for(size_t i = 0; i < cgf->num_constrained_faces(); ++i){
					Face* cdf = cgf->constrained_face(i);
					if(marked_refine(cdf))
						add_hmark(cdf, HNRM_TO_CONSTRAINING);
					else
						add_hmark(cdf, HNRM_TO_NORMAL);
				}
			}
		}

		for(typename selector_t::template traits<ConstrainingQuadrilateral>::iterator
				iter = m_selMarkedElements.template begin<ConstrainingQuadrilateral>();
			iter != m_selMarkedElements.template end<ConstrainingQuadrilateral>(); ++iter)
		{
			ConstrainingQuadrilateral* cgf = *iter;
			if(marked_refine(cgf)){
				add_hmark(cgf, HNRM_TO_NORMAL);
				for(size_t i = 0; i < cgf->num_constrained_vertices(); ++i)
					add_hmark(cgf->constrained_vertex(i), HNRM_TO_NORMAL);
				for(size_t i = 0; i < cgf->num_constrained_edges(); ++i){
					Edge* cde = cgf->constrained_edge(i);
					if(marked_refine(cde))
						add_hmark(cde, HNRM_TO_CONSTRAINING);
					else
						add_hmark(cde, HNRM_TO_NORMAL);
				}
				for(size_t i = 0; i < cgf->num_constrained_faces(); ++i){
					Face* cdf = cgf->constrained_face(i);
					if(marked_refine(cdf))
						add_hmark(cdf, HNRM_TO_CONSTRAINING);
					else
						add_hmark(cdf, HNRM_TO_NORMAL);
				}
			}
		}
	}
	
	if(grid.num<Face>() > 0){
		for(sel_edge_iter iter = m_selMarkedElements.template begin<Edge>();
			iter != m_selMarkedElements.template end<Edge>(); ++iter)
		{
			Edge* e = *iter;
			CollectAssociated(faces, grid, e);
			for(size_t i = 0; i < faces.size(); ++i){
				if(marked_to_constraining(faces[i])){
					add_hmark(e, HNRM_TO_CONSTRAINING);
					break;
				}
				else if(refinement_is_allowed(faces[i])
						&& (!m_selMarkedElements.is_selected(faces[i])))
				{
					add_hmark(e, HNRM_TO_CONSTRAINING);
					break;
				}
			}
		}

		for(typename selector_t::template traits<ConstrainingEdge>::iterator
				iter = m_selMarkedElements.template begin<ConstrainingEdge>();
			iter != m_selMarkedElements.template end<ConstrainingEdge>(); ++iter)
		{
			ConstrainingEdge* cge = *iter;
			if(marked_refine(cge)){
				add_hmark(cge, HNRM_TO_NORMAL);
				for(size_t i = 0; i < cge->num_constrained_vertices(); ++i)
					add_hmark(cge->constrained_vertex(i), HNRM_TO_NORMAL);
				for(size_t i = 0; i < cge->num_constrained_edges(); ++i){
					Edge* cde = cge->constrained_edge(i);
					if(marked_refine(cde))
						add_hmark(cde, HNRM_TO_CONSTRAINING);
					else
						add_hmark(cde, HNRM_TO_NORMAL);
				}
			}
		}
	}
}

template <class TSelector>
template <class TElem>
void HangingNodeRefinerBase<TSelector>::
add_hmark(TElem* elem, HNodeRefMarks mark)
{
//	we have to consider periodic boundaries
	PeriodicBoundaryManager* pbm = m_pGrid->periodic_boundary_manager();
	typedef typename TElem::grid_base_object	base_elem_t;
	base_elem_t* e = elem;
	if(pbm && pbm->is_periodic(e)){
		base_elem_t* master = pbm->master(e);
		m_selMarkedElements.select(master,
					m_selMarkedElements.get_selection_status(master) | mark);

		typedef typename PeriodicBoundaryManager::Group<base_elem_t>::SlaveContainer
						 slave_container_t;
		slave_container_t* slaves = pbm->slaves(master);
		for(typename slave_container_t::iterator i = slaves->begin();
			i != slaves->end(); ++i)
		{
			m_selMarkedElements.select(*i,
					m_selMarkedElements.get_selection_status(*i) | mark);
		}
	}
	else{
		m_selMarkedElements.select(elem,
					m_selMarkedElements.get_selection_status(elem) | mark);
	}
}

template <class TSelector>
template <class TElem>
void HangingNodeRefinerBase<TSelector>::
remove_hmark(TElem* elem, uint mark)
{
	//	we have to consider periodic boundaries
	PeriodicBoundaryManager* pbm = m_pGrid->periodic_boundary_manager();
	typedef typename TElem::grid_base_object	base_elem_t;
	base_elem_t* e = elem;
	if(pbm && pbm->is_periodic(e)){
		base_elem_t* master = pbm->master(e);
		m_selMarkedElements.select(master,
					m_selMarkedElements.get_selection_status(master) & (~mark));

		typedef typename PeriodicBoundaryManager::Group<base_elem_t>::SlaveContainer
						 slave_container_t;
		slave_container_t* slaves = pbm->slaves(master);
		for(typename slave_container_t::iterator i = slaves->begin();
			i != slaves->end(); ++i)
		{
			m_selMarkedElements.select(*i,
					m_selMarkedElements.get_selection_status(*i) & (~mark));
		}
	}
	else{
		m_selMarkedElements.select(elem,
					m_selMarkedElements.get_selection_status(elem) & (~mark));
	}
}

////////////////////////////////////////////////////////////////////////
//	implementation of refine-methods.
template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
process_constrained_vertex(ConstrainedVertex* cdv)
{

	if(marked_to_normal(cdv)){
		Edge* parentEdge = NULL;
		Face* parentFace = NULL;
		GridObject* constrObj = cdv->get_constraining_object();

		if(constrObj){
			switch(cdv->get_parent_base_object_id()){
				case EDGE:{
						ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(constrObj);
						UG_ASSERT(cge, "Constraining edge has invalid type!");
						UG_ASSERT(marked_to_normal(cge),
								  "Constraining edge has to be converted to a normal edge!");
						cge->clear_constrained_vertices();
						parentEdge = cge;
						UG_ASSERT(get_center_vertex(cge) == cdv,
								  "Center vertex and constrained vertex do not match!");
					}break;
				case FACE:{
						ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(constrObj);
						UG_ASSERT(cgf, "Constraining face has invalid type!");
						UG_ASSERT(marked_to_normal(cgf),
								  "Constraining edge has to be converted to a normal edge!");
						cgf->clear_constrained_vertices();
						parentFace = cgf;
						UG_ASSERT((get_center_vertex(cgf) == NULL) || (get_center_vertex(cgf) == cdv),
								  "Center vertex and constrained vertex do not match!");
					}break;
				default:
					break;
			}
		}

		Grid& grid = *m_pGrid;
		Vertex* nVrt = *grid.create_and_replace<RegularVertex>(cdv);
		if(parentEdge)
			set_center_vertex(parentEdge, nVrt);
		else if(parentFace)
			set_center_vertex(parentFace, nVrt);
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
process_constrained_edge(ConstrainedEdge* cde)
{
	GridObject* constrObj = cde->get_constraining_object();
	if(constrObj && (marked_to_normal(cde) || marked_to_constraining(cde))){
		switch(cde->get_parent_base_object_id()){
			case EDGE:{
					ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(constrObj);
					UG_ASSERT(cge, "Constraining edge has invalid type!");
					UG_ASSERT(marked_to_normal(cge),
							  "Constraining edge has to be converted to a normal edge!");
					cge->clear_constrained_edges();
				}break;
			case FACE:{
					ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(constrObj);
					UG_ASSERT(cgf, "Constraining face has invalid type!");
					UG_ASSERT(marked_to_normal(cgf),
							  "Constraining face has to be converted to a normal face!");
					cgf->clear_constrained_edges();
				}break;
			default:
				break;
		}
	}

	Grid& grid = *m_pGrid;
	if(marked_to_normal(cde)){
		grid.create_and_replace<RegularEdge>(cde);
	}
	else if(marked_to_constraining(cde)){
		refine_edge_with_hanging_vertex(cde);
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
process_constraining_edge(ConstrainingEdge* cge)
{
/*	HNODE_PROFILE_FUNC();

//	make sure that there is one hanging vertex and two constrained edges.
	assert(cge->num_constrained_vertices() == 1 && "bad number of constrained vertices. There has to be exactly 1.");
	assert(cge->num_constrained_edges() == 2 && "bad number of constrained edges. There have to be exactly 2.");

//	the grid
	Grid& grid = *m_pGrid;

//	the central hanging vertex has to be transformed into a normal vertex
	ConstrainedVertex* centralHV = NULL;
	if(cge->num_constrained_vertices() > 0)
		centralHV = dynamic_cast<ConstrainedVertex*>(cge->constrained_vertex(0));

	if(!centralHV){
		UG_LOG("The central hanging vertex of a constraining edge is missing. ignoring edge.\n");
		return;
	}

//	replace the central vertex with a normal vertex
	RegularVertex* centerVrt = *grid.create_and_replace<RegularVertex>(centralHV);

//	Iterate over the constrained edges.
//	Unmarked constrained edges will be replaced by a normal edge.
//	Marked ones will be replaced by a ConstrainingEdge. Additionally
//	associated constrained edges will be created together with the
//	new central vertex.
	for(size_t i = 0; i < cge->num_constrained_edges(); ++i){
		Edge* cde = cge->constrained_edge(i);
		if(is_marked(cde)){
			refine_edge_with_hanging_vertex(cde);
		}
		else{
		//	the constrained-edge can be transformed to a normal edge
			grid.create_and_replace<RegularEdge>(cde);
		}
	}

	cge->clear_constrained_objects();
	set_center_vertex(cge, centerVrt);*/
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
refine_edge_with_normal_vertex(Edge* e, Vertex** newCornerVrts)
{
	UG_ASSERT(refinement_is_allowed(e), "Refinement of given edge not allowed!");

//	the grid
	Grid& grid = *m_pGrid;

	if((marked_copy(e) || marked_adaptive(e)) && newCornerVrts){
		grid.create_by_cloning(e, EdgeDescriptor(newCornerVrts[0], newCornerVrts[1]), e);
		return;
	}

	RegularVertex* nVrt = *grid.create<RegularVertex>(e);
	set_center_vertex(e, nVrt);

//	allow projector to calculate a new position
	if(m_projector.valid())
		m_projector->new_vertex(nVrt, e);

//	split the edge
	vector<Edge*> vEdges(2);
	e->refine(vEdges, nVrt, newCornerVrts);
	assert((vEdges.size() == 2) && "Edge::refine - produced wrong number of edges.");
	grid.register_element(vEdges[0], e);
	grid.register_element(vEdges[1], e);
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
refine_edge_with_hanging_vertex(Edge* e, Vertex** newCornerVrts)
{
	UG_ASSERT(refinement_is_allowed(e), "Refinement of given edge not allowed!");

	Grid& grid = *m_pGrid;
//	we have to insert a hanging node.
//	e has to be transformed to a constraining edge at the same time.
	assert(!ConstrainingEdge::type_match(e) && "invalid operation. e is a constraining edge.");
	//assert(!ConstrainedEdge::type_match(e) && "invalid operation. e is a constrained edge.");

	ConstrainingEdge* ce = *grid.create_and_replace<ConstrainingEdge>(e);

	ConstrainedVertex* hv = *grid.create<ConstrainedVertex>(ce);

//	allow projector to calculate a new position
	if(m_projector.valid())
		m_projector->new_vertex(hv, ce);

	set_center_vertex(ce, hv);
	hv->set_constraining_object(ce);
	ce->add_constrained_object(hv);

	hv->set_local_coordinate_1(0.5);

//	two constrained edges have to be created.
	ConstrainedEdge* nEdge[2];
	if(newCornerVrts){
		nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(newCornerVrts[0], hv), ce);
		nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, newCornerVrts[1]), ce);
	}
	else{
		nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(ce->vertex(0), hv), ce);
		nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, ce->vertex(1)), ce);
	}

	for(uint i = 0; i < 2; ++i)
	{
		ce->add_constrained_object(nEdge[i]);
		nEdge[i]->set_constraining_object(ce);
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
refine_face_with_normal_vertex(Face* f, Vertex** newCornerVrts)
{
	UG_ASSERT(refinement_is_allowed(f), "Refinement of given face not allowed!");

//UG_LOG("refine_face_with_normal_vertex\n");
	Grid& grid = *m_pGrid;

	Vertex* vNewEdgeVertices[MAX_FACE_VERTICES];
	vector<Face*>		vFaces(f->num_vertices());// heuristic

	size_t numVrts = f->num_vertices();
	size_t numEdges = f->num_edges();
	bool noEdgeVrts = true;
	for(size_t i = 0; i < numEdges; ++i){
		Edge* e = grid.get_edge(f, i);

	//	if the face is refined with a regular rule, then every edge has to have
	//	an associated center vertex
		assert(marked_adaptive(f) ||
			   (get_mark(f) == RM_REFINE && get_center_vertex(e) != NULL));

	//	assign the center vertex
		vNewEdgeVertices[i] = get_center_vertex(e);
		if(vNewEdgeVertices[i])
			noEdgeVrts = false;
	}

	if((marked_copy(f) || (marked_adaptive(f) && noEdgeVrts)) && newCornerVrts){
		FaceDescriptor desc(numVrts);
		for(size_t i = 0; i < numVrts; ++i)
			desc.set_vertex(i, newCornerVrts[i]);
		grid.create_by_cloning(f, desc, f);
		return;
	}

//	we'll perform a regular refine
	Vertex* nVrt = NULL;
	/*f->refine_regular(vFaces, &nVrt, vNewEdgeVertices, NULL,
					  RegularVertex(), newCornerVrts);*/
	f->refine(vFaces, &nVrt, vNewEdgeVertices, NULL, newCornerVrts);

//	if a new vertex has been created during refine, then register it at the grid.
	if(nVrt)
	{
		grid.register_element(nVrt, f);

	//	allow projector to calculate a new position
		if(m_projector.valid())
			m_projector->new_vertex(nVrt, f);
	}

	for(uint i = 0; i < vFaces.size(); ++i)
	{
		grid.register_element(vFaces[i], f);
	}

	set_center_vertex(f, nVrt);
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
refine_face_with_hanging_vertex(Face* f, Vertex** newCornerVrts)
{
	UG_ASSERT(refinement_is_allowed(f), "Refinement of given face not allowed!");

//UG_LOG("refine_face_with_hanging_vertex\n");
	Grid& grid = *m_pGrid;

	size_t numVrts = f->num_vertices();
/*
	vector<Edge*> 	vEdges(f->num_edges());
	vector<Vertex*> vNewEdgeVertices(f->num_edges());
	vector<Face*>		vFaces(numVrts);// heuristic

//todo: iterate over edges directly
//	collect all associated edges.
	CollectEdges(vEdges, grid, f);
	size_t numEdges = vEdges.size();

	assert(numEdges == f->num_edges() && "ERROR in RefineFaceWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(size_t i = 0; i < numEdges; ++i)
	{
		Edge* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(f, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineFaceWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		//assert((get_center_vertex(e) != NULL) && "ERROR in RefineFaceWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}
*/
	Vertex* vNewEdgeVertices[MAX_FACE_VERTICES];
	vector<Face*>		vFaces(f->num_vertices());// heuristic

	size_t numEdges = f->num_edges();
	for(size_t i = 0; i < numEdges; ++i){
		Edge* e = grid.get_edge(f, i);

	//	if the face is refined with a regular rule, then every edge has to have
	//	an associated center vertex
		assert(marked_adaptive(f) ||
				((get_mark(f) == RM_REFINE) && (get_center_vertex(e) != NULL)));

	//	assign the center vertex
		vNewEdgeVertices[i] = get_center_vertex(e);
	}

	ConstrainingFace* cgf = NULL;
	ConstrainedVertex* hv = NULL;

//	the face has to be replaced by a constraining face.
//	we'll perform a switch here depending on the number of vertices
	switch(numVrts)
	{
		case 3:
			{
			//	create the constraining triangle and replace the old face.
				cgf = *grid.create_and_replace<ConstrainingTriangle>(f);

			//	create the constrained faces.
			//	the following triangle will not be registered at the grid. Just a temporary one.
				ConstrainedTriangle constrainedTri(cgf->vertex(0),
													cgf->vertex(1),
													cgf->vertex(2));

			//	refine the constrained tri
				Vertex* tmpVrt;
				constrainedTri.refine(vFaces, &tmpVrt, vNewEdgeVertices,
									  NULL, newCornerVrts);
			}
			break;
		case 4:
			{
				cgf = *grid.create_and_replace<ConstrainingQuadrilateral>(f);

			//	a central hanging vertex is required
				hv = *grid.create<ConstrainedVertex>(cgf);

			//	allow projector to calculate a new position
				if(m_projector.valid())
					m_projector->new_vertex(hv, cgf);

				set_center_vertex(cgf, hv);
				hv->set_constraining_object(cgf);
				cgf->add_constrained_object(hv);
				hv->set_local_coordinates(0.5, 0.5);

			//	create the constrained faces.
			//	the following quadrilateral will not be registered at the grid. Just a temporary one.
				ConstrainedQuadrilateral cdf(cgf->vertex(0), cgf->vertex(1),
											 cgf->vertex(2), cgf->vertex(3));

			//	refine the constrained quad
				Vertex* tmpVrt;
				cdf.refine(vFaces, &tmpVrt, vNewEdgeVertices, hv, newCornerVrts);
			}
			break;
		default:
			assert(!"unsupported element type.");
			break;
	}

	if(cgf)
	{
	//	register the new faces
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(vFaces[i]);
			assert(cdf && "constrained face refine did produce faces which are not constrained.");
			if(cdf)
			{
				grid.register_element(cdf, cgf);
				cdf->set_constraining_object(cgf);
				cgf->add_constrained_object(cdf);
			}
		}

	//	we have to link the new constrained edges which have been auto-generated between the constrained faces.
	//	Since only edges that lie inside of the constraining face are newly created, and since only those
	//	have to be linked with the constraining face, the following algorithm will be ok for
	//	triangles and quadrilaterals.
	//	Check for each new edge-vertex, if an edge exists with the new center vertex or with it's next neighbor.
		for(size_t i = 0; i < numEdges; ++i)
		{
			if(hv)
			{
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], hv));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(cgf);
					cgf->add_constrained_object(e);
				}
			}
			else{
			//	check if a constrained edge exists between the vertex and its next neighbor
				Vertex* vNext = vNewEdgeVertices[(i + 1) % numEdges];
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], vNext));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(cgf);
					cgf->add_constrained_object(e);
				}
			}
		}
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
process_constrained_face(ConstrainedFace* cdf)
{
	GridObject* constrObj = cdf->get_constraining_object();

	if(constrObj && (marked_to_normal(cdf) || marked_to_constraining(cdf))){
		switch(cdf->get_parent_base_object_id()){
			case FACE:{
					ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(constrObj);
					UG_ASSERT(cgf, "Constraining face has invalid type!");
					UG_ASSERT(marked_to_normal(cgf),
							  "Constraining face has to be converted to a normal face!");
					cgf->clear_constrained_faces();
				}break;
			default:
				UG_ASSERT(cdf->get_constraining_object() == NULL,
						  "Invalid constraining object encountered!");
				break;
		}
	}

	Grid& grid = *m_pGrid;
	if(marked_to_normal(cdf)){
		if(cdf->num_vertices() == 3)
			grid.create_and_replace<Triangle>(cdf);
		else
			grid.create_and_replace<Quadrilateral>(cdf);
	}
	else if(marked_to_constraining(cdf)){
		refine_face_with_hanging_vertex(cdf);
	}
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
process_constraining_face(ConstrainingFace* cgf)
{
	Grid& grid = *m_pGrid;
/*
	size_t numVrts = cgf->num_vertices();

//	make sure that there is one hanging vertex and two constrained edges.
	UG_ASSERT(cgf->num_constrained_edges() == numVrts,
			 "bad number of constrained edges: " << cgf->num_constrained_edges()
			 << ". There have to be as many as vertices: " << numVrts << "."
			 << "At face with center " << GetGridObjectCenter(grid, cgf));
	UG_ASSERT(cgf->num_constrained_faces() == 4,
			  "bad number of constrained faces. There have to be exactly 4. "
			  << "At face with center " << GetGridObjectCenter(grid, cgf));

	ConstrainedVertex* centralHV = NULL;
	RegularVertex* centerVrt = NULL;

	if(numVrts == 4){
	//	the central hanging vertex has to be transformed into a normal vertex
		centralHV = NULL;
		if(cgf->num_constrained_vertices() > 0)
			centralHV = dynamic_cast<ConstrainedVertex*>(cgf->constrained_vertex(0));

	//	replace the central vertex with a normal vertex
		if(centralHV)
			centerVrt = *grid.create_and_replace<RegularVertex>(centralHV);
	}

//	Iterate over the constrained edges.
//	Unmarked constrained edges will be replaced by a normal edge.
//	Marked ones will be replaced by a ConstrainingEdge. Additionally
//	associated constrained edges will be created together with the
//	new central vertex.
	for(size_t i = 0; i < cgf->num_constrained_edges(); ++i){
		Edge* cde = cgf->constrained_edge(i);
		if(is_marked(cde)){
			refine_edge_with_hanging_vertex(cde);
		}
		else{
		//	the constrained-edge can be transformed to a normal edge
			grid.create_and_replace<RegularEdge>(cde);
		}
	}

//	iterate over the constrained faces.
//	If it is marked, we'll replace it by a constraining face and create
//	associated constrained faces.
//	if not, it will simply be transformed to a normal face.
//	To ease implementation we will transform it anyway and if required we will
//	call refine_face_with_hanging_vertex(...).
	for(size_t i = 0; i < cgf->num_constrained_faces(); ++i){
		Face* f = cgf->constrained_face(i);
		if(is_marked(f)){
		//	refine it using hanging_node_refinement.
			refine_face_with_hanging_vertex(f);
		}
		else{
		//	replace it by a normal face
			if(f->num_vertices() == 3)
				f = *grid.create_and_replace<Triangle>(f);
			else
				f = *grid.create_and_replace<Quadrilateral>(f);
		}
	}
*/
//	cgf->clear_constrained_objects();
//	cgf itself now has to be transformed to a normal face
	UG_ASSERT(marked_to_normal(cgf), "A constraining face has to be converted to"
									 " a normal face when it is refined.");
	Vertex* centerVrt = get_center_vertex(cgf);
	Face* nFace;
	if(cgf->num_vertices() == 3)
		nFace = *grid.create_and_replace<Triangle>(cgf);
	else
		nFace = *grid.create_and_replace<Quadrilateral>(cgf);

	if(centerVrt)
		set_center_vertex(nFace, centerVrt);
}

template <class TSelector>
void HangingNodeRefinerBase<TSelector>::
refine_volume_with_normal_vertex(Volume* v, Vertex** newCornerVrts)
{
	UG_ASSERT(refinement_is_allowed(v), "Refinement of given volume not allowed!");

	Grid& grid = *m_pGrid;

	//vector<Edge*> 	vEdges(v->num_edges());
	vector<Vertex*> vNewEdgeVertices(v->num_edges());
	//vector<Face*>		vFaces(v->num_faces());
	vector<Vertex*>	vNewFaceVertices(v->num_faces());
	vector<Volume*>		vVolumes(8);// heuristic
//	collect all associated edges.

	size_t numEdges = v->num_edges();
	bool noEdgeVrts = true;
	for(size_t i = 0; i < numEdges; ++i){
		Edge* e = grid.get_edge(v, i);
		vNewEdgeVertices[i] = get_center_vertex(e);
		if(vNewEdgeVertices[i])
			noEdgeVrts = false;
		UG_ASSERT(marked_adaptive(v) || vNewEdgeVertices[i],
					"In order to fully refine a volume, all edges have"
					"to contain a new vertex!");
	}

	if((marked_copy(v) || (marked_adaptive(v) && noEdgeVrts)) && newCornerVrts){
		size_t numVrts = v->num_vertices();
		VolumeDescriptor desc(numVrts);
		for(size_t i = 0; i < numVrts; ++i)
			desc.set_vertex(i, newCornerVrts[i]);
		grid.create_by_cloning(v, desc, v);
		return;
	}

	size_t numFaces = v->num_faces();
	for(size_t i = 0; i < numFaces; ++i){
		Face* f = grid.get_face(v, i);

		/*if(!VolumeContains(v, f))
		{
			UG_LOG("Grid::get_face(vol, ind) returned bad face.");
			MultiGrid* pmg = dynamic_cast<MultiGrid*>(m_pGrid);
			UG_LOG("Vol in level " << pmg->get_level(v));
			UG_LOG(", face in level " << pmg->get_level(f) << endl);
			UG_LOG("positions of volume vertices:");
			for(size_t i_c = 0; i_c < v->num_vertices(); ++i_c){
				UG_LOG(" " << GetGridObjectCenter(grid, v->vertex(i_c)));
			}
			UG_LOG("\npositions of face vertices:");
			for(size_t i_c = 0; i_c < f->num_vertices(); ++i_c){
				UG_LOG(" " << GetGridObjectCenter(grid, f->vertex(i_c)));
			}
			UG_LOG(endl);
		}*/

		if(f->num_vertices() == 3)
			vNewFaceVertices[i] = NULL;
		else{
			vNewFaceVertices[i] = get_center_vertex(f);
			UG_ASSERT(marked_adaptive(v) || vNewFaceVertices[i],
					"In order to fully refine a volume, all quadrilateral sides have"
					" to contain a new vertex!"
					<< ", vol-center " << GetGridObjectCenter(grid, v)
					<< ", face-center " << GetGridObjectCenter(grid, f));
		//todo:remove this!!!
			/*if(!vNewFaceVertices[i]){
				UG_LOG("missing face-vertex: " << GetGridObjectCenter(grid, f) << endl);
				UG_LOG("corners of face:");
				for(size_t i_c = 0; i_c < f->num_vertices(); ++i_c){
					UG_LOG(" " << GetGridObjectCenter(grid, f->vertex(i_c)));
				}
				UG_LOG(endl);
				UG_LOG("during refinement of volume: " << GetGridObjectCenter(grid, v) << endl);
			}*/
		}
	}

//	if we're performing tetrahedral refinement, we have to collect
//	the corner coordinates, so that the refinement algorithm may choose
//	the best interior diagonal.
	vector3 corners[4];
	vector3* pCorners = NULL;
	if((v->num_vertices() == 4) && m_projector.valid()){
		for(size_t i = 0; i < 4; ++i){
			corners[i] = m_projector->geometry()->pos(v->vertex(i));
		}
		pCorners = corners;
	}

//	refine the volume and register new volumes at the grid.
	Vertex* createdVrt = NULL;
	v->refine(vVolumes, &createdVrt, &vNewEdgeVertices.front(),
			  &vNewFaceVertices.front(), NULL, RegularVertex(), newCornerVrts, pCorners);

	if(createdVrt){
	//	register the new vertex
		grid.register_element(createdVrt, v);

	//	allow projector to calculate a new position
		if(m_projector.valid())
			m_projector->new_vertex(createdVrt, v);
	}

	for(uint i = 0; i < vVolumes.size(); ++i)
		grid.register_element(vVolumes[i], v);
}

template class HangingNodeRefinerBase<Selector>;
template void HangingNodeRefinerBase<Selector>::add_hmark(Vertex*, HNodeRefMarks);
template void HangingNodeRefinerBase<Selector>::add_hmark(Edge*, HNodeRefMarks);
template void HangingNodeRefinerBase<Selector>::add_hmark(Face*, HNodeRefMarks);
template void HangingNodeRefinerBase<Selector>::add_hmark(Volume*, HNodeRefMarks);

template class HangingNodeRefinerBase<MGSelector>;
template void HangingNodeRefinerBase<MGSelector>::add_hmark(Vertex*, HNodeRefMarks);
template void HangingNodeRefinerBase<MGSelector>::add_hmark(Edge*, HNodeRefMarks);
template void HangingNodeRefinerBase<MGSelector>::add_hmark(Face*, HNodeRefMarks);
template void HangingNodeRefinerBase<MGSelector>::add_hmark(Volume*, HNodeRefMarks);

}// end of namespace
