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
#include <string>
#include <sstream>

// include bridge
#include "bindings/lua/lua_user_data.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "lib_disc/domain.h"
#include "lib_disc/domain_traits.h"
#include "lib_grid/lib_grid.h"
//todo: include this in algorithms.h
#include "lib_grid/algorithms/refinement_mark_util.h"
#include "lib_grid/refinement/adaptive_regular_mg_refiner.h"
#include "lib_grid/refinement/global_fractured_media_refiner.h"
#include "lib_grid/refinement/global_multi_grid_refiner.h"
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
#include "lib_grid/refinement/ref_mark_adjusters/horizontal_anisotropy_adjuster.h"
#include "lib_grid/algorithms/subdivision/subdivision_volumes.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#include "lib_grid/parallelization/parallel_refinement/parallel_refinement.h"
#endif

#include "common/math/misc/math_util.h"

using namespace std;

namespace ug{

/**
 * \defgroup refinement_bridge Refinement Bridge
 * \ingroup domain_bridge
 * \{
 */

////////////////////////////////////////////////////////////////////////////////
///	Returns a refinement mark, given a string describing it.
/**	Valid parameters are: "refine", "coarsen"*/
static RefinementMark StringToRefinementMark(std::string markType)
{
	TrimString(markType);
	std::transform(markType.begin(), markType.end(), markType.begin(), ::tolower);
	if(markType == "refine") return RM_REFINE;
	if(markType == "coarsen") return RM_COARSEN;

	UG_THROW("StringToRefinementMark: non-supported type: "<<markType);
	return RM_NONE;
}


////////////////////////////////////////////////////////////////////////////////
///	Creates a global domain refiner.
/**	Automatically chooses whether a parallel refiner is required.*/
template <typename TDomain>
static SmartPtr<IRefiner> GlobalDomainRefiner(TDomain* dom)
{
//todo: support normal grids, too!

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			return SmartPtr<IRefiner>(
					new ParallelGlobalRefiner_MultiGrid(
							*dom->distributed_grid_manager(),
							dom->refinement_projector()));
		}
	#endif

	return SmartPtr<IRefiner>(
				new GlobalMultiGridRefiner(
						*dom->grid(),
						dom->refinement_projector()));
}

////////////////////////////////////////////////////////////////////////////////
///	Creates an adaptive hanging node domain refiner.
/**	Automatically chooses whether a parallel refiner is required.*/
template <typename TDomain>
static SmartPtr<IRefiner> HangingNodeDomainRefiner(TDomain* dom)
{
	if(!dom->is_adaptive()){
		UG_THROW("Can't create an adaptive refiner for the given domain. "
				 	   "Construct the domain with isAdaptive enabled.");
	}

//todo: support normal grids, too!
	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			return SmartPtr<IRefiner>(
					new ParallelHangingNodeRefiner_MultiGrid(
							*dom->distributed_grid_manager(),
							dom->refinement_projector()));
		}
	#endif

	return SmartPtr<IRefiner>(
				new HangingNodeRefiner_MultiGrid(
						*dom->grid(),
						dom->refinement_projector()));
}

////////////////////////////////////////////////////////////////////////////////
///	Creates an adaptive regular domain refiner.
/**	Automatically chooses whether a parallel refiner is required.*/
template <typename TDomain>
static SmartPtr<IRefiner> CreateAdaptiveRegularDomainRefiner(TDomain* dom)
{
	if(!dom->is_adaptive()){
		UG_THROW("Can't create an adaptive refiner for the given domain. "
				 	   "Construct the domain with isAdaptive enabled.");
	}

//todo: support normal grids, too!
//todo: support parallelism, too!
	/*#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			return SmartPtr<IRefiner>(new ParallelHangingNodeRefiner_MultiGrid(*dom->distributed_grid_manager()));
		}
	#endif*/

	return SmartPtr<IRefiner>(
				new AdaptiveRegularRefiner_MultiGrid(
						*dom->grid(),
						dom->refinement_projector()));
}

////////////////////////////////////////////////////////////////////////////////
///	Creates a global fractured domain refiner.
template <class TDomain>
static SmartPtr<GlobalFracturedMediaRefiner>
CreateGlobalFracturedDomainRefiner(TDomain* dom)
{
	if(!dom->is_adaptive()){
		UG_THROW("Can't create an fractured domain refiner for the given domain. "
				 	   "Construct the domain with isAdaptive enabled.");
	}

	GlobalFracturedMediaRefiner* ref = NULL;
	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			ref = new ParallelGlobalFracturedMediaRefiner(
							*dom->distributed_grid_manager(),
							dom->refinement_projector());
		}
	#endif

	if(!ref)
		ref = new GlobalFracturedMediaRefiner(
						*dom->grid(),
						dom->refinement_projector());

	ref->set_subset_handler(*dom->subset_handler());

	return SmartPtr<GlobalFracturedMediaRefiner>(ref);
}


///	Adds a horizontal-anisotropy-adjuster to the given refiner.
/** Adds a horizontal-anisotropy-adjuster to the given refiner. The provided
 * domain is used to obtain the used position attachment.
 *
 * \note	ref has to be derived from HangingNodeRefinerBase*/
template <class TDomain>
static SPIRefMarkAdjuster
AddHorizontalAnisotropyAdjuster(IRefiner* ref, TDomain* dom)
{
	HangingNodeRefiner_MultiGrid* href = dynamic_cast<HangingNodeRefiner_MultiGrid*>(ref);
	UG_COND_THROW(!href, "A horizontal anisotropy adjuster can only be added to an instance of HangingNodeRefiner_MultiGrid.");

	SmartPtr<HorizontalAnisotropyAdjuster<
			typename TDomain::position_attachment_type> >
		haa = make_sp(new HorizontalAnisotropyAdjuster<
						typename TDomain::position_attachment_type>
						(dom->position_attachment()));

	href->add_ref_mark_adjuster(haa);
	return haa;
}



// ////////////////////////////////////////////////////////////////////////////////
// /// marks face for anisotropic refinement, if it contains edges below given sizeRatio
// /**
//  * marks face and for anisotropic refinement, if it contains edges
//  * below given sizeRatio. These edges are also marked.
//  * @return true, if face has been marked for anisotropic refinement.
//  * This is the case, when one of its edges has been marked for refinement.*/
// template <class TAAPos> bool MarkIfAnisotropic(Face* f, IRefiner* refiner, number sizeRatio, TAAPos& aaPos)
// {
// 	bool marked = false;
// 	uint num_edges = f->num_edges();
// 	vector<Edge*> edges(num_edges);
// 	// collect associated edges
// 	CollectAssociated(edges, *refiner->grid(), f);

// 	//	find the shortest edge
// 	Edge* minEdge = FindShortestEdge(edges.begin(), edges.end(), aaPos);
// 	UG_ASSERT(minEdge,
// 			"Associated edges of each face have to exist at this point.");
// 	number minLen = EdgeLength(minEdge, aaPos);

// 	//	compare all associated edges of f against minEdge (even minEdge itself,
// 	//	if somebody sets edgeRatio to 1 or higher)
// 	for (uint i_edge = 0; i_edge < num_edges; ++i_edge) {
// 		Edge* e = edges[i_edge];
// 		number len = EdgeLength(e, aaPos);
// 		//	to avoid division by 0, we only consider edges with a length > 0
// 		if(len > 0) {
// 			if(minLen / len <= sizeRatio) {
// 				//	the edge will be refined
// 				refiner->mark(e);
// 				marked = true;
// 			}
// 		}
// 	}

// 	if(marked) {
// 		//	if a vertex has been marked, also mark the face, or else just a hanging
// 		//	node would be inserted.
// 		//	We do not mark other associated objects here since this would
// 		//	cause creation of a closure and would also behave differently
// 		//	in a parallel environment, compared to a serial environment.
// 		//	By using RM_ANISOTROPIC, we avoid that associated edges of
// 		//	the marked face will automatically be marked, too.
// 		refiner->mark(f, RM_ANISOTROPIC);
// 	}

// 	return marked;
// }

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements from refinement.
/**	If used in a parallel environment only elements on the calling procs
 * are marked.
 */
static void MarkForRefinement_All(SmartPtr<IRefiner> ref)
{
	PROFILE_FUNC_GROUP("grid");
	Grid* g = ref->get_associated_grid();
	if(!g){
		UG_LOG("Refiner is not registered at a grid. Aborting.\n");
		return;
	}
	ref->mark(g->vertices_begin(), g->vertices_end());
	ref->mark(g->edges_begin(), g->edges_end());
	ref->mark(g->faces_begin(), g->faces_end());
	ref->mark(g->volumes_begin(), g->volumes_end());
}

template <class TElem>
static void MarkForRefinement_AllAnisotropic(SmartPtr<IRefiner> ref)
{
	PROFILE_FUNC_GROUP("grid");
	Grid* g = ref->get_associated_grid();
	if(!g){
		UG_LOG("Refiner is not registered at a grid. Aborting.\n");
		return;
	}
	ref->mark(g->begin<TElem>(), g->end<TElem>(), RM_ANISOTROPIC);
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all vertices in the given d-dimensional sphere.
template <class TDomain>
void MarkForAdaption_VerticesInSphere(TDomain& dom, SmartPtr<IRefiner> refiner,
                                      const typename TDomain::position_type& center,
                                      number radius, std::string markType)
{
	PROFILE_FUNC_GROUP("grid");
	typedef typename TDomain::position_accessor_type	position_accessor_type;

	RefinementMark rmMark = StringToRefinementMark(markType);

//	make sure that the refiner was created for the given domain
	if(refiner->get_associated_grid() != dom.grid().get()){
		throw(UGError("ERROR in MarkForRefinement_VerticesInSphere: "
					"Refiner was not created for the specified domain. Aborting."));
	}

	Grid& grid = *refiner->get_associated_grid();
	position_accessor_type& aaPos = dom.position_accessor();

//	we'll compare against the square radius.
	number radiusSq = radius * radius;

//	we'll store associated edges, faces and volumes in those containers
	vector<Edge*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all vertices of the grid. If a vertex is inside the given sphere,
//	then we'll mark all associated elements.
	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter)
	{
		if(VecDistanceSq(center, aaPos[*iter]) <= radiusSq){
			CollectAssociated(vEdges, grid, *iter);
			CollectAssociated(vFaces, grid, *iter);
			CollectAssociated(vVols, grid, *iter);

			refiner->mark(vEdges.begin(), vEdges.end(), rmMark);
			refiner->mark(vFaces.begin(), vFaces.end(), rmMark);
			refiner->mark(vVols.begin(), vVols.end(), rmMark);
		}
	}
}

template <class TDomain>
void MarkForRefinement_VerticesInSphere(TDomain& dom, SmartPtr<IRefiner> refiner,
                                        const typename TDomain::position_type& center,
                                        number radius)
{
	MarkForAdaption_VerticesInSphere(dom, refiner, center, radius, "refine");
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements which lie completely in the given d-dimensional sphere.
/**	Valid parameters for TElem are Edge, Face, Volume.*/
template <class TDomain, class TElem>
void MarkForRefinement_ElementsInSphere(TDomain& dom, SmartPtr<IRefiner> refiner,
									const typename TDomain::position_type& center,
									number radius)
{
	PROFILE_FUNC_GROUP("grid");
	typedef typename TDomain::position_accessor_type	position_accessor_type;
	typedef typename geometry_traits<TElem>::iterator	ElemIter;

//	make sure that the refiner was created for the given domain
	if(refiner->get_associated_grid() != dom.grid().get()){
		throw(UGError("ERROR in MarkForRefinement_VerticesInCube: "
					"Refiner was not created for the specified domain. Aborting."));
	}

	Grid& grid = *refiner->get_associated_grid();
	position_accessor_type& aaPos = dom.position_accessor();

//	we'll compare against the square radius.
	number radiusSq = radius * radius;

//	we'll store associated edges, faces and volumes in those containers
	vector<Edge*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all elements of the grid. If all vertices of an element are inside
//	the given sphere, then we'll mark those elements.
	for(ElemIter iter = grid.begin<TElem>();
		iter != grid.end<TElem>(); ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	bool flag to check whether all vertices are in the sphere
		bool bInSphere = true;

	//	loop all vertices
		for(size_t i = 0; i < elem->num_vertices(); ++i)
		{
		//	check if vertex is in sphere
			if(VecDistanceSq(center, aaPos[elem->vertex(i)]) > radiusSq)
				bInSphere = false;
		}

	//	mark the element
		if(bInSphere)
			refiner->mark(elem);
	}
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements which have vertices in the given d-dimensional cube.
/**	Make sure that TAPos is an attachment of vector_t position types.*/
template <class TDomain>
void MarkForAdaption_VerticesInCube(TDomain& dom, SmartPtr<IRefiner> refiner,
									const typename TDomain::position_type& min,
									const typename TDomain::position_type& max,
									std::string markType)
{
	PROFILE_FUNC_GROUP("grid");
	typedef typename TDomain::position_type 			position_type;
	typedef typename TDomain::position_accessor_type	position_accessor_type;

	RefinementMark rmMark = StringToRefinementMark(markType);

//	make sure that the refiner was created for the given domain
	if(refiner->get_associated_grid() != dom.grid().get()){
		throw(UGError("ERROR in MarkForRefinement_VerticesInCube: "
					"Refiner was not created for the specified domain. Aborting."));
	}

	Grid& grid = *refiner->get_associated_grid();
	position_accessor_type& aaPos = dom.position_accessor();

//	we'll store associated edges, faces and volumes in those containers
	vector<Edge*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all vertices of the grid. If a vertex is inside the given cube,
//	then we'll mark all associated elements.
	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter)
	{
	//	Position
		position_type& pos = aaPos[*iter];

	//	check flag
		bool bRefine = true;

	//	check node
		for(size_t d = 0; d < pos.size(); ++d)
			if(pos[d] < min[d] || max[d] < pos[d])
				bRefine = false;

		if(bRefine)
		{
			CollectAssociated(vEdges, grid, *iter);
			CollectAssociated(vFaces, grid, *iter);
			CollectAssociated(vVols, grid, *iter);

			refiner->mark(vEdges.begin(), vEdges.end(), rmMark);
			refiner->mark(vFaces.begin(), vFaces.end(), rmMark);
			refiner->mark(vVols.begin(), vVols.end(), rmMark);
		}
	}
}

template <class TDomain>
void MarkForRefinement_VerticesInCube(TDomain& dom, SmartPtr<IRefiner> refiner,
									const typename TDomain::position_type& min,
									const typename TDomain::position_type& max)
{
	MarkForAdaption_VerticesInCube(dom, refiner, min, max, "refine");
}


///	Marks all elements for anisotropic refienment and also marks all edges > minLen.
template <class TDomain>
void MarkAnisotropic_LongEdges(TDomain& dom, IRefiner& refiner, number minLen)
{
	UG_ASSERT(dom.grid().get() == refiner.get_associated_grid(),
			  "Grids in domain and in refiner have to match!");

	typedef typename domain_traits<TDomain::dim>::element_type	elem_t;
	typedef typename MultiGrid::traits<elem_t>::iterator iter_t;

	typename TDomain::position_accessor_type aaPos = dom.position_accessor();
	MultiGrid& mg = *dom.grid();
	MultiGrid::edge_traits::secure_container	edges;
	MultiGrid::face_traits::secure_container	faces;

	number minLenSq = sq(minLen);

	for(iter_t e_iter = mg.begin<elem_t>(); e_iter != mg.end<elem_t>(); ++e_iter){
		elem_t* elem = *e_iter;
		if(mg.has_children(elem))
			continue;

	//	we'll mark all volumes as anisotropic, since we currently have to use
	//	copy elements during anisotropic refinement.
		refiner.mark(elem, RM_ANISOTROPIC);

		mg.associated_elements(edges, elem);

		for(size_t i = 0; i < edges.size(); ++i){
			if(EdgeLengthSq(edges[i], aaPos) >= minLenSq){
				refiner.mark(edges[i]);
			}
		}
	}
}


// ////////////////////////////////////////////////////////////////////////////////
// ///	Marks the long edges in anisotropic faces and faces with a big area in anisotropic volumes.
// /**
//  * THE NOT HIGHLY SUCCESSFUL IMPLEMENTATION OF THIS METHOD LEAD TO NEW INSIGHT.
//  * A NEW ANISOTROPIC REFINEMENT STRATEGY WILL BE IMPLEMENTED, WHICH WILL ALSO
//  * LEAD TO A REIMPLEMENTATION OF THIS METHOD. BEST NOT TO USE IT UNTIL THEN.
//  *
//  * The sizeRatio determines Whether a face or a volume is considered anisotropic.
//  * Make sure that the sizeRatio is in the interval [0, 1]. If the
//  * ratio of the shortest edge to an other edge falls below the given threshold,
//  * then the associated face is considered anisotropic and the longer edge will be
//  * marked. The face itself will then be marked for anisotropic refinement.
//  * The same technique is applied to volumes, this time however the ratio between
//  * face-areas is considered when judging whether a volume is anisotropic.
//  *
//  * VOLUME MARKS ARE CURRENTLY DISABLED
//  *
//  * Note that this algorithm only really works for a serial environment.
//  * \todo	improve it so that it works in parallel, too. A specialization of
//  * 			HangingNodeRefiner may be required for this.
//  *
//  * \todo	activate and improve volume marks
//  **/
// template <class TDomain>
// void MarkForRefinement_AnisotropicElements(TDomain& dom, SmartPtr<IRefiner> refiner,
// 											number sizeRatio)
// {
// 	PROFILE_FUNC_GROUP("grid");
// 	typedef typename TDomain::position_accessor_type	position_accessor_type;

// //	make sure that the refiner was created for the given domain
// 	if(refiner->get_associated_grid() != dom.grid().get()){
// 		throw(UGError("ERROR in MarkForRefinement_VerticesInCube: "
// 					"Refiner was not created for the specified domain. Aborting."));
// 	}

// //	access the grid and the position attachment
// 	Grid& grid = *refiner->get_associated_grid();
// 	position_accessor_type& aaPos = dom.position_accessor();

// //	If the grid is a multigrid, we want to avoid marking of elements, which do
// //	have children
// 	MultiGrid* pmg = dynamic_cast<MultiGrid*>(&grid);

// //	make sure that the grid automatically generates sides for each element
// 	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES)){
// 		UG_LOG("WARNING in MarkForRefinement_AnisotropicElements: "
// 				"Enabling GRIDOPT_AUTOGENERATE_SIDES.\n");
// 		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
// 	}

// //	we'll store associated edges and faces in those containers
// 	vector<Edge*> edges;
// 	vector<Face*> faces;

// //	iterate over all faces of the grid.
// 	for(FaceIterator iter = grid.begin<Face>();
// 		iter != grid.end<Face>(); ++iter)
// 	{
// 		Face* f = *iter;

// 	//	ignore faces with children
// 		if(pmg && pmg->has_children(f))
// 			continue;

// 	//	collect associated edges
// 		CollectAssociated(edges, grid, f);

// 	//	find the shortest edge
// 		Edge* minEdge = FindShortestEdge(edges.begin(), edges.end(), aaPos);
// 		UG_ASSERT(minEdge, "Associated edges of each face have to exist at this point.");
// 		number minLen = EdgeLength(minEdge, aaPos);

// 	//	compare all associated edges of f against minEdge (even minEdge itself,
// 	//	if somebody sets edgeRatio to 1 or higher)
// 		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
// 			Edge* e = edges[i_edge];
// 			number len = EdgeLength(e, aaPos);
// 		//	to avoid division by 0, we only consider edges with a length > 0
// 			if(len > 0){
// 				if(minLen / len <= sizeRatio){
// 				//	the edge will be refined
// 					refiner->mark(e);

// 				//	we'll also mark the current face, or else just a hanging
// 				//	node would be inserted.
// 				//	We do not mark other associated objects here since this would
// 				//	cause creation of a closure and would also behave differently
// 				//	in a parallel environment, compared to a serial environment.
// 				//	By using RM_ANISOTROPIC, we avoid that associated edges of
// 				//	the marked face will automatically be marked, too.
// 					refiner->mark(f, RM_ANISOTROPIC);
// 				}
// 			}
// 		}
// 	}

// //	iterate over all faces again. We have to make sure that faces which have
// //	a marked short edge are refined regular.
// //	first push all marked faces into a queue
// //	we're using grid::mark to make sure that each face lies only once on the queue.
// //	Grid::mark has nothing to do with refinement. It is just a helper for us.
// 	grid.begin_marking();

// 	queue<Face*> queFaces;
// 	for(FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>(); ++iter){
// 		queFaces.push(*iter);
// 		grid.mark(*iter);
// 	}

// 	while(!queFaces.empty()){
// 		Face* f = queFaces.front();
// 		queFaces.pop();

// 		if(pmg && pmg->has_children(f)){
// 			grid.unmark(f);
// 			continue;
// 		}

// 	//	collect associated edges
// 		CollectAssociated(edges, grid, f);
// /*
// 		if(refiner->get_mark(f) == RM_NONE){
// 			bool gotOne = false;
// 			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
// 				Edge* e = edges[i_edge];
// 				if(refiner->get_mark(e) != RM_NONE){
// 					gotOne = true;
// 					break;
// 				}
// 			}

// 			if(gotOne){
// 				for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
// 					Edge* e = edges[i_edge];
// 					if(refiner->get_mark(e) == RM_NONE){
// 						refiner->mark(e);
// 						CollectFaces(faces, grid, e);
// 						for(size_t i_face = 0; i_face < faces.size(); ++i_face){
// 							Face* nbr = faces[i_face];
// 							if(!grid.is_marked(nbr)
// 							   && (refiner->get_mark(nbr) == RM_ANISOTROPIC))
// 							{
// 								grid.mark(nbr);
// 								queFaces.push(nbr);
// 							}
// 						}
// 					}
// 				}
// 				refiner->mark(f);
// 			}
// 		}
// 		else */if(refiner->get_mark(f) == RM_ANISOTROPIC){
// 		//	find the shortest edge
// 			Edge* minEdge = FindShortestEdge(edges.begin(), edges.end(), aaPos);
// 			UG_ASSERT(minEdge, "Associated edges of each face have to exist at this point.");
// 			number minLen = EdgeLength(minEdge, aaPos);

// 		//	check if a short edge and a long edge is selected
// 			bool longEdgeSelected = false;
// 			bool shortEdgeSelected = false;

// 			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
// 				Edge* e = edges[i_edge];
// 				if(refiner->get_mark(e) == RM_NONE)
// 					continue;

// 				number len = EdgeLength(e, aaPos);
// 			//	to avoid division by 0, we only consider edges with a length > 0
// 				if(len > 0){
// 					if(minLen / len <= sizeRatio){
// 						longEdgeSelected = true;
// 					}
// 					else{
// 						shortEdgeSelected = true;
// 					}
// 				}
// 			}

// 		//	if a short edge and a long edge was selected, we'll have to mark all
// 		//	edges and push associated faces with anisotropic mark to the queue
// 			if(longEdgeSelected && shortEdgeSelected){
// 				for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
// 					Edge* e = edges[i_edge];
// 					if(refiner->get_mark(e) == RM_NONE){
// 					//	mark it and push associated anisotropic faces to the queue
// 						refiner->mark(e);
// 	//!!!

// 	if(ConstrainingEdge::type_match(e)){
// 		UG_LOG("CONSTRAINING EDGE MARKED (2)\n");
// 	}

// 						CollectFaces(faces, grid, e);
// 						for(size_t i_face = 0; i_face < faces.size(); ++i_face){
// 							Face* nbr = faces[i_face];
// 							if(!grid.is_marked(nbr)
// 							   && (refiner->get_mark(nbr) == RM_ANISOTROPIC))
// 							{
// 								grid.mark(nbr);
// 								queFaces.push(nbr);
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	//	now unmark the face
// 		grid.unmark(f);
// 	}

// 	grid.end_marking();

// //	mark unmarked neighbors of marked edges for regular refinement
// /*
// 	for(FaceIterator iter = grid.begin<Face>();
// 		iter != grid.end<Face>(); ++iter)
// 	{
// 		Face* f = *iter;

// 	//	if it is already marked, leave it as it is
// 		if(refiner->get_mark(f) != RM_NONE)
// 			continue;

// 	//	ignore faces with children
// 		if(pmg && pmg->has_children(f))
// 			continue;

// 		for(size_t i = 0; i < f->num_edges(); ++i){
// 			if(refiner->get_mark(grid.get_side(f, i)) != RM_NONE)
// 				refiner->mark(f);
// 		}
// 	}*/

// /*
// //	now that all faces are marked, we can process volumes. We consider a
// //	volume which has an anisotropic side as an anisotropic volume
// 	for(VolumeIterator iter = grid.begin<Volume>();
// 		iter != grid.end<Volume>(); ++iter)
// 	{
// 		Volume* v = *iter;

// 	//	collect associated faces
// 		CollectAssociated(faces, grid, v);

// 	//	find the smallest face
// 		Face* minFace = FindSmallestFace(faces.begin(), faces.end(), aaPos);
// 		UG_ASSERT(minFace, "Associated faces of each volume have to exist at this point.");
// 		number minArea = FaceArea(minFace, aaPos);

// 	//	compare all associated faces of v against minArea
// 		for(size_t i_face = 0; i_face < faces.size(); ++i_face){
// 			Face* f = faces[i_face];
// 			number area = FaceArea(f, aaPos);
// 		//	avoid division by 0
// 			if(area > 0){
// 				if(minArea / area <= sizeRatio){
// 				//	the face will be refined. If it is already marked, we'll
// 				//	leave it at that, to keep the anisotropy.
// 				//	If it wasn't marked, we'll mark it for full refinement
// 				//	(all anisotropic faces have already been marked).
// 					if(refiner->get_mark(f) == RM_NONE)
// 						refiner->mark(f);

// 				//	the volume itself now has to be marked, too.
// 					refiner->mark(v, RM_ANISOTROPIC);
// 				}
// 			}
// 		}
// 	}*/
// }

// ///////
// /**
//  *
//  */
// template <class TDomain>
// void MarkForRefinement_AnisotropicElements2(TDomain& dom, SmartPtr<IRefiner> refiner,
// 												number sizeRatio)
// {
// 	PROFILE_FUNC_GROUP("grid");
// 	typedef typename TDomain::position_accessor_type	position_accessor_type;

// //	make sure that the refiner was created for the given domain
// 	if(refiner->get_associated_grid() != dom.grid().get()){
// 		throw(UGError("ERROR in MarkForRefinement_VerticesInCube: "
// 					"Refiner was not created for the specified domain. Aborting."));
// 	}

// //	access the grid and the position attachment
// 	Grid& grid = *refiner->get_associated_grid();
// 	position_accessor_type& aaPos = dom.position_accessor();
// 	IRefiner* ref = refiner.get_nonconst();

// //	If the grid is a multigrid, we want to avoid marking of elements, which do
// //	have children
// 	MultiGrid* pmg = dynamic_cast<MultiGrid*>(&grid);

// //	make sure that the grid automatically generates sides for each element
// 	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES)){
// 		UG_LOG("WARNING in MarkForRefinement_AnisotropicElements: "
// 				"Enabling GRIDOPT_AUTOGENERATE_SIDES.\n");
// 		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
// 	}

// //	we'll store associated edges, faces and volumes in those containers
// 	vector<Edge*> edges;
// 	vector<Face*> faces;
// 	vector<Volume*> volumes;
// //	iterate over all faces of the grid.
// 	for(FaceIterator iter = grid.begin<Face>();
// 		iter != grid.end<Face>(); ++iter)
// 	{
// 		Face* f = *iter;
// 		// ignore faces with children
// 		if(pmg && pmg->has_children(f))
// 			continue;

// 		// if faces has been marked, store it for later marking of its neighbours
// 		if(MarkIfAnisotropic(f, ref, sizeRatio, aaPos))
// 			faces.push_back(f);
// 		else // fixme: mark for regular refinement should not be needed!
// 			refiner->mark(f);
// 	}

// // if a face is marked for anisotropic refinement (AR),
// // also mark associated volumes for AR
// 	for(vector<Face*>::iterator iter = faces.begin(); iter != faces.end(); ++iter)
// 		CollectVolumes(volumes, grid, *iter, false);

// 	refiner->mark(volumes.begin(), volumes.end(), RM_ANISOTROPIC);
// }

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements which have vertices in the given d-dimensional cube.
/**	Make sure that TAPos is an attachment of vector_t position types.*/
template <class TDomain>
void AssignSubset_VerticesInCube(TDomain& dom, 
									const typename TDomain::position_type& min,
									const typename TDomain::position_type& max, int si)
{
	typedef typename TDomain::position_type 			position_type;
	typedef typename TDomain::position_accessor_type	position_accessor_type;

	Grid& grid = *dom.grid();
	position_accessor_type& aaPos = dom.position_accessor();
	MGSubsetHandler& sh = *dom.subset_handler();

//	iterate over all vertices of the grid. If a vertex is inside the given cube,
//	then we'll mark all associated elements.
	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter)
	{
	//	Position
		position_type& pos = aaPos[*iter];

	//	check flag
		bool bInside = true;

	//	check node
		for(size_t d = 0; d < pos.size(); ++d)
			if(pos[d] < min[d] || max[d] < pos[d])
				bInside = false;

		if(bInside)
		{
			sh.assign_subset(*iter, si);
		}
	}
}


#ifdef UG_DIM_1
static number DistanceToSurfaceDomain(MathVector<1>& tpos, const Vertex* vrt,
									  Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPosSurf)
{
	UG_THROW("Not implemented");
}
#endif

#ifdef UG_DIM_2
static number DistanceToSurfaceDomain(MathVector<2>& tpos, const Edge* edge,
									  Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPosSurf)
{
	return DistancePointToLine(tpos , aaPosSurf[edge->vertex(0)], aaPosSurf[edge->vertex(1)]);
}
#endif

#ifdef UG_DIM_3
static number DistanceToSurfaceDomain(MathVector<3>& tpos, const Face* face,
									  Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPosSurf)
{

	MathVector<3> a, b, n;
	VecSubtract(a, aaPosSurf[face->vertex(1)], aaPosSurf[face->vertex(0)]);
	VecSubtract(b, aaPosSurf[face->vertex(2)], aaPosSurf[face->vertex(0)]);
	VecCross(n, a,b);

	number bc1Out, bc2Out;
	return DistancePointToTriangle(	a, bc1Out, bc2Out, 
									tpos, 
									aaPosSurf[face->vertex(0)], aaPosSurf[face->vertex(1)], aaPosSurf[face->vertex(2)], n);
}
#endif

template <class TDomain>
void MarkForRefinement_CloseToSurface(TDomain& dom, SmartPtr<IRefiner> refiner,
									  double time, size_t maxLvl,
										const TDomain& surfaceDomain)
{
	PROFILE_FUNC();
	typedef typename TDomain::grid_type TGrid;
	typedef typename TDomain::subset_handler_type TSubsetHandler;
	typedef typename TDomain::position_type TPos;
	typedef typename TDomain::position_accessor_type TAAPos;
	typedef typename domain_traits<TDomain::dim>::element_type TElem;
	typedef typename domain_traits<TDomain::dim>::side_type TSide;
	typedef typename TGrid::template traits<TElem>::iterator TIter;

	TGrid& g = *dom.grid();
	//TSubsetHandler& sh = *dom.subset_handler();
	TAAPos aaPos = dom.position_accessor();
	TAAPos aaPosSurf = surfaceDomain.position_accessor();
	const TGrid& surface = *surfaceDomain.grid();
	Grid::edge_traits::secure_container	edges;

	for(TIter iter = g.template begin<TElem>(); iter != g.template end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		size_t lvl = g.get_level(e);
		if(lvl >= maxLvl)
			continue;

		if(!g.has_children(e)){
		//	calculate the element center and the length of the longest edge
			TPos tpos = CalculateCenter(e, aaPos);
			vector3 pos;
			VecCopy(pos, tpos, 0);


		//	the minimal edge defines h for this element
			number h = numeric_limits<number>::max();
			g.associated_elements(edges, e);
			for(size_t i = 0; i < edges.size(); ++i){
				number l = EdgeLengthSq(edges[i], aaPos);
				if(l < h){
					h = l;
				}
			}
			h = sqrt(h);


			int refine = 0;

			for( typename TGrid::template traits<TSide>::const_iterator it = surface.template begin<TSide>(); 
					it != surface.template end<TSide>(); ++it){
				const TSide* side = *it;

				number distance = DistanceToSurfaceDomain(tpos, side, aaPosSurf);

				if(distance < h)
					refine = 1;
			}

			if(refine)
				refiner->mark(e);
		}
	}
}


template <class TDomain>
void MarkForRefinement_ContainsSurfaceNode(TDomain& dom, SmartPtr<IRefiner> refiner,
									  	double time, size_t maxLvl,
										const TDomain& surfaceDomain)
{
	PROFILE_FUNC();
	typedef typename TDomain::grid_type TGrid;
	typedef typename TDomain::subset_handler_type TSubsetHandler;
	typedef typename TDomain::position_type TPos;
	typedef typename TDomain::position_accessor_type TAAPos;
	typedef typename domain_traits<TDomain::dim>::element_type TElem;
	typedef typename domain_traits<TDomain::dim>::side_type TSide;
	typedef typename TGrid::template traits<TElem>::iterator TIter;

	TGrid& g = *dom.grid();
	//TSubsetHandler& sh = *dom.subset_handler();
	TAAPos aaPos = dom.position_accessor();
	TAAPos aaPosSurf = surfaceDomain.position_accessor();
	const TGrid& surface = *surfaceDomain.grid();
	Grid::edge_traits::secure_container	edges;

	for(TIter iter = g.template begin<TElem>(); iter != g.template end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		size_t lvl = g.get_level(e);
		if(lvl >= maxLvl)
			continue;

		if(!g.has_children(e)){
		//	calculate the element center and the length of the longest edge

			int refine = 0;

			for( typename TGrid::template traits<TSide>::const_iterator it = surface.template begin<TSide>(); 
					it != surface.template end<TSide>(); ++it){
				TSide* side = const_cast<TSide*>(*it);

				for(size_t co = 0; co < NumVertices(side); ++co){

					if(ContainsPoint(e, aaPosSurf[GetVertex(side,co)], aaPos))
						refine = 1;
				}
			}

			if(refine)
				refiner->mark(e);
		}
	}
}



template <class TDomain>
void MarkForRefinement_ElementsByLuaCallback(TDomain& dom, SmartPtr<IRefiner> refiner,
											 double time, size_t maxLvl,
											 const char* luaCallbackName)
{
	PROFILE_FUNC();
	typedef typename TDomain::grid_type TGrid;
	typedef typename TDomain::subset_handler_type TSubsetHandler;
	typedef typename TDomain::position_type TPos;
	typedef typename TDomain::position_accessor_type TAAPos;
	typedef typename domain_traits<TDomain::dim>::element_type TElem;
	typedef typename TGrid::template traits<TElem>::iterator TIter;

	TGrid& g = *dom.grid();
	TSubsetHandler& sh = *dom.subset_handler();
	TAAPos aaPos = dom.position_accessor();
//	Grid::edge_traits::secure_container	edges;

	LuaFunction<int, number>	callback;
//	we'll pass the following arguments: x, y, z, /*h,*/ lvl, si, time
	callback.set_lua_callback(luaCallbackName, /*7*/6);

	for(TIter iter = g.template begin<TElem>(); iter != g.template end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		size_t lvl = g.get_level(e);
		if(lvl >= maxLvl)
			continue;

		if(!g.has_children(e)){
		//	calculate the element center and the length of the longest edge
			TPos tpos = CalculateCenter(e, aaPos);
			vector3 pos;
			VecCopy(pos, tpos, 0);

			// note: this part had been inserted, but the lua callbacks rely on 6 arguments 
			//       if you need the mesh size in the lua callback you may add it AT THE END of the parameter list
			//       HOWEVER, please contact the author of this file BEFORE inserting it
			/*
		//	the longest edge defines h for this element
			number h = numeric_limits<number>::max();
			g.associated_elements(edges, e);
			for(size_t i = 0; i < edges.size(); ++i){
				number l = EdgeLengthSq(edges[i], aaPos);
				if(l < h){
					h = l;
				}
			}
			h = sqrt(h);
			*/

			int refine = 0;
			callback(refine, 6, pos.x(), pos.y(), pos.z(), /*h,*/ (number)lvl,
					 (number)sh.get_subset_index(e), (number)time);
			if(refine)
				refiner->mark(e);
		}
	}
}

template <class TDomain>
void MarkForCoarsen_ElementsByLuaCallback(TDomain& dom, SmartPtr<IRefiner> refiner,
										  double time, const char* luaCallbackName)
{
	PROFILE_FUNC();
	if(!refiner->coarsening_supported()){
		UG_LOG("WARNING in MarkForCoarsen_ElementsByLuaCallback: "
			   "Refiner doesn't support coarsening!\n");
		return;
	}

	typedef typename TDomain::grid_type TGrid;
	typedef typename TDomain::subset_handler_type TSubsetHandler;
	typedef typename TDomain::position_type TPos;
	typedef typename TDomain::position_accessor_type TAAPos;
	typedef typename domain_traits<TDomain::dim>::element_type TElem;
	typedef typename TGrid::template traits<TElem>::iterator TIter;

	TGrid& g = *dom.grid();
	TSubsetHandler& sh = *dom.subset_handler();
	TAAPos aaPos = dom.position_accessor();

	LuaFunction<int, number>	callback;
//	we'll pass the following arguments: x, y, z, lvl, si, time
	callback.set_lua_callback(luaCallbackName, 6);

	for(TIter iter = g.template begin<TElem>(); iter != g.template end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		if(!g.has_children(e)){
			int coarsen = 0;
			TPos tpos = CalculateCenter(e, aaPos);
			vector3 pos;
			VecCopy(pos, tpos, 0);
			callback(coarsen, 6, pos.x(), pos.y(), pos.z(), (number)g.get_level(e),
					 (number)sh.get_subset_index(e), (number)time);
			if(coarsen){
				refiner->mark(e, RM_COARSEN);
			}
		}
	}
}


template <class TDomain, class TSubsetHandler, class TElem>
void MarkForRefinement_ElementsInSubset(TDomain& dom, IRefiner& refiner,
										TSubsetHandler& sh, int subsetIndex)
{
	PROFILE_FUNC();
	typedef typename GridObjectCollection::traits<TElem>::iterator iterator_t;

	GridObjectCollection goc = sh.get_grid_objects_in_subset(subsetIndex);

	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
		for(iterator_t iter = goc.begin<TElem>(lvl);
			iter != goc.end<TElem>(lvl); ++iter)
		{
			refiner.mark(*iter);
		}
	}
}


template <class TDomain, class TElem>
void MarkForAdaption_ElementsContainingPoint(TDomain& dom, IRefiner& refiner,
											   number x, number y, number z,
											   std::string markType)
{
	PROFILE_FUNC();
	typedef typename TDomain::grid_type TGrid;
	typedef typename TDomain::position_type TPos;
	typedef typename TDomain::position_accessor_type TAAPos;
	typedef typename TGrid::template traits<TElem>::iterator TIter;

	RefinementMark rmMark = StringToRefinementMark(markType);

	vector3 p3(x, y, z);
	TPos p;
	VecCopy(p, p3, 0);

	TGrid& g = *dom.grid();
	TAAPos aaPos = dom.position_accessor();

	for(TIter iter = g.template begin<TElem>();
		iter != g.template end<TElem>(); ++iter)
	{
		if(ContainsPoint(*iter, p, aaPos))
			refiner.mark(*iter, rmMark);
	}
}


template <class TDomain>
void MarkForAdaption_ElementsTouchingSubset(TDomain& dom, IRefiner& refiner,
											ISubsetHandler& sh, int subsetIndex,
											std::string markType)
{
	PROFILE_FUNC();
	typedef typename TDomain::grid_type TGrid;
	typedef typename domain_traits<TDomain::dim>::element_type TElem;
	typedef typename TGrid::template traits<TElem>::iterator TIter;

	RefinementMark rmMark = StringToRefinementMark(markType);

	TGrid& g = *dom.grid();

	Grid::vertex_traits::secure_container	vrts;
	Grid::edge_traits::secure_container		edges;
	Grid::face_traits::secure_container		faces;

	for(TIter iter = g.template begin<TElem>();
		iter != g.template end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		bool gotOne = false;
		if(VERTEX < TElem::BASE_OBJECT_ID){
			g.associated_elements(vrts, e);
			for(size_t i = 0; i < vrts.size(); ++i){
				if(sh.get_subset_index(vrts[i]) == subsetIndex){
					refiner.mark(e, rmMark);
					gotOne = true;
					break;
				}
			}
			if(gotOne)
				continue;
		}
		if(EDGE < TElem::BASE_OBJECT_ID){
			g.associated_elements(edges, e);
			for(size_t i = 0; i < edges.size(); ++i){
				if(sh.get_subset_index(edges[i]) == subsetIndex){
					refiner.mark(e, rmMark);
					gotOne = true;
					break;
				}
			}
			if(gotOne)
				continue;
		}
		if(FACE < TElem::BASE_OBJECT_ID){
			g.associated_elements(faces, e);
			for(size_t i = 0; i < faces.size(); ++i){
				if(sh.get_subset_index(faces[i]) == subsetIndex){
					refiner.mark(e, rmMark);
					gotOne = true;
					break;
				}
			}
			if(gotOne)
				continue;
		}
	}
}

template <class TDomain>
void MarkForAdaption_ElementsTouchingSubsets(TDomain& dom, IRefiner& refiner,
											 const char* subsets,
											 std::string markType)
{
	ISubsetHandler& sh = *dom.subset_handler();
	vector<string> vSubsets = TokenizeTrimString(std::string(subsets));
	for(size_t i = 0; i < vSubsets.size(); ++i){
		const int si = sh.get_subset_index(vSubsets[i].c_str());
		if(si < 0)
			UG_THROW("MarkForAdaption_ElementsTouchingSubsets: Subset '"<<vSubsets[i]<<"' not found");

		MarkForAdaption_ElementsTouchingSubset(dom, refiner, sh, si, markType);
	}
}


template <class elem_t>
void MarkForRefinement_SubsetInterfaceElements(IRefiner& refiner, ISubsetHandler& sh)
{
//	ALGORITHM
//	- store max subset index of associated higher-dim elements in elem_t.
//	- AllReduce (max)
//	- Check if an associated higher-dim elem exists with a different subset index
//	- if so: set attachment value to -2
//	- AllReduce (min)
//	- Mark all elems with attachment value -2
//
//	NOTE: a serial implementation could do without the attachment by iterating
//		  over associated higher-dim elements of elem_t and directly checking
//		  for different subsets!

	AInt aSI;
	UG_COND_THROW(!refiner.grid(), "A grid has to be associated with the given refiner.");
	UG_COND_THROW(refiner.grid() != sh.grid(), "Grid of supplied refiner and subset handler have to match");
	Grid& grid = *refiner.grid();

	grid.attach_to<elem_t> (aSI);
	Grid::AttachmentAccessor<elem_t, AInt> aaSI(grid, aSI);

	SetAttachmentValues (aaSI, grid.begin<elem_t>(), grid.end<elem_t>(), -1);

	typedef typename elem_t::sideof						hdelem_t;
	typedef typename Grid::traits<elem_t>::iterator		iter_t;
	typedef typename Grid::traits<hdelem_t>::iterator	hditer_t;
	
	typename Grid::traits<elem_t>::secure_container	sides;

	for(hditer_t i_hd = grid.begin<hdelem_t>();
		i_hd != grid.end<hdelem_t>(); ++i_hd)
	{
		hdelem_t* elem = *i_hd;
		int si = sh.get_subset_index(elem);
		grid.associated_elements(sides, elem);

		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			elem_t* side = sides[i_side];
			aaSI[side] = max(aaSI[side], si);
		}
	}

	#ifdef UG_PARALLEL
		AttachmentAllReduce<elem_t>(grid, aSI, PCL_RO_MAX);
	#endif

	for(hditer_t i_hd = grid.begin<hdelem_t>();
		i_hd != grid.end<hdelem_t>(); ++i_hd)
	{
		hdelem_t* elem = *i_hd;
		int si = sh.get_subset_index(elem);
		grid.associated_elements(sides, elem);

		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			elem_t* side = sides[i_side];
			if(aaSI[side] != si){
				aaSI[side] = -2;
			}
		}
	}

	#ifdef UG_PARALLEL
		AttachmentAllReduce<elem_t>(grid, aSI, PCL_RO_MIN);
	#endif

	for(iter_t i_elem = grid.begin<elem_t>();
		i_elem != grid.end<elem_t>(); ++i_elem)
	{
		if(aaSI[*i_elem] == -2)
			refiner.mark(*i_elem);
	}

	grid.detach_from<elem_t> (aSI);
}

template <class TDomain, class TElem>
void MarkForRefinement_SubsetInterfaceElements(TDomain& dom, IRefiner& refiner)
{
	MarkForRefinement_SubsetInterfaceElements<TElem>(refiner, *dom.subset_handler());
}


template <class TDomain>
void MarkForRefinement_AnisotropicElements(TDomain& dom, IRefiner& refiner,
										 number minEdgeRatio)
{
	typedef typename domain_traits<TDomain::dim>::element_type	TElem;

	Grid& grid = *dom.grid();
	MarkForAnisotropicRefinement(grid, refiner, minEdgeRatio,
								 grid.begin<TElem>(), grid.end<TElem>(),
								 dom.position_accessor());
}


void MarkNeighborsForFullRefinement(IRefiner& refiner, bool sideNbrsOnly)
{
	refiner.mark_neighborhood(1, RM_REFINE, sideNbrsOnly);
}

void MarkNeighborsForAnisotropicRefinement(IRefiner& refiner, bool sideNbrsOnly)
{
	refiner.mark_neighborhood(1, RM_ANISOTROPIC, sideNbrsOnly);
}


template <class TDomain>
void MarkForRefinement_AnisotropicDirection (
		TDomain& dom,
		IRefiner& refiner,
		MathVector<TDomain::dim>& dir,
		number minEdgeRatio)
{
	using std::min;
	using std::max;

	typedef MathVector<TDomain::dim> 							vector_t;
	typedef typename domain_traits<TDomain::dim>::element_type	TElem;

	vector_t ndir;
	VecNormalize(ndir, dir);

	typename TDomain::position_accessor_type aaPos = dom.position_accessor();
	MultiGrid& mg = *dom.grid();
	
	MultiGrid::edge_traits::secure_container	assEdges;

	vector<Edge*> anisoEdges;
	vector<Edge*> normalEdges;

	lg_for_each_template(TElem, elem, mg){
		if(mg.has_children(elem))
			continue;

		number shortestAnisoEdgeSq = numeric_limits<number>::max();
		number longestNormalEdgeSq = 0;

	//	we'll mark all elements as closure since we use copy elements anyway
		refiner.mark(elem, RM_CLOSURE);

		anisoEdges.clear();
		normalEdges.clear();

		mg.associated_elements(assEdges, elem);
		for_each_in_vec(Edge* e, assEdges){
			vector_t edgeDir;
			VecSubtract(edgeDir, aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
			VecNormalize(edgeDir, edgeDir);
			number dot = VecDot(edgeDir, ndir);
			if((dot + SMALL >= 1) || (dot - SMALL <= -1)){
				anisoEdges.push_back(e);
				shortestAnisoEdgeSq = min(shortestAnisoEdgeSq, EdgeLengthSq(e, aaPos));
			}
			else{
				normalEdges.push_back(e);
				longestNormalEdgeSq = max(longestNormalEdgeSq, EdgeLengthSq(e, aaPos));
			}
		}end_for;

		if(longestNormalEdgeSq > 0){
			if(shortestAnisoEdgeSq / longestNormalEdgeSq <= sq(minEdgeRatio)){
			//	mark all normal edges for full refinement
				for_each_in_vec(Edge* e, normalEdges){
					refiner.mark(e, RM_REFINE);	
				}end_for;
			}
		}
	}lg_end_for;
}

template <class TDomain>
void MarkForRefinement_EdgeDirection (
		TDomain& dom,
		IRefiner& refiner,
		MathVector<TDomain::dim>& dir,
		number minDeviationAngle,
		number maxDeviationAngle,
		bool selectFlipped)
{
	MultiGrid& mg = *dom.grid();
	MarkForRefinementByDirection(	refiner,
	                             	dom.position_accessor(),
	                             	mg.begin<Edge>(),
	                             	mg.end<Edge>(),
	                             	dir,
	                             	minDeviationAngle,
	                             	maxDeviationAngle,
	                             	selectFlipped);
}

// end group refinement_bridge
/// \}

////////////////////////////////////////////////////////////////////////////////
//	REFINEMENT PROJECTORS
// template <class TDomain>
// SmartPtr<ProjectionHandler>
// DomainProjectionHandler (TDomain& dom)
// {
// 	return make_sp(new ProjectionHandler(dom.geometry3d(), dom.subset_handler()));
// }

// template <class TDomain>
// SmartPtr<IRefinementCallback>
// LinearProjectorFactory(TDomain* dom)
// {
// 	return	make_sp(new LinearProjector(
// 				make_sp(new Geometry<3, TDomain::dim>(
// 					*dom->grid(), dom->position_attachment()))));
// }

// template <class vector_t>
// static
// vector_t StdVecToMathVec(const std::vector<number>& v)
// {
// 	vector_t mv;
// 	VecSet(mv, 0);
// 	for(size_t i = 0; (i < vector_t::Size) && (i < v.size()); ++i)
// 		mv[i] = v[i];
// 	return mv;
// }

// ///	Creates a refinement projector which projects new vertices onto a sphere
// /** Specify a domain, the center of the sphere (cx, cy, cz), and the sphere's radius.
//  */
// template <class TDomain>
// SmartPtr<IRefinementCallback>
// SphereProjectorFactory(TDomain* dom, std::vector<number> center)
// {
// 	typedef SphereProjector<typename TDomain::position_attachment_type>	TRefProj;
// 	return SmartPtr<TRefProj>(
// 			new TRefProj(*dom->grid(), dom->position_attachment(),
// 						 StdVecToMathVec<typename TDomain::position_type>(center)));
// }

// ///	Creates a refinement projector which projects new vertices onto a sphere
// /** An outer radius can also be specified. Vertices outside this outer radius will
//  * be projected linear.
//  * Specify a domain, the center of the sphere (cx, cy, cz), an inner and an outer radius.
//  */
// template <class TDomain>
// SmartPtr<IRefinementCallback>
// SphericalFalloffProjectorFactory(TDomain* dom, std::vector<number> center,
// 						  number innerRadius, number outerRadius)
// {
// 	typedef SphericalFalloffProjector<typename TDomain::position_attachment_type>	TRefProj;
// 	return SmartPtr<TRefProj>(
// 			new TRefProj(*dom->grid(), dom->position_attachment(),
// 						 StdVecToMathVec<typename TDomain::position_type>(center),
// 						 innerRadius, outerRadius));
// }

// ///	Creates a refinement projector which projects new vertices onto a cylinder
// /** Specify a domain, a point on the cylinder's axis c, the direction
//  * of the axis
//  */
// template <class TDomain>
// SmartPtr<IRefinementCallback>
// CylinderProjectorFactory(TDomain* dom, std::vector<number> c, std::vector<number> axis)
// {
// 	typedef CylinderProjector<typename TDomain::position_attachment_type>	TRefProj;
// 	return SmartPtr<TRefProj>(
// 			new TRefProj(*dom->grid(), dom->position_attachment(),
// 						StdVecToMathVec<typename TDomain::position_type>(c),
// 						StdVecToMathVec<typename TDomain::position_type>(axis)));
// }

// template <class TDomain>
// SmartPtr<IRefinementCallback>
// CylindricalFalloffProjectorFactory(TDomain* dom, std::vector<number> c,
// 				  	  	  	  	   std::vector<number> a,
// 				  	  	  	  	   number innerRadius, number outerRadius)
// {
// 	typedef CylindricalFalloffProjector<typename TDomain::position_attachment_type>	TRefProj;
// 	return SmartPtr<TRefProj>(
// 			new TRefProj(*dom->grid(), dom->position_attachment(),
// 						StdVecToMathVec<typename TDomain::position_type>(c),
// 						StdVecToMathVec<typename TDomain::position_type>(a),
// 						innerRadius, outerRadius));
// }

// template <class TDomain>
// SmartPtr<IRefinementCallback>
// SubdivisionLoopProjectorFactory(TDomain* dom)
// {
// 	typedef SubdivisionLoopProjector<typename TDomain::position_attachment_type>	TRefProj;
// 	return SmartPtr<TRefProj>(
// 			new TRefProj(*dom->grid(), dom->position_attachment(),
// 						 dom->position_attachment()));
// }


namespace domain_wrappers {

}//	end of namespace

void SetTetRefinementRule(std::string ruleName)
{
	ruleName = ToLower(ruleName);
	if(ruleName.compare("standard") == 0)
		tet_rules::SetRefinementRule(tet_rules::STANDARD);
	else if(ruleName.compare("hybrid_tet_oct") == 0)
			tet_rules::SetRefinementRule(tet_rules::HYBRID_TET_OCT);
	else{
		UG_THROW("Unknown refinement rule! Known rules are: standard, hybrid_tet_oct");
	}
}

void SetSmoothSubdivisionVolumesBoundaryRefinementRule(std::string bndRefRule)
{
	bndRefRule = ToLower(bndRefRule);
	if(bndRefRule.compare("linear") == 0)
		SetBoundaryRefinementRule(LINEAR);
	else if(bndRefRule.compare("subdiv_surf_loop_scheme") == 0)
		SetBoundaryRefinementRule(SUBDIV_SURF_LOOP_SCHEME);
	else if(bndRefRule.compare("subdiv_surf_averaging_scheme") == 0)
		SetBoundaryRefinementRule(SUBDIV_SURF_AVERAGING_SCHEME);
	else if(bndRefRule.compare("subdiv_surf_butterfly_scheme") == 0)
			SetBoundaryRefinementRule(SUBDIV_SURF_BUTTERFLY_SCHEME);
	else if(bndRefRule.compare("subdiv_vol") == 0)
		SetBoundaryRefinementRule(SUBDIV_VOL);
	else
		UG_THROW("ERROR in SetBoundaryRefinementRule: Unknown boundary refinement rule! Known rules are: 'linear', 'subdiv_surf_loop_scheme', 'subdiv_surf_averaging_scheme', 'subdiv_surf_butterfly_scheme' or 'subdiv_vol'.");
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace bridge{
namespace Refinement{

/// \addtogroup refinement_bridge
/// \{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
//	register domain independent mark methods
	reg.add_function("MarkForRefinement_All", &MarkForRefinement_All, grp, "", "ref");
	reg.add_function("MarkForRefinement_AllVolumesAnisotropic", &MarkForRefinement_AllAnisotropic<Volume>, grp, "", "ref");
//	register refinement rule switch function
	reg.add_function("SetTetRefinementRule", &SetTetRefinementRule, grp, "", "refRuleName",
			"Sets the refinement rule which is used to refine tetrahedrons. Possible parameters: 'standard', 'hybrid_tet_oct");
//	register boundary refinement rule switch function for Subdivision Volumes smoothing
	reg.add_function("SetSmoothSubdivisionVolumesBoundaryRefinementRule", &SetSmoothSubdivisionVolumesBoundaryRefinementRule, grp, "", "bndRefRule",
			"Sets the boundary refinement rule used during Subdivision Volumes smoothing. Possible parameters: 'linear', 'subdiv_surf_loop_scheme', 'subdiv_surf_averaging_scheme' or 'subdiv_vol'.");
//	smooth volume/surface subdivision
	reg.add_function("ApplySmoothSubdivisionVolumesToTopLevel", (void (*)(ug::MultiGrid&, ug::MGSubsetHandler&, ug::MGSubsetHandler&, const char*)) (&ug::ApplySmoothSubdivisionVolumesToTopLevel), grp);
	reg.add_function("ApplyConstrainedSmoothSubdivisionVolumesToTopLevel", &ApplyConstrainedSmoothSubdivisionVolumesToTopLevel, grp);
	reg.add_function("ApplySmoothSubdivisionSurfacesToTopLevel2d", &ApplySmoothSubdivisionSurfacesToTopLevel2d, grp);
	reg.add_function("ApplySmoothSubdivisionSurfacesToTopLevel3d", &ApplySmoothSubdivisionSurfacesToTopLevel3d, grp);
	reg.add_function("ProjectHierarchyToLimitSubdivisionSurface", &ProjectHierarchyToLimitSubdivisionSurface, grp);
	reg.add_function("ProjectHierarchyToLimitSubdivisionVolume", &ProjectHierarchyToLimitSubdivisionVolume, grp);
	reg.add_function("TetrahedralizeHybridTetOctGrid", &TetrahedralizeHybridTetOctGrid, grp);
	reg.add_function("CheckValences", &CheckValences, grp);

	reg.add_function("MarkNeighborsForFullRefinement",
				&MarkNeighborsForFullRefinement,
				grp, "", "refiner#sideNeighborsOnly")
		.add_function("MarkNeighborsForAnisotropicRefinement",
				&MarkNeighborsForAnisotropicRefinement,
				grp, "", "refiner#sideNeighborsOnly");
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	typedef TDomain 							domain_type;
	typedef typename TDomain::position_attachment_type apos_type;

	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	refiner factory-method registration
//	Note that the refiners themselves have already been registered in lib_grid_bridge.
	reg.add_function("GlobalDomainRefiner",
					 &GlobalDomainRefiner<domain_type>, grp, "GlobalDomainRefiner", "dom");
	reg.add_function("HangingNodeDomainRefiner",
					 &HangingNodeDomainRefiner<domain_type>, grp, "HangingNodeDomainRefiner", "dom");
	reg.add_function("GlobalFracturedDomainRefiner",
					 &CreateGlobalFracturedDomainRefiner<domain_type>, grp, "GlobalFracturedDomainRefiner", "dom");
	reg.add_function("AdaptiveRegularDomainRefiner",
					 &CreateAdaptiveRegularDomainRefiner<domain_type>, grp, "AdaptiveRegularDomainRefiner", "dom");
	reg.add_function("AddHorizontalAnisotropyAdjuster",
					&AddHorizontalAnisotropyAdjuster<domain_type>, grp, "", "refiner # dom");

//	register domain dependent mark methods
	reg.add_function("MarkForRefinement_VerticesInSphere",
				&MarkForRefinement_VerticesInSphere<domain_type>, grp,
				"", "dom#refiner#center#radius")
		.add_function("MarkForAdaption_VerticesInSphere",
				&MarkForAdaption_VerticesInSphere<domain_type>, grp,
				"", "dom#refiner#center#radius#adaption_type")
		.add_function("MarkForRefinement_EdgesInSphere",
				&MarkForRefinement_ElementsInSphere<domain_type, Edge>, grp,
				"", "dom#refiner#center#radius")
		.add_function("MarkForRefinement_FacesInSphere",
				&MarkForRefinement_ElementsInSphere<domain_type, Face>, grp,
				"", "dom#refiner#center#radius")
		.add_function("MarkForRefinement_VolumesInSphere",
				&MarkForRefinement_ElementsInSphere<domain_type, Volume>, grp,
				"", "dom#refiner#center#radius")
		.add_function("MarkForRefinement_VerticesInCube",
				&MarkForRefinement_VerticesInCube<domain_type>, grp,
				"", "dom#refiner#min#max")
		.add_function("MarkForAdaption_VerticesInCube",
				&MarkForAdaption_VerticesInCube<domain_type>, grp,
				"", "dom#refiner#min#max#adaption_type")
		.add_function("MarkAnisotropic_LongEdges",
					&MarkAnisotropic_LongEdges<domain_type>, grp,
					"", "dom#refiner#maxEdgeLen")
		// .add_function("MarkForRefinement_AnisotropicElements",
		// 		&MarkForRefinement_AnisotropicElements<domain_type>, grp,
		// 		"", "dom#refiner#sizeRatio")
		// .add_function("MarkForRefinement_AnisotropicElements2",
		// 		&MarkForRefinement_AnisotropicElements2<domain_type>, grp,
		// 		"", "dom#refiner#sizeRatio")
		.add_function("MarkForRefinement_ElementsByLuaCallback",
				&MarkForRefinement_ElementsByLuaCallback<domain_type>, grp,
				"", "dom#refiner#time#callbackName")
		.add_function("MarkForCoarsen_ElementsByLuaCallback",
				&MarkForCoarsen_ElementsByLuaCallback<domain_type>, grp,
				"", "dom#refiner#time#callbackName")
		.add_function("MarkForRefinement_ElementsInSubset",
				&MarkForRefinement_ElementsInSubset<domain_type, MGSubsetHandler,
							typename domain_traits<TDomain::dim>::element_type>,
				grp, "", "dom#refiner#subsetHandler#subsetIndex")
		.add_function("MarkForRefinement_VerticesInSubset",
				&MarkForRefinement_ElementsInSubset<domain_type, MGSubsetHandler, Vertex>,
				grp, "", "dom#refiner#subsetHandler#subsetIndex")
		.add_function("MarkForRefinement_EdgesInSubset",
				&MarkForRefinement_ElementsInSubset<domain_type, MGSubsetHandler, Edge>,
				grp, "", "dom#refiner#subsetHandler#subsetIndex")
		.add_function("MarkForRefinement_FacesInSubset",
				&MarkForRefinement_ElementsInSubset<domain_type, MGSubsetHandler, Face>,
				grp, "", "dom#refiner#subsetHandler#subsetIndex")
		.add_function("MarkForRefinement_VolumesInSubset",
				&MarkForRefinement_ElementsInSubset<domain_type, MGSubsetHandler, Volume>,
				grp, "", "dom#refiner#subsetHandler#subsetIndex")
		.add_function("MarkForAdaption_ElementsTouchingSubset",
				&MarkForAdaption_ElementsTouchingSubset<domain_type>,
				grp, "", "dom#refiner#subsetHandler#subsetIndex#strMark")
		.add_function("MarkForAdaption_ElementsTouchingSubsets",
				&MarkForAdaption_ElementsTouchingSubsets<domain_type>,
				grp, "", "dom#refiner#subsetHandler#subsetNames#strMark")
		.add_function("MarkForRefinement_SubsetInterfaceVertices",
					&MarkForRefinement_SubsetInterfaceElements<TDomain, Vertex>,
					grp, "dom#refiner")
		.add_function("MarkForRefinement_SubsetInterfaceEdges",
					&MarkForRefinement_SubsetInterfaceElements<TDomain, Edge>,
					grp, "dom#refiner")
		.add_function("MarkForRefinement_SubsetInterfaceFaces",
					&MarkForRefinement_SubsetInterfaceElements<TDomain, Face>,
					grp, "dom#refiner")
		.add_function("MarkForRefinement_AnisotropicElements",
				&MarkForRefinement_AnisotropicElements<domain_type>,
				grp, "", "dom#refiner#minEdgeRatio")
		.add_function("MarkForRefinement_AnisotropicDirection",
				&MarkForRefinement_AnisotropicDirection<domain_type>,
				grp, "", "dom#refiner#dir#minEdgeRatio")
		.add_function("MarkForRefinement_EdgeDirection",
				&MarkForRefinement_EdgeDirection<domain_type>,
				grp, "", "dom#refiner#dir#minDeviationAngle#maxDeviationAngle#selectFlipped")
		.add_function("MarkForRefinement_CloseToSurface",
				&MarkForRefinement_CloseToSurface<domain_type>,
				grp, "", "")
		.add_function("MarkForRefinement_ContainsSurfaceNode",
				&MarkForRefinement_ContainsSurfaceNode<domain_type>,
				grp, "", "")
		.add_function("AssignSubset_VerticesInCube",
				&AssignSubset_VerticesInCube<domain_type>,
				grp, "", "");


//		.add_function("MarkForAdaption_EdgesContainingPoint",
//				&MarkForAdaption_ElementsContainingPoint<domain_type, Edge>,
//				grp, "", "dom#refiner#x#y#z#markType")
//		.add_function("MarkForAdaption_FacesContainingPoint",
//				&MarkForAdaption_ElementsContainingPoint<domain_type, Face>,
//				grp, "", "dom#refiner#x#y#z#markType");
//		.add_function("MarkForRefinement_VolumesContainingPoint",
//				&MarkForAdaption_ElementsContainingPoint<domain_type, Volume>,
//				grp, "", "dom#refiner#x#y#z#markType");


//	register refinement projection handler and factories
	// {
	// 	typedef RefinementProjectionHandler<apos_type> T;
	// 	string name = string("RefinementProjectionHandler").append(suffix);
	// 	reg.add_class_<T, IRefinementCallback>(name, grp)
	// 			.add_method("set_default_callback", &T::set_default_callback, grp)
	// 			.add_method("set_callback",
	// 					static_cast<void (T::*)(int, SmartPtr<IRefinementCallback>) >
	// 						(&T::set_callback), grp)
	// 			.add_method("set_callback",
	// 					static_cast<void (T::*)(std::string, SmartPtr<IRefinementCallback>) >
	// 						(&T::set_callback), grp);
	// 	reg.add_class_to_group(name, "RefinementProjectionHandler", tag);
	// }

	// reg.add_function("DomainRefinementProjectionHandler",
	// 				&DomainRefinementProjectionHandler<TDomain>, grp,
	// 				"RefinementProjectionHandler", "domain")
	// 	.add_function("LinearProjector", &LinearProjectorFactory<TDomain>, grp,
	// 				"IRefinementCallback", "domain")
	// 	.add_function("SphereProjector", &SphereProjectorFactory<TDomain>, grp,
	// 				"IRefinementCallback", "domain#centerX#centerY#centerZ#radius")
	// 	.add_function("SphericalFalloffProjector", &SphericalFalloffProjectorFactory<TDomain>, grp,
	// 				"IRefinementCallback", "domain#centerX#centerY#centerZ#innerRadius#outerRadius")
	// 	.add_function("CylinderProjector", &CylinderProjectorFactory<TDomain>, grp,
	// 				"IRefinementCallback", "domain#centerX#centerY#centerZ#axisX#axisY#axisZ")
	// 	.add_function("CylindricalFalloffProjector", &CylindricalFalloffProjectorFactory<TDomain>, grp,
	// 				"IRefinementCallback", "domain#centerX#centerY#centerZ#axisX#axisY#axisZ#innerRadius#outerRadius")
	// 	.add_function("SubdivisionLoopProjector", &SubdivisionLoopProjectorFactory<TDomain>, grp,
	// 				"IRefinementCallback", "domain");
}

}; // end Functionality

// end group refinement_bridge
/// \}

}// end Refinement

/// \addtogroup refinement_bridge
void RegisterBridge_Refinement(Registry& reg, string grp)
{
	grp.append("/Refinement");
	typedef Refinement::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace bridge
}// end of namespace ug
