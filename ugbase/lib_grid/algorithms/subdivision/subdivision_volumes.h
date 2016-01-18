/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__UG__SUBDIVISION_VOLUMES__
#define __H__UG__SUBDIVISION_VOLUMES__

#include <vector>
#include <cassert>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#include "lib_grid/algorithms/attachment_util.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/util/compol_attachment_reduce.h"
	#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#endif

#include "lib_grid/parallelization/distributed_grid.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
//	BOUNDARY REFINEMENT RULE

/// identification of boundary refinement rule to be used
enum GlobalBoundaryRefinementRule
{
	LINEAR,
	SUBDIV_LOOP,
	SUBDIV_VOL
};

/// global boundary refinement rule information switching between linear and subdivision Loop refinement
static GlobalBoundaryRefinementRule g_boundaryRefinementRule = LINEAR;

void SetBoundaryRefinementRule(GlobalBoundaryRefinementRule refRule)
{
	g_boundaryRefinementRule = refRule;
}

GlobalBoundaryRefinementRule GetBoundaryRefinementRule()
{
	return g_boundaryRefinementRule;
}

/// Function for splitting an octahedron to 4 sub-tetrahedrons
/** Recall the refinement of a tetrahedron (s. tetrahdron_rules.cpp). A tetrahedron is
 *	refined into 4 outer tetrahedrons and 4 inner tetrahedrons. After the 4 outer
 *	tetrahedrons are created the remaining inner cavity corresponds to an octahedron.
 *	This octahedron can be split into 4 tetrahedrons in 3 different ways, depending
 *	on the length of the following diagonals:
 *	Based on the original tetrahedron we look at the three diagonals between the
 *	following edge-centers: 0-5, 1-3, 2-4
 *
 *	The diagonal between edge-centers 0-5 of the tetrahedron equals
 *	a segment between vertices 1 and 3 of the octahedron
 *
 *	The diagonal between edge-centers 1-3 of the tetrahedron equals
 *	a segment between vertices 0 and 5 of the octahedron
 *
 *	the diagonal between edge-centers 2-4 of the tetrahedron equals
 *	a segment between vertices 2 and 4 of the octahedron
 *
 *	@param grid			reference to grid
 * 	@param oct			pointer to octahedron
 * 	@param parentVol	pointer to parent volume
 * 	@param vTetsOut		reference to vector with pointers to new tetrahedrons

**/
void SplitOctahedronToTetrahedrons(	Grid& grid, Octahedron* oct, Volume* parentVol,
									std::vector<Tetrahedron*>& vTetsOut)
{
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	Determine the shortest diagonal to split upon the octahedron
	int bestDiag = 2;

	number d05 = VecDistanceSq(aaPos[oct->vertex(1)], aaPos[oct->vertex(3)]);
	number d13 = VecDistanceSq(aaPos[oct->vertex(0)], aaPos[oct->vertex(5)]);
	number d   = VecDistanceSq(aaPos[oct->vertex(2)], aaPos[oct->vertex(4)]);

	if(d13 < d){
		bestDiag = 1;
		d = d13;
	}
	if(d05 < d){
		bestDiag = 0;
	}

	Tetrahedron* tet1;
	Tetrahedron* tet2;
	Tetrahedron* tet3;
	Tetrahedron* tet4;

	switch(bestDiag){

		case 0:// diag: 0-5
		//	Remark: element creation without father element specification
			tet1 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(1),
															oct->vertex(0),
															oct->vertex(4),
															oct->vertex(3)), parentVol);

			tet2 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(0),
															oct->vertex(2),
															oct->vertex(3),
															oct->vertex(1)), parentVol);

			tet3 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(4),
															oct->vertex(3),
															oct->vertex(5),
															oct->vertex(1)), parentVol);

			tet4 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(1),
															oct->vertex(5),
															oct->vertex(2),
															oct->vertex(3)), parentVol);

			vTetsOut.push_back(tet1);
			vTetsOut.push_back(tet2);
			vTetsOut.push_back(tet3);
			vTetsOut.push_back(tet4);

		//	compare case 0 tetrahedron pattern from tetrahedron_rules.cpp:
			/*
			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 1;
			inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 4;
			inds[fi++] = NUM_VERTICES + 5;	inds[fi++] = NUM_VERTICES + 0;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;
			inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 0;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 3;
			inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 5;
			*/

			break;

		case 1:// diag: 1-3
			tet1 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(1),
															oct->vertex(0),
															oct->vertex(4),
															oct->vertex(5)), parentVol);

			tet2 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(0),
															oct->vertex(2),
															oct->vertex(3),
															oct->vertex(5)), parentVol);

			tet3 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(4),
															oct->vertex(3),
															oct->vertex(5),
															oct->vertex(0)), parentVol);

			tet4 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(1),
															oct->vertex(5),
															oct->vertex(2),
															oct->vertex(0)), parentVol);

			vTetsOut.push_back(tet1);
			vTetsOut.push_back(tet2);
			vTetsOut.push_back(tet3);
			vTetsOut.push_back(tet4);

		//	compare case 1 tetrahedron pattern from tetrahedron_rules.cpp:
			/*
			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 1;
			inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 3;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 4;
			inds[fi++] = NUM_VERTICES + 5;	inds[fi++] = NUM_VERTICES + 3;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;
			inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 1;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 3;
			inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 1;
			*/

			break;

		case 2:// diag 2-4
			tet1 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(1),
															oct->vertex(4),
															oct->vertex(5),
															oct->vertex(2)), parentVol);

			tet2 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(0),
															oct->vertex(4),
															oct->vertex(1),
															oct->vertex(2)), parentVol);

			tet3 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(4),
															oct->vertex(5),
															oct->vertex(2),
															oct->vertex(3)), parentVol);

			tet4 = *grid.create<Tetrahedron>(TetrahedronDescriptor(	oct->vertex(2),
															oct->vertex(0),
															oct->vertex(4),
															oct->vertex(3)), parentVol);

			vTetsOut.push_back(tet1);
			vTetsOut.push_back(tet2);
			vTetsOut.push_back(tet3);
			vTetsOut.push_back(tet4);

		//	compare case 2 tetrahedron pattern from tetrahedron_rules.cpp:
			/*
			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 2;
			inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 4;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 2;
			inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 4;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 3;
			inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 5;

			inds[fi++] = 4;
			inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 1;
			inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;
			*/

			break;
	}

//	Erase original octahedron
	grid.erase(oct);
}


/// Projection function for subdivision volumes refinement
/** This function projects the vertices of all levels to their smooth limit
 * 	positions determined by the subdivision volumes refinement.
 *
 * 	@param mg			reference to MultiGrid
**/
void ProjectToLimitSubdivisionVolume(MultiGrid& mg)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

//	Loop all vertices of top level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* v = *vrtIter;

	//	Reposition all parents according to the position on the finest level
		while(v)
		{
			Vertex* parent = dynamic_cast<Vertex*>(mg.get_parent(v));
			if(parent){
				aaPos[parent] = aaPos[v];
			}
			v = parent;
		}
	}
}


/// Vertex repositioning function for subdivision volumes refinement
/** This function repositions a vertex to its smoothed position in the current
 * 	level determined by the subdivision volumes refinement.
 *
 * 	@param mg			reference to MultiGrid
 * 	@param markSH		reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param vol			pointer to Volume
 * 	@param aaPos		reference to VertexAttachmentAccessor accessing aPos
 * 	@param aaSmoothVolPos	reference to VertexAttachmentAccessor accessing aSmoothVolPos
 * 	@param aaNumElems	reference to VertexAttachmentAccessor accessing aNumElems
**/
void ApplySmoothSubdivisionToVolume(MultiGrid& mg, MGSubsetHandler& markSH, Volume* vol,
									Grid::VertexAttachmentAccessor<APosition>& aaPos,
									Grid::VertexAttachmentAccessor<APosition>& aaSmoothVolPos,
									Grid::VertexAttachmentAccessor<AInt>& aaNumElems)
{
//	Declare volume centroid coordinate vector
	typedef APosition::ValueType pos_type;
	pos_type p;

//	Iterate over all volume vertices, calculate and apply local centroid masks
	for(size_t i = 0; i < vol->num_vertices(); ++i)
	{
	//	Init
		Vertex* vrt = vol->vertex(i);
		VecSet(p, 0);

	//	In case of linear or subdivision Loop boundary manifold refinement:
	//	handle vertices of separating manifolds separately
		if(markSH.get_subset_index(vrt) != -1 && g_boundaryRefinementRule != SUBDIV_VOL)
		{
			continue;
		}

	//	TETRAHEDRON CASE
		if(vol->reference_object_id() == ROID_TETRAHEDRON)
		{
		//	Summate coordinates of neighbor vertices to vrt inside tetrahedron
			for(size_t j = 0; j < vol->num_vertices(); ++j)
			{
				if(j != i)
				{
					VecAdd(p, p, aaPos[vol->vertex(j)]);
				}
			}

		//	Smooth vertex position
			VecScaleAppend(aaSmoothVolPos[vrt], -1.0/16, aaPos[vrt], 17.0/48, p);
		}

	//	OCTAHEDRON CASE
		else if(vol->reference_object_id() == ROID_OCTAHEDRON)
		{
		//	Get cell-adjacent vertex
			Vertex* oppVrt = vol->vertex(vol->get_opposing_object(vrt).second);

		//	Summate coordinates of DIRECT neighbor vertices to vrt inside octahedron
			for(size_t j = 0; j < vol->num_vertices(); ++j)
			{
				if(GetVertexIndex(vol, oppVrt) == -1)
				{
					UG_THROW("ERROR in ApplySmoothSubdivisionToVolume: identified opposing vertex actually not included in current volume.");
				}

				if(j != i && j != (size_t)GetVertexIndex(vol, oppVrt))
				{
					VecAdd(p, p, aaPos[vol->vertex(j)]);
				}
			}

		//	Smooth vertex position
			VecScaleAppend(aaSmoothVolPos[vrt], 3.0/8, aaPos[vrt], 1.0/12, p, 7.0/24, aaPos[oppVrt]);
		}

	//	UNSUPPORTED VOLUME ELEMENT CASE
		else
		{
			UG_THROW("ERROR in ApplySmoothSubdivisionToVolume: Volume type not supported for subdivision volumes refinement.");
		}
	}
}


/// Vertex repositioning function for subdivision surfaces refinement (by C. Loop, 1987)
/** This function repositions a vertex to its smoothed position in the current
 * 	level determined by the subdivision surfaces refinement.
 *
 * 	@param mg				reference to MultiGrid
 * 	@param markSH			reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param vrt				pointer to Vertex
 * 	@param aaPos			reference to VertexAttachmentAccessor accessing aPos
 * 	@param aaSmoothBndPos	reference to VertexAttachmentAccessor accessing aSmoothBndPos
**/
void ApplySmoothSubdivisionToManifoldVertex(MultiGrid& mg, MGSubsetHandler& markSH, Vertex* vrt,
											Grid::VertexAttachmentAccessor<APosition>& aaPos,
											Grid::VertexAttachmentAccessor<APosition>& aaSmoothBndPos,
											Grid::VertexAttachmentAccessor<AInt> aaNumManifoldEdges)
{
	#ifdef UG_PARALLEL
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
	#endif

//	Declare centroid coordinate vector
	typedef APosition::ValueType pos_type;
	pos_type p;
	VecSet(p, 0);

//	Load subdivision surfaces rules
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	EVEN VERTEX
	if(mg.get_parent(vrt)->reference_object_id() == ROID_VERTEX)
	{
		Vertex* parentVrt = static_cast<Vertex*>(mg.get_parent(vrt));

	//	perform loop subdivision on even surface vertices
	//	first get neighbored manifold vertices
		for(Grid::AssociatedEdgeIterator eIter = mg.associated_edges_begin(parentVrt); eIter != mg.associated_edges_end(parentVrt); ++eIter)
		{
			Edge* e = *eIter;

			if(markSH.get_subset_index(e) != -1)
			{
			//	Exclude ghost and horizontal slave neighbor vertices from contributing to centroid
				#ifdef UG_PARALLEL
					if(dgm.is_ghost(e))
						continue;

					if(dgm.contains_status(e, ES_H_SLAVE))
						continue;
				#endif

				VecAdd(p, p, aaPos[GetConnectedVertex(e, parentVrt)]);
			}
		}

		number centerWgt 	= subdiv.ref_even_inner_center_weight(aaNumManifoldEdges[vrt]);
		number nbrWgt 		= subdiv.ref_even_inner_nbr_weight(aaNumManifoldEdges[vrt]);

	//	Exclude horizontal slaves of the currently smoothed vertex to avoid multiple contributions to centroid
		#ifdef UG_PARALLEL
			if(dgm.contains_status(vrt, ES_H_SLAVE))
			{
				VecScaleAppend(aaSmoothBndPos[vrt], nbrWgt, p);
				return;
			}
		#endif

		VecScaleAppend(aaSmoothBndPos[vrt], centerWgt, aaPos[parentVrt], nbrWgt, p);
	}

/*
 * Smoothing of odd vertices x
 *
 * Weights of face adjacent vertices to x: 1/8
 * Weights of edge adjacent vertices to x: 3/8
 *
			1/8
			/ \
		3/8--x--3/8
			\ /
			1/8
 *
 */

//	ODD VERTEX
	else if(mg.get_parent(vrt)->reference_object_id() == ROID_EDGE)
	{
	//	Get parent edge
		Edge* parentEdge = static_cast<Edge*>(mg.get_parent(vrt));

	//	apply loop-subdivision on inner elements
	//	get the neighbored manifold triangles
		std::vector<Face*> associatedFaces;
		std::vector<Face*> associatedManifoldFaces;

		CollectAssociated(associatedFaces, mg, parentEdge);

		for(size_t i = 0; i < associatedFaces.size(); ++i)
		{
			if(markSH.get_subset_index(associatedFaces[i]) != -1)
			{
			//	Exclude ghost and horizontal slave manifold faces
				#ifdef UG_PARALLEL
					if(dgm.is_ghost(associatedFaces[i]))
						continue;

					if(dgm.contains_status(associatedFaces[i], ES_H_SLAVE))
						continue;
				#endif

				if(associatedManifoldFaces.size() < 2)
				{
					associatedManifoldFaces.push_back(associatedFaces[i]);
				}
			}
		}

		if(associatedManifoldFaces.size() <= 2)
		{
		//	Check, if all faces are triangles
			for(size_t i = 0; i < associatedManifoldFaces.size(); ++i)
			{
				if(associatedManifoldFaces[i]->num_vertices() != 3)
				{
					UG_THROW("ERROR in ApplySmoothSubdivisionToManifoldVertex: Non triangular faces included in grid.");
				}
			}

		//	Summate centroid of face adjacent vertices
			for(size_t i = 0; i < associatedManifoldFaces.size(); ++i)
			{
				VecAdd(p, p, aaPos[GetConnectedVertex(parentEdge, associatedManifoldFaces[i])]);
			}

		//	Exclude ghost and horizontal slaves of the parent edge vertices of the currently smoothed vertex
		//	to avoid multiple contributions to the centroid of the edge adjacent vertices
			#ifdef UG_PARALLEL
				if(dgm.is_ghost(parentEdge))
				{
					//return;
				}

				if(dgm.contains_status(parentEdge, ES_H_SLAVE))
				{
					VecScaleAppend(aaSmoothBndPos[vrt], 0.125, p);
					return;
				}
			#endif

			VecScaleAppend(aaSmoothBndPos[vrt], 0.375, aaPos[parentEdge->vertex(0)], 0.375, aaPos[parentEdge->vertex(1)], 0.125, p);
		}
		else
			UG_THROW("ERROR in ApplySmoothSubdivisionToManifoldVertex: numAssociatedManifoldFaces > 2.");
	}
}

/// Function to create a smooth subdivision volumes hierarchy
/** This function transforms a linearly refined hybrid tetra-/octahedral volume
 * 	grid hierarchy into a smoothed subdivision volumes hierarchy
 * 	(s. Schaefer et al, "Smooth subdivision of tetrahedral meshes", 2004)
 *
 * 	@param mg					reference to MultiGrid
 * 	@param markSH				reference to SubsetHandler containing marked (inner) boundary manifold
**/
void ApplySmoothSubdivisionToTopLevel(MultiGrid& mg, MGSubsetHandler& markSH)
{
//	Ensure, that hybrid tet-/oct refinement is used as refinement rule for tetrahedrons
	if(tet_rules::GetRefinementRule() != tet_rules::HYBRID_TET_OCT)
		UG_THROW("ERROR in ApplySubdivisionVolumesToTopLevel: Set necessary refinement rule by SetTetRefinementRule('hybrid_tet_oct').");

//	Vertex attachments for associated number of elements, manifold edges and smooth position
	AInt aNumElems;
	AInt aNumManifoldEdges;
	APosition aSmoothVolPos;
	APosition aSmoothBndPos;

//	attach previously declared vertex attachments with initial value 0
	mg.attach_to_vertices_dv(aNumElems, 0);
	mg.attach_to_vertices_dv(aNumManifoldEdges, 0);
	mg.attach_to_vertices_dv(aSmoothVolPos, vector3(0, 0, 0));
	mg.attach_to_vertices_dv(aSmoothBndPos, vector3(0, 0, 0));

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems(mg, aNumElems);
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldEdges(mg, aNumManifoldEdges);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos(mg, aSmoothVolPos);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothBndPos(mg, aSmoothBndPos);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policies for the above attachments
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
#ifdef UG_PARALLEL
	ComPol_AttachmentReduce<VertexLayout, AInt> comPolNumElems(mg, aNumElems, PCL_RO_SUM);
	ComPol_AttachmentReduce<VertexLayout, AInt> comPolNumManifoldEdges(mg, aNumManifoldEdges, PCL_RO_SUM);
	ComPol_AttachmentReduce<VertexLayout, AVector3> comPolSmoothVolPos(mg, aSmoothVolPos, PCL_RO_SUM);
	ComPol_AttachmentReduce<VertexLayout, AVector3> comPolSmoothBndPos(mg, aSmoothBndPos, PCL_RO_SUM);

	ComPol_CopyAttachment<VertexLayout, AInt> copyNumElems(mg, aNumElems);
	ComPol_CopyAttachment<VertexLayout, AInt> copyNumManifoldEdges(mg, aNumManifoldEdges);
	ComPol_CopyAttachment<VertexLayout, AVector3> copySmoothVolPos(mg, aSmoothVolPos);
	ComPol_CopyAttachment<VertexLayout, AVector3> copySmoothBndPos(mg, aSmoothBndPos);

	pcl::InterfaceCommunicator<VertexLayout> com;
	DistributedGridManager& dgm = *mg.distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();
#endif

/*****************************************
 *
 *	(1) DETERMINE aNumElems
 *
 *****************************************/

//	Loop all volumes of top level and calculate number of volumes each vertex is contained by
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

	//	Skip ghosts
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(vol))
				continue;
		#endif

		for(size_t i = 0; i < vol->num_vertices(); ++i)
		{
			++aaNumElems[vol->vertex(i)];
		}
	}

/*****************************************
 *
 *	(2) DETERMINE aNumManifoldEdges
 *
 *****************************************/

//	Loop all edges of top level and calculate number of associated manifold edges of each vertex
	for(EdgeIterator eIter = mg.begin<Edge>(mg.top_level()); eIter != mg.end<Edge>(mg.top_level()); ++eIter)
	{
		Edge* e = *eIter;

	// 	Check, if edge is contained in subset with marked manifold elements
		if (markSH.get_subset_index(e) != -1)
		{
		//	Skip ghost and horizontal slave edges
			#ifdef UG_PARALLEL
				if(dgm.is_ghost(e))
					continue;

				if(dgm.contains_status(e, ES_H_SLAVE))
					continue;
			#endif

		   ++aaNumManifoldEdges[e->vertex(0)];
		   ++aaNumManifoldEdges[e->vertex(1)];
		}
	}

//	Manage vertex attachment communication in parallel case -> COMMUNICATE (1) aNumElems and (2) aNumManifoldEdges
#ifdef UG_PARALLEL
//	Reduce operations (sum)
	com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, comPolNumElems);
	com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, comPolNumManifoldEdges);  // TODO: ??? effectively HM = HM + 0, since HSLAVES are left out in (2)
	com.communicate();

//	Copy operations
	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, copyNumElems);
	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, copyNumManifoldEdges);
	com.communicate();

	com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, copyNumElems);			// copy to ghosts = VMASTER for step (3) APPLY SUBDIVISION
	com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, copyNumManifoldEdges);	// copy to ghosts = VMASTER for step (3) APPLY SUBDIVISION
	com.communicate();
#endif

/*****************************************
 *
 *	(3) CALCULATE aSmoothVolPos
 *
 *****************************************/

//	VOLUME GRID SMOOTHING
//	Loop all volumes of top_level
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

	//	Skip ghost volumes
		#ifdef UG_PARALLEL
	    	if(dgm.is_ghost(vol))
	    		continue;
		#endif

	    ApplySmoothSubdivisionToVolume(mg, markSH, vol, aaPos, aaSmoothVolPos, aaNumElems);
	}

//	MANIFOLD GRID SMOOTHING
//	Loop all vertices of top level
	if(g_boundaryRefinementRule != SUBDIV_VOL)
	{
		for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
		{
			Vertex* vrt = *vrtIter;

		// 	Skip vertical slave vertices as they don't have their paren on the same processor
			#ifdef UG_PARALLEL
				if(dgm.contains_status(vrt, ES_V_SLAVE))
					continue;
			#endif

		//	In case of marked manifold vertices and activated subdivision Loop refinement apply subdivision surfaces smoothing,
		//	else linear refinement
			if(markSH.get_subset_index(vrt) != -1)
			{
				if(g_boundaryRefinementRule == SUBDIV_LOOP)
				{
					ApplySmoothSubdivisionToManifoldVertex(mg, markSH, vrt, aaPos, aaSmoothBndPos, aaNumManifoldEdges);
				}
				else if(g_boundaryRefinementRule == LINEAR)
				{
					continue;
				}
				else
				{
					UG_THROW("ERROR in ApplySmoothSubdivisionToTopLevel: Unknown boundary refinement rule. Known rules are 'LINEAR', 'SUBDIV_LOOP' or 'SUBDIV_VOL'.");
				}
			}
		}
	}

//	Manage vertex attachment communication in parallel case -> COMMUNICATE (3) aSmoothVolPos and aSmoothBndPos
#ifdef UG_PARALLEL
//------------------------------------------------------------------------------
//	@ VOLUME GRID SMOOTHING
//------------------------------------------------------------------------------
	com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, comPolSmoothVolPos);
	com.communicate();

	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, copySmoothVolPos);
	com.communicate();

	com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, copySmoothVolPos);
	com.communicate();

//------------------------------------------------------------------------------
//	@ MANIFOLD GRID SMOOTHING
//------------------------------------------------------------------------------
	//TODO: WRONG ORDER PREVIOUSLY, NOW FINALLY!!! CORRECT
	//com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, comPolSmoothBndPos);
	com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, comPolSmoothBndPos);	//raus
	com.communicate();

	//com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, copySmoothBndPos);
	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, copySmoothBndPos);	//raus
	com.communicate();

	com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, copySmoothBndPos);	// stattdessen: com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, comPolSmoothBndPos); aber vorher VSLAVES die nicht VMASTER sind mit 0 initialisieren
	com.communicate();
	// dann copy auf VMASTER
#endif

/*****************************************
 *
 *	(3) APPLY SUBDIVISION (aPos = aSmooth(Vol/Bnd)Pos)
 *
 *****************************************/

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		if(g_boundaryRefinementRule == SUBDIV_VOL)
		{
		//	Scale smooth vertex position by the number of associated volume elements (SubdivisionVolumes smoothing)
			VecScale(aaSmoothVolPos[vrt],  aaSmoothVolPos[vrt], 1.0/aaNumElems[vrt]);
			VecScale(aaPos[vrt], aaSmoothVolPos[vrt], 1.0);
		}
		else if(g_boundaryRefinementRule == SUBDIV_LOOP)
		{
		//	Only in case of inner vertices
			if(markSH.get_subset_index(vrt) == -1)
			{
			//	Scale smooth vertex position by the number of associated volume elements
				VecScale(aaSmoothVolPos[vrt],  aaSmoothVolPos[vrt], 1.0/aaNumElems[vrt]);
				VecScale(aaPos[vrt], aaSmoothVolPos[vrt], 1.0);
			}
			else
			{
				VecScale(aaPos[vrt], aaSmoothBndPos[vrt], 1.0);
			}
		}
		else if(g_boundaryRefinementRule == LINEAR)
		{
			if(markSH.get_subset_index(vrt) == -1)
			{
			//	Scale smooth vertex position by the number of associated volume elements
				VecScale(aaSmoothVolPos[vrt],  aaSmoothVolPos[vrt], 1.0/aaNumElems[vrt]);
				VecScale(aaPos[vrt], aaSmoothVolPos[vrt], 1.0);
			}
			else
				continue;
		}
		else
		{
			UG_THROW("ERROR in ApplySmoothSubdivisionToTopLevel: Unknown boundary refinement rule. Known rules are 'LINEAR', 'SUBDIV_LOOP' or 'SUBDIV_VOL'.");
		}
	}

//	detach vertex attachments
	mg.detach_from_vertices(aNumElems);
	mg.detach_from_vertices(aNumManifoldEdges);
	mg.detach_from_vertices(aSmoothVolPos);
	mg.detach_from_vertices(aSmoothBndPos);
}


}//	end of namespace

#endif
