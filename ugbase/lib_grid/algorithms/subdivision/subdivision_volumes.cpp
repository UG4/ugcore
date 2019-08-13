/*
 * Copyright (c) 2014-2019:  G-CSC, Goethe University Frankfurt
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


#include  "subdivision_volumes.h"

namespace ug
{


////////////////////////////////////////////////////////////////////////////////
void SetBoundaryRefinementRule(GlobalBoundaryRefinementRule refRule)
{
	g_boundaryRefinementRule = refRule;
}


////////////////////////////////////////////////////////////////////////////////
GlobalBoundaryRefinementRule GetBoundaryRefinementRule()
{
	return g_boundaryRefinementRule;
}


////////////////////////////////////////////////////////////////////////////////
void CheckValences(MultiGrid& mg, MGSubsetHandler& markSH, const char* filename)
{
//	Output file SubsetHandler
	SubsetHandler sh(mg);

//	Attachments for associated number of elements
	AInt aEdgeValence;
	AInt aVertexValence_toc;
	AInt aVertexValence_prism;
	AInt aVertexValence_hex;

//	attach previously declared attachments with initial value 0
	mg.attach_to_edges_dv(aEdgeValence, 0);
	mg.attach_to_vertices_dv(aVertexValence_toc, 0);
	mg.attach_to_vertices_dv(aVertexValence_prism, 0);
	mg.attach_to_vertices_dv(aVertexValence_hex, 0);

//	Define attachment accessors
	Grid::EdgeAttachmentAccessor<AInt> aaEdgeValence(mg, aEdgeValence);
	Grid::VertexAttachmentAccessor<AInt> aaVertexValence_toc(mg, aVertexValence_toc);
	Grid::VertexAttachmentAccessor<AInt> aaVertexValence_prism(mg, aVertexValence_prism);
	Grid::VertexAttachmentAccessor<AInt> aaVertexValence_hex(mg, aVertexValence_hex);

//	Calculate edge valence
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* v = *vIter;

		Grid::edge_traits::secure_container edges;
		mg.associated_elements(edges, v);

		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(markSH.get_subset_index(edges[i]) != -1)
				continue;

			++aaEdgeValence[edges[i]];
		}
	}

//	Calculate vertex valence
	CalculateNumElemsVertexAttachmentInTopLevel(mg, aVertexValence_toc, aVertexValence_prism, aVertexValence_hex);

//	Assign subset indices by valence
	for(EdgeIterator eIter = mg.begin<Edge>(mg.top_level()); eIter != mg.end<Edge>(mg.top_level()); ++eIter)
	{
		Edge* e = *eIter;

		if(markSH.get_subset_index(e) != -1)
			continue;

		if(aaEdgeValence[e] < 4)
		{
			sh.assign_subset(e, aaEdgeValence[e]);
		}
	}

//	Assign subset indices by valence
	for(VertexIterator vIter = mg.begin<Vertex>(mg.top_level()); vIter != mg.end<Vertex>(mg.top_level()); ++vIter)
	{
		Vertex* v = *vIter;

		if(markSH.get_subset_index(v) != -1)
			continue;

		if(aaVertexValence_toc[v] < 4)
		{
			sh.assign_subset(v, aaVertexValence_toc[v]);
		}
	}

//	Output grid with valence SubsetHandler
	SaveGridToFile(mg, sh, filename);

//	Clean up
	mg.detach_from_edges(aEdgeValence);
	mg.detach_from_vertices(aVertexValence_toc);
	mg.detach_from_vertices(aVertexValence_prism);
	mg.detach_from_vertices(aVertexValence_hex);
}


////////////////////////////////////////////////////////////////////////////////
void PrintSubdivisionVolumesRefinementMask()
{
//	size_t n = 3;
//
//	double lin[n][n][n];
//	double m[2*n-1][2*n-1][2*n-1];
//
//	for(size_t i = 0; i < n; ++i)
//	{
//		for(size_t j = 0; j < n; ++j)
//		{
//			for(size_t k = 0; k < n; ++k)
//			{
//				lin[i][j][k] = 0.0;
//			}
//		}
//	}
//
//	for(size_t i = 0; i < 2*n-1; ++i)
//	{
//		for(size_t j = 0; j < 2*n-1; ++j)
//		{
//			for(size_t k = 0; k < 2*n-1; ++k)
//			{
//				m[i][j][k] = 0.0;
//			}
//		}
//	}
//
//	lin[0][2][0] = 1.0;
//	lin[1][1][0] = 3.0;
//	lin[1][2][0] = 3.0;
//	lin[2][0][0] = 1.0;
//	lin[2][1][0] = 3.0;
//	lin[2][2][0] = 1.0;
//
//	lin[0][1][1] = 3.0;
//	lin[0][2][1] = 3.0;
//	lin[1][0][1] = 3.0;
//	lin[1][1][1] = 6.0;
//	lin[1][2][1] = 3.0;
//	lin[2][0][1] = 3.0;
//	lin[2][1][1] = 3.0;
//
//	lin[0][0][2] = 1.0;
//	lin[0][1][2] = 3.0;
//	lin[0][2][2] = 1.0;
//	lin[1][0][2] = 3.0;
//	lin[1][1][2] = 3.0;
//	lin[2][0][2] = 1.0;
//
//	for(size_t i = 0; i < n; ++i)
//	{
//		for(size_t j = 0; j < n; ++j)
//		{
//			for(size_t k = 0; k < n; ++k)
//			{
//				for(size_t a = 0; a < n; ++a)
//				{
//					for(size_t b = 0; b < n; ++b)
//					{
//						for(size_t c = 0; c < n; ++c)
//						{
//							m[i+a][j+b][k+c] += lin[i][j][k]*lin[a][b][c];
//						}
//					}
//				}
//			}
//		}
//	}
//
//	for(size_t i = 0; i < 2*n-1; ++i)
//	{
//		for(size_t k = 0; k < 2*n-1; ++k)
//		{
//			for(size_t j = 0; j < 2*n-1; ++j)
//			{
//				UG_LOG(m[i][j][k] << "\t ");
//			}
//
//			UG_LOG(" | ");
//		}
//
//		UG_LOG(std::endl);
//	}

	const size_t n = 3;

	double lin[n][n];
	double m[2*n-1][2*n-1];

	for(size_t i = 0; i < n; ++i)
	{
		for(size_t j = 0; j < n; ++j)
		{
			lin[i][j] = 0.0;
		}
	}

	for(size_t i = 0; i < 2*n-1; ++i)
	{
		for(size_t j = 0; j < 2*n-1; ++j)
		{
			m[i][j] = 0.0;
		}
	}

	lin[0][1] = 1.0;
	lin[0][2] = 1.0;
	lin[1][0] = 1.0;
	lin[1][1] = 2.0;
	lin[1][2] = 1.0;
	lin[2][0] = 1.0;
	lin[2][1] = 1.0;

	for(size_t i = 0; i < n; ++i)
	{
		for(size_t j = 0; j < n; ++j)
		{
			for(size_t k = 0; k < n; ++k)
			{
				for(size_t l = 0; l < n; ++l)
				{
					m[i+k][j+l] += lin[i][j]*lin[k][l];
				}
			}
		}
	}

	for(size_t i = 0; i < 2*n-1; ++i)
	{
		for(size_t j = 0; j < 2*n-1; ++j)
		{
			UG_LOG(m[i][j] << "\t");
		}

		UG_LOG(std::endl);
	}
}


////////////////////////////////////////////////////////////////////////////////
void SplitOctahedronToTetrahedrons(	Grid& grid, Octahedron* oct, Volume* parentVol,
									std::vector<Tetrahedron*>& vTetsOut, int bestDiag)
{
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	Determine the shortest diagonal to split upon the octahedron
	if(bestDiag != 0 && bestDiag != 1 && bestDiag != 2)
	{
		bestDiag = 2;

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
}


////////////////////////////////////////////////////////////////////////////////
void TetrahedralizeHybridTetOctGrid(MultiGrid& mg, int bestDiag)
{
	PROFILE_FUNC_GROUP("subdivision_volumes");

	DistributedGridManager* dgm = mg.distributed_grid_manager();

	if(bestDiag != 0 && bestDiag != 1 && bestDiag != 2)
	{
		bestDiag = -1;
	}

//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

	std::vector<Tetrahedron*> vTetsOut;

//	Loop over all levels and split octahedrons
	for(size_t i = mg.num_levels()-1; i > 0; --i)
	{
	//	Loop over all octahedrons on each level
		for(VolumeIterator octIter = mg.begin<Octahedron>(i); octIter != mg.end<Octahedron>(i); ++octIter)
		{
			Octahedron* oct 	= dynamic_cast<Octahedron*>(*octIter);
			Volume* parentVol 	= dynamic_cast<Volume*>(mg.get_parent(oct));

			#ifdef UG_PARALLEL
				/*
				 * In case of parallel distribution e.g. in level 2
				 * there can be subsequently be elements which
				 * locally don’t have a parent (esp. V_SLAVES),
				 * i.e. parentVol = NULL. And elements with
				 * parent = NULL are associated to mg.level = 0.
				 * Therefore, if(dgm->is_ghost(oct)) is not
				 * sufficient.
				 */
				if(dgm->contains_status(oct, ES_V_SLAVE))
					continue;
			#endif

			SplitOctahedronToTetrahedrons(mg, oct, parentVol, vTetsOut, bestDiag);
		}
	}

//	Erase all octahedrons in multigrid
	while(mg.begin<Octahedron>() != mg.end<Octahedron>())
	{
		Octahedron* oct = *mg.begin<Octahedron>();
		mg.erase(oct);
	}
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void ProjectHierarchyToSubdivisionLimit(MultiGrid& mg, TAPosition& aPos)
{
	PROFILE_FUNC_GROUP("subdivision_volumes");

	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in ProjectHierarchyToSubdivisionLimit:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in ProjectHierarchyToSubdivisionLimit:\n"
				 "Procedure only to be used for MultiGrids with more than one level.");
	}

	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, TAPosition> comPolCopyAPosition(mg, aPos);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif

	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);

//	AttachmentCopy VSLAVE->VMASTER on mg.top_level() in case not all VMasters in toplevel don't have correct position already
	#ifdef UG_PARALLEL
	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyAPosition);
		com.communicate();
	#endif

//	Loop all levels from toplevel down to base level
	for(int lvl = (int)mg.top_level(); lvl > 0; --lvl)
	{
	//	Loop all vertices of current level and submit positions to parent vertices
		for(VertexIterator vrtIter = mg.begin<Vertex>(lvl); vrtIter != mg.end<Vertex>(lvl); ++vrtIter)
		{
			Vertex* v = *vrtIter;
			Vertex* parent = dynamic_cast<Vertex*>(mg.get_parent(v));

		//	Only, if parent vertex exists
			if(parent)
			{
				aaPos[parent] = aaPos[v];
			}
		}

	#ifdef UG_PARALLEL
	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyAPosition);
		com.communicate();
	#endif
	}
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CalculateSmoothManifoldPosInParentLevelButterflyScheme(MultiGrid& mg, TAPosition& aPos,
															MGSubsetHandler& markSH,
															MGSubsetHandler& linearManifoldSH,
															TAPosition& aSmoothBndPosOddVrt,
															AInt& aNumManifoldEdges)
{
	/*
	 * Scheme reference:
	 *
	 * D. N. Zorin, Interpolating Subdivision for Meshes with Arbitrary Topology,
	 * SIGGRAPH '96 Proceedings of the 23rd annual conference on Computer graphics
	 * and interactive techniques, 1996.
	 */

	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in CalculateSmoothManifoldPosInParentLevelButterflyScheme:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	WARNING: Parallel implementation has to be fixed
	#ifdef UG_PARALLEL
		UG_LOG("WARNING: CalculateSmoothManifoldPosInParentLevelButterflyScheme: "
					"Parallel implementation has to be fixed." << std::endl);
	#endif

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in CalculateSmoothManifoldPosInParentLevelButterflyScheme: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);
	Grid::EdgeAttachmentAccessor<TAPosition> aaSmoothBndPosOddVrt(mg, aSmoothBndPosOddVrt);
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldEdges(mg, aNumManifoldEdges);

	#ifdef UG_PARALLEL
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
	#endif

//	Declare centroid coordinate vector
	typedef typename TAPosition::ValueType position_type;
	position_type p;
	position_type q;
	VecSet(p, 0);
	VecSet(q, 0);

	/*
	 * Smoothing of odd vertices x
	 *
		 -1/16	1/8  -1/16
			  \	/ \  /
			1/2--x--1/2
			  /	\ /  \
		 -1/16	1/8  -1/16
	 *
	 */

//	Calculate smooth position for ODD vertices (EVEN vertices will be interpolated by this scheme)
	for(EdgeIterator eIter = mg.begin<Edge>(mg.top_level()-1); eIter != mg.end<Edge>(mg.top_level()-1); ++eIter)
	{
	//	Reset centroids
		VecSet(p, 0);
		VecSet(q, 0);

		Edge* e = *eIter;

	//	Skip ghost edges
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(e))
				continue;
		#endif

	//	In case of marked manifold edges, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated subdivision Butterfly refinement calculate subdivision surfaces smooth position
		if(markSH.get_subset_index(e) != -1 && linearManifoldSH.get_subset_index(e) == -1)
		{
		//	REGULAR CASE: both edge vertices are of valence 6
			if(aaNumManifoldEdges[e->vertex(0)] == 6 && aaNumManifoldEdges[e->vertex(1)] == 6)
			{
			//	perform Butterfly subdivision on odd manifold vertices
			//	get the neighbored manifold triangles
				std::vector<Face*> associatedFaces;
				std::vector<Face*> associatedButterflyFaces;
				std::vector<Face*> associatedManifoldFaces;
				std::vector<Face*> associatedButterflyManifoldFaces;

				CollectAssociated(associatedFaces, mg, e);

				for(size_t i = 0; i < associatedFaces.size(); ++i)
				{
				//	Only consider associated faces, which are marked as manifold faces
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

			//	THROW, if more then 2 associated manifold faces have been found
				if(associatedManifoldFaces.size() <= 2)
				{
				//	Check, if all faces are triangles
					for(size_t i = 0; i < associatedManifoldFaces.size(); ++i)
					{
						if(associatedManifoldFaces[i]->num_vertices() != 3)
						{
							UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme3d: "
								"Non triangular faces included in grid: " << ElementDebugInfo(mg, associatedManifoldFaces[i]));
						}
					}

				//	Summate centroid of face adjacent vertices (with corresponding weights 1/8)
					for(size_t i = 0; i < associatedManifoldFaces.size(); ++i)
					{
						VecAdd(p, p, aaPos[GetConnectedVertex(e, associatedManifoldFaces[i])]);
					}

				//	Extend original "Loop's neighborhood diamond" to include 'BUTTERFLY VERTICES' (with corresponding weights -1/16)
					for(size_t i = 0; i < associatedManifoldFaces.size(); ++i)
					{
					//	Iterate over edges of original "Loop's neighborhood diamond" to extend to Butterfly neighborhood
						for(size_t j = 0; j < associatedManifoldFaces[i]->num_edges(); ++j)
						{
						//	Clear face container
							associatedButterflyFaces.clear();
							associatedButterflyManifoldFaces.clear();

						//	Exclude edge e currently being edited
							if(j != (size_t)GetEdgeIndex(associatedManifoldFaces[i], e))
							{
							//	Collect associated Butterfly face adjacent to edge j
								GetNeighbours(associatedButterflyFaces, mg, associatedManifoldFaces[i], j);

								for(size_t k = 0; k < associatedButterflyFaces.size(); ++k)
								{
								//	Only consider associated butterfly faces, which are marked as manifold faces
									if(markSH.get_subset_index(associatedButterflyFaces[k]) != -1)
									{
									//	Exclude ghost and horizontal slave manifold faces
										#ifdef UG_PARALLEL
											if(dgm.is_ghost(associatedButterflyFaces[k]))
												continue;

											if(dgm.contains_status(associatedButterflyFaces[k], ES_H_SLAVE))
												continue;
										#endif

										if(associatedButterflyFaces[k]->num_vertices() != 3)
										{
											UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme: "
												"Non triangular faces included in grid: " << ElementDebugInfo(mg, associatedButterflyFaces[k]));
										}

										associatedButterflyManifoldFaces.push_back(associatedButterflyFaces[k]);
									}
								}

								if(associatedButterflyManifoldFaces.size() != 1)
									UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme: number of edge associated Butterfly Manifold faces != 1.");

							//	Summate centroid of butterly face adjacent vertex
								VecAdd(q, q, aaPos[GetConnectedVertex(mg.get_edge(associatedManifoldFaces[i]->edge_desc(j)), associatedButterflyManifoldFaces[0])]);
							}
						}
					}

				//	Exclude ghost and horizontal slaves of the parent edge vertices of the currently smoothed vertex
				//	to avoid multiple contributions to the centroid of the edge adjacent vertices
					#ifdef UG_PARALLEL
						if(dgm.is_ghost(e))
						{
							continue;
						}

						if(dgm.contains_status(e, ES_H_SLAVE))
						{
							VecScaleAppend(aaSmoothBndPosOddVrt[e], 0.125, p, -1.0/16, q);
							continue;
						}
					#endif

					VecScaleAppend(aaSmoothBndPosOddVrt[e], 0.5, aaPos[e->vertex(0)], 0.5, aaPos[e->vertex(1)], 0.125, p, -1.0/16, q);
				}
				else
					UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme: numAssociatedManifoldFaces > 2.");
			}

		//	IRREGULAR CASE: at least one edge vertex irregular
			if(aaNumManifoldEdges[e->vertex(0)] != 6 || aaNumManifoldEdges[e->vertex(1)] != 6)
			{
			//	Number of centroids to calculate (1 or 2, depending on the valence of the vertices of e)
				size_t numLoops;

				if((aaNumManifoldEdges[e->vertex(0)] != 6 && aaNumManifoldEdges[e->vertex(1)] == 6) ||
				   (aaNumManifoldEdges[e->vertex(0)] == 6 && aaNumManifoldEdges[e->vertex(1)] != 6))
				{
					numLoops = 1;
				}
			//	case aaNumManifoldEdges[e->vertex(0)] != 6 && aaNumManifoldEdges[e->vertex(1)] != 6
				else
					numLoops = 2;

			//	Loop centroids to calculate
				for(size_t n = 0; n < numLoops; ++n)
				{
				//	Reset centroids
					VecSet(p, 0);
					VecSet(q, 0);

					Vertex* vrt;
					Vertex* butterflyVertex;

					std::vector<Face*> associatedFaces;
					std::vector<Face*> associatedManifoldFaces;

				//	Determine smoothing case
					if(aaNumManifoldEdges[e->vertex(0)] != 6 && aaNumManifoldEdges[e->vertex(1)] == 6)
					{
						vrt = e->vertex(0);

					//	Push back e->vertex(1) as vertex s_0
						butterflyVertex = e->vertex(1);
					}
					else if (aaNumManifoldEdges[e->vertex(0)] == 6 && aaNumManifoldEdges[e->vertex(1)] != 6)
					{
						vrt = e->vertex(1);

					//	Push back e->vertex(0) as vertex s_0
						butterflyVertex = e->vertex(0);
					}
				//	case aaNumManifoldEdges[e->vertex(0)] != 6 && aaNumManifoldEdges[e->vertex(1)] != 6
					else
					{
						vrt = e->vertex(n);

					//	Push back the other vertex as vertex s_0
						butterflyVertex = e->vertex(1 - (n % 2));
					}

				//	perform Butterfly subdivision on odd manifold vertices
				//	get the neighbored manifold triangles
					CollectAssociated(associatedFaces, mg, e);

					for(size_t i = 0; i < associatedFaces.size(); ++i)
					{
					//	Only consider associated faces, which are marked as manifold faces
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
								UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme3d: "
									"Non triangular faces included in grid: " << ElementDebugInfo(mg, associatedManifoldFaces[i]));
							}
						}

					// 	Start with one triangle of "Loop's neighborhood diamond"
						Face* f = associatedManifoldFaces[0];
						size_t k = (size_t)aaNumManifoldEdges[vrt];

						number butterflyWeight = 0.0;
						std::vector<number> butterflyWeights;

						if(k == 3)
						{
							butterflyWeights.push_back(5.0/12);
							butterflyWeights.push_back(-1.0/12);
							butterflyWeights.push_back(-1.0/12);
						}

						if(k == 4)
						{
							butterflyWeights.push_back(3.0/8);
							butterflyWeights.push_back(0.0);
							butterflyWeights.push_back(-1.0/8);
							butterflyWeights.push_back(0.0);
						}

					//	Special parallel treatment for s_0 in case e is ghost or horizontal slave
						if(k != 3 && k != 4)
							butterflyWeight = 1.0/k * 7.0/4;
						else
							butterflyWeight = butterflyWeights[0];

						VecScaleAppend(p, butterflyWeight, aaPos[butterflyVertex]);

					//	Ordered traversing of associated edges of currently considered irregular vertex vrt
						for(size_t i = 1; i < k; ++i)
						{
						//	Clear face containers
							associatedFaces.clear();
							associatedManifoldFaces.clear();

						//	Get connecting edge of next butterfly vertex s_i in line to vrt
							if(f->get_opposing_object(butterflyVertex).first != EDGE)
							{
								UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme3d: "
										"Opposing object of butterfly vertex in manifold face is not an edge: "
										<< ElementDebugInfo(mg, f))
							}

							Edge* egdeOfNextFace = mg.get_edge(f, f->get_opposing_object(butterflyVertex).second);

						//	Store butterfly vertex s_i
							egdeOfNextFace->get_opposing_side(vrt, &butterflyVertex);

						//	butterfly weight s_i
							if(k != 3 && k != 4)
								butterflyWeight = 1.0/k * (1.0/4 + cos(2*i*PI/k) + 1.0/2*cos(4*i*PI/k));
							else
								butterflyWeight = butterflyWeights[i];

						//	Summate centroid of butterly vertices s_i
							VecScaleAppend(q, butterflyWeight, aaPos[butterflyVertex]);

						//	Get next face to traverse
							GetNeighbours(associatedFaces, mg, f, GetEdgeIndex(f, egdeOfNextFace));

							for(size_t j = 0; j < associatedFaces.size(); ++j)
							{
							//	Only consider associated faces, which are marked as manifold faces
								if(markSH.get_subset_index(associatedFaces[j]) != -1)
								{
								//	Exclude ghost and horizontal slave manifold faces
									#ifdef UG_PARALLEL
										if(dgm.is_ghost(associatedFaces[j]))
											continue;

										if(dgm.contains_status(associatedFaces[j], ES_H_SLAVE))
											continue;
									#endif

									associatedManifoldFaces.push_back(associatedFaces[j]);
								}
							}

							if(associatedManifoldFaces.size() != 1)
								UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme: number of edge associated Butterfly Manifold faces != 1.");

						//	Store next face in line to traverse in next iteration
							f = associatedManifoldFaces[0];
						}

					//	Exclude ghost and horizontal slaves of the parent edge vertices of the currently smoothed vertex
					//	to avoid multiple contributions to the centroid of the edge adjacent vertices
						#ifdef UG_PARALLEL
							if(dgm.is_ghost(e))
							{
								continue;
							}

							if(dgm.contains_status(e, ES_H_SLAVE))
							{
								VecScaleAppend(aaSmoothBndPosOddVrt[e], 1.0/numLoops, q);
								continue;
							}
						#endif

						VecScaleAppend(aaSmoothBndPosOddVrt[e], 1.0/numLoops*3.0/4, aaPos[vrt], 1.0/numLoops, p, 1.0/numLoops, q);
					}
					else
						UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelButterflyScheme: numAssociatedManifoldFaces > 2.");
				}
			}
		}
	}

//	Manage vertex and edge attachment communication in parallel case -> COMMUNICATE aSmoothBndPosEvenVrt, aSmoothBndPosOddVrt
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Edge>(mg, aSmoothBndPosOddVrt, PCL_RO_SUM);
	#endif
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CalculateSmoothManifoldPosInParentLevelLoopScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
											 	 	   MGSubsetHandler& linearManifoldSH,
													   TAPosition& aSmoothBndPosEvenVrt,
													   TAPosition& aSmoothBndPosOddVrt,
													   AInt& aNumManifoldEdges)
{
	/*
	 * Scheme reference:
	 *
	 * C. Loop, Smooth subdivision surfaces based on triangles,
	 * master’s thesis, University of Utah, 1987.
	 */

	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in CalculateSmoothManifoldPosInParentLevelLoopScheme:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in CalculateSmoothManifoldPosInParentLevelLoopScheme: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);
	Grid::VertexAttachmentAccessor<TAPosition> aaSmoothBndPosEvenVrt(mg, aSmoothBndPosEvenVrt);
	Grid::EdgeAttachmentAccessor<TAPosition> aaSmoothBndPosOddVrt(mg, aSmoothBndPosOddVrt);
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldEdges(mg, aNumManifoldEdges);

	#ifdef UG_PARALLEL
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
	#endif

//	Declare centroid coordinate vector
	typedef typename TAPosition::ValueType position_type;
	position_type p;
	VecSet(p, 0);

//	Load subdivision surfaces rules
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	Calculate smooth position for EVEN vertices
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()-1); vrtIter != mg.end<Vertex>(mg.top_level()-1); ++vrtIter)
	{
		VecSet(p, 0);

		Vertex* vrt = *vrtIter;

	//	Skip ghost vertices
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(vrt))
				continue;
		#endif

	//	In case of marked manifold vertices, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated subdivision Loop refinement calculate subdivision surfaces smooth position
		if(markSH.get_subset_index(vrt) != -1 && linearManifoldSH.get_subset_index(vrt) == -1)
		{
		//	perform loop subdivision on even manifold vertices
		//	first get neighbored manifold vertices
			for(Grid::AssociatedEdgeIterator eIter = mg.associated_edges_begin(vrt); eIter != mg.associated_edges_end(vrt); ++eIter)
			{
				Edge* e = *eIter;

			//	Only consider associated edges, which are marked as manifold edges
				if(markSH.get_subset_index(e) != -1)
				{
				//	Exclude ghost and horizontal slave neighbor vertices from contributing to centroid
					#ifdef UG_PARALLEL
						if(dgm.is_ghost(e))
							continue;

						if(dgm.contains_status(e, ES_H_SLAVE))
							continue;
					#endif

					VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
				}
			}

			number centerWgt 	= subdiv.ref_even_inner_center_weight(aaNumManifoldEdges[vrt]);
			number nbrWgt 		= subdiv.ref_even_inner_nbr_weight(aaNumManifoldEdges[vrt]);

		//	Exclude horizontal slaves of the currently smoothed vertex to avoid multiple contributions to centroid
			#ifdef UG_PARALLEL
				if(dgm.is_ghost(vrt))
					continue;

				if(dgm.contains_status(vrt, ES_H_SLAVE))
				{
					VecScaleAppend(aaSmoothBndPosEvenVrt[vrt], nbrWgt, p);
					continue;
				}
			#endif

			VecScaleAppend(aaSmoothBndPosEvenVrt[vrt], centerWgt, aaPos[vrt], nbrWgt, p);
		}
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

//	Calculate smooth position for ODD vertices
	for(EdgeIterator eIter = mg.begin<Edge>(mg.top_level()-1); eIter != mg.end<Edge>(mg.top_level()-1); ++eIter)
	{
		VecSet(p, 0);

		Edge* e = *eIter;

	//	Skip ghost edges
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(e))
				continue;
		#endif

	//	In case of marked manifold edges, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated subdivision Loop refinement calculate subdivision surfaces smooth position
		if(markSH.get_subset_index(e) != -1 && linearManifoldSH.get_subset_index(e) == -1)
		{
		//	perform loop subdivision on odd manifold vertices
		//	get the neighbored manifold triangles
			std::vector<Face*> associatedFaces;
			std::vector<Face*> associatedManifoldFaces;

			CollectAssociated(associatedFaces, mg, e);

			for(size_t i = 0; i < associatedFaces.size(); ++i)
			{

			//	Only consider associated faces, which are marked as manifold faces
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
						UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelLoopScheme: "
							"Non triangular faces included in grid: " << ElementDebugInfo(mg, associatedManifoldFaces[i]));
					}
				}

			//	Summate centroid of face adjacent vertices
				for(size_t i = 0; i < associatedManifoldFaces.size(); ++i)
				{
					VecAdd(p, p, aaPos[GetConnectedVertex(e, associatedManifoldFaces[i])]);
				}

			//	Exclude ghost and horizontal slaves of the parent edge vertices of the currently smoothed vertex
			//	to avoid multiple contributions to the centroid of the edge adjacent vertices
				#ifdef UG_PARALLEL
					if(dgm.is_ghost(e))
					{
						continue;
					}

					if(dgm.contains_status(e, ES_H_SLAVE))
					{
						VecScaleAppend(aaSmoothBndPosOddVrt[e], 0.125, p);
						continue;
					}
				#endif

				VecScaleAppend(aaSmoothBndPosOddVrt[e], 0.375, aaPos[e->vertex(0)], 0.375, aaPos[e->vertex(1)], 0.125, p);
			}
			else
				UG_THROW("ERROR in CalculateSmoothManifoldPosInParentLevelLoopScheme: numAssociatedManifoldFaces > 2.");
		}
	}

//	Manage vertex and edge attachment communication in parallel case -> COMMUNICATE aSmoothBndPosEvenVrt, aSmoothBndPosOddVrt
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aSmoothBndPosEvenVrt, PCL_RO_SUM);
		AttachmentAllReduce<Edge>(mg, aSmoothBndPosOddVrt, PCL_RO_SUM);
	#endif
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CalculateSmoothManifoldPosInTopLevelAveragingScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
														 MGSubsetHandler& linearManifoldSH,
														 TAPosition& aSmoothBndPos_tri,
														 TAPosition& aSmoothBndPos_quad)
{
	/*
	 * Scheme reference:
	 *
	 * J. Warren and H. Weimer, Subdivision Methods for Geometric Design: A Constructive Approach,
	 * Morgan Kaufmann Publishers Inc., San Francisco, CA, USA, 1st ed., 2001.
	 */

	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in CalculateSmoothManifoldPosInTopLevelAveragingScheme:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);
	Grid::VertexAttachmentAccessor<TAPosition> aaSmoothBndPos_tri(mg, aSmoothBndPos_tri);
	Grid::VertexAttachmentAccessor<TAPosition> aaSmoothBndPos_quad(mg, aSmoothBndPos_quad);

	#ifdef UG_PARALLEL
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
	#endif

//	Declare centroid coordinate vector
	typedef typename TAPosition::ValueType position_type;
	position_type p;
	VecSet(p, 0);

//	Loop all manifold faces of top_level
	for(FaceIterator fIter = mg.begin<Face>(mg.top_level()); fIter != mg.end<Face>(mg.top_level()); ++fIter)
	{
		Face* f = *fIter;

	//	Skip ghost volumes
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(f))
				continue;
		#endif

	//	In case of marked manifold faces, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated Averaging scheme calculate subdivision surfaces smooth position
		if(markSH.get_subset_index(f) != -1 && linearManifoldSH.get_subset_index(f) == -1)
		{
		//	TRIANGLE CASE
			if(f->reference_object_id() == ROID_TRIANGLE)
			{
			//	Iterate over all face vertices, calculate and apply local centroid masks
				for(size_t i = 0; i < f->num_vertices(); ++i)
				{
				//	Init
					Vertex* vrt = f->vertex(i);
					VecSet(p, 0);

				//	Summate coordinates of neighbor vertices to vrt inside face
					for(size_t j = 0; j < f->num_vertices(); ++j)
					{
						if(j != i)
						{
							VecAdd(p, p, aaPos[f->vertex(j)]);
						}
					}

				//	Smooth vertex position
					VecScaleAppend(aaSmoothBndPos_tri[vrt], 2.0/8, aaPos[vrt], 3.0/8, p);
				}
			}

		//	QUADRILATERAL CASE
			else if(f->reference_object_id() == ROID_QUADRILATERAL)
			{
			//	Iterate over all face vertices, calculate and apply local centroid masks
				for(size_t i = 0; i < f->num_vertices(); ++i)
				{
				//	Init
					Vertex* vrt = f->vertex(i);
					VecSet(p, 0);

				//	Summate coordinates of neighbor vertices to vrt inside face
					for(size_t j = 0; j < f->num_vertices(); ++j)
					{
						if(j != i)
						{
							VecAdd(p, p, aaPos[f->vertex(j)]);
						}
					}

				//	Smooth vertex position
					VecScaleAppend(aaSmoothBndPos_quad[vrt], 1.0/4, aaPos[vrt], 1.0/4, p);
				}
			}

		//	UNSUPPORTED MANIFOLD ELEMENT CASE
			else
				UG_THROW("ERROR in CalculateSmoothManifoldPosInTopLevelAveragingScheme: Non triangular/quadrilateral faces included in grid.");
		}
	}

//	Manage vertex attachment communication in parallel case -> COMMUNICATE aaSmoothBndPos
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aSmoothBndPos_tri, PCL_RO_SUM);
		AttachmentAllReduce<Vertex>(mg, aSmoothBndPos_quad, PCL_RO_SUM);
	#endif
}


////////////////////////////////////////////////////////////////////////////////
void CalculateSmoothVolumePosInTopLevel(MultiGrid& mg, MGSubsetHandler& markSH,
										APosition& aSmoothVolPos_toc,
										APosition& aSmoothVolPos_prism,
										APosition& aSmoothVolPos_hex)
{
	/*
	 * Scheme references:
	 *
	 * S. Schaefer, J. Hakenberg, and J. Warren, Smooth subdivision of tetrahedral meshes,
	 * Proceedings of the 2004 Eurographics/ACM Symposium on Geometry Processing.
	 *
	 * J. Hakenberg, Smooth Subdivision for Mixed Volumetric Meshes, thesis, 2004
	 */

	#ifdef UG_PARALLEL
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
	#endif

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_toc(mg, aSmoothVolPos_toc);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_prism(mg, aSmoothVolPos_prism);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_hex(mg, aSmoothVolPos_hex);

//	Declare volume centroid coordinate vector
	typedef APosition::ValueType pos_type;
	pos_type p;

//	Loop all volumes of top_level
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

	//	Skip ghost volumes
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(vol))
				continue;
		#endif

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
				VecScaleAppend(aaSmoothVolPos_toc[vrt], -1.0/16, aaPos[vrt], 17.0/48, p);
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
						UG_THROW("ERROR in CalculateSmoothVolumePosInTopLevel: identified opposing vertex actually not included in current volume.");
					}

					if(j != i && j != (size_t)GetVertexIndex(vol, oppVrt))
					{
						VecAdd(p, p, aaPos[vol->vertex(j)]);
					}
				}

			//	Smooth vertex position
				VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 1.0/12, p, 7.0/24, aaPos[oppVrt]);
			}

		//	PRISM CASE
			else if(vol->reference_object_id() == ROID_PRISM)
			{
			//	Get cell-adjacent vertex
				Vertex* oppVrt;

				if(i == 0)
					oppVrt = vol->vertex(3);
				else if(i == 1)
					oppVrt = vol->vertex(4);
				else if(i == 2)
					oppVrt = vol->vertex(5);
				else if(i == 3)
					oppVrt = vol->vertex(0);
				else if(i == 4)
					oppVrt = vol->vertex(1);
				else
					oppVrt = vol->vertex(2);

			//	Summate coordinates of DIRECT neighbor vertices to vrt inside octahedron
				for(size_t j = 0; j < vol->num_vertices(); ++j)
				{
					if(j != i && j != (size_t)GetVertexIndex(vol, oppVrt))
					{
						VecAdd(p, p, aaPos[vol->vertex(j)]);
					}
				}

			//	Smooth vertex position
				VecScaleAppend(aaSmoothVolPos_prism[vrt], 1.0/8, aaPos[vrt], 3.0/16, p, 1.0/8, aaPos[oppVrt]);
			}

		//	HEXAHEDRON CASE
			else if(vol->reference_object_id() == ROID_HEXAHEDRON)
			{
			//	Summate coordinates of neighbor vertices to vrt inside hexahedron
				for(size_t j = 0; j < vol->num_vertices(); ++j)
				{
					if(j != i)
					{
						VecAdd(p, p, aaPos[vol->vertex(j)]);
					}
				}

			//	Smooth vertex position
				VecScaleAppend(aaSmoothVolPos_hex[vrt], 1.0/8, aaPos[vrt], 1.0/8, p);
			}

		//	UNSUPPORTED VOLUME ELEMENT CASE
			else
			{
				UG_THROW("ERROR in CalculateSmoothVolumePosInTopLevel: Volume type not supported for subdivision volumes refinement.");
			}
		}
	}


//	Manage vertex attachment communication in parallel case -> COMMUNICATE aSmoothVolPos
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aSmoothVolPos_toc, PCL_RO_SUM);
		AttachmentAllReduce<Vertex>(mg, aSmoothVolPos_prism, PCL_RO_SUM);
		AttachmentAllReduce<Vertex>(mg, aSmoothVolPos_hex, PCL_RO_SUM);
	#endif
}


////////////////////////////////////////////////////////////////////////////////
void CalculateConstrainedSmoothVolumePosInTopLevel(MultiGrid& mg, MGSubsetHandler& markSH,
												   APosition& aSmoothVolPos_toc)
{
	#ifdef UG_PARALLEL
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
	#endif

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_toc(mg, aSmoothVolPos_toc);

//	Declare volume centroid coordinate vector
	typedef APosition::ValueType pos_type;
	pos_type p;

//	boundary neighbor counter
	size_t bndNbrCnt;

//	Loop all volumes of top_level
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

	//	Skip ghost volumes
		#ifdef UG_PARALLEL
			if(dgm.is_ghost(vol))
				continue;
		#endif

	//	Iterate over all volume vertices, calculate and apply local centroid masks
		for(size_t i = 0; i < vol->num_vertices(); ++i)
		{
		//	Init
			Vertex* vrt = vol->vertex(i);
			VecSet(p, 0);
			bndNbrCnt = 0;

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
					//	Only consider non-manifold neighbors
						if(markSH.get_subset_index(vol->vertex(j)) == -1)
						{
							VecAdd(p, p, aaPos[vol->vertex(j)]);
						}
						else
						{
							bndNbrCnt++;
						}
					}
				}

			//	Smooth vertex position
				//VecScaleAppend(aaSmoothVolPos_toc[vrt], -1.0/16 + bndNbrCnt*17.0/48, aaPos[vrt], 17.0/48, p);

				if(bndNbrCnt == 3)
					VecScaleAppend(aaSmoothVolPos_toc[vrt], 1.0, aaPos[vrt]);
				else if(bndNbrCnt == 0)
					VecScaleAppend(aaSmoothVolPos_toc[vrt], -1.0/16, aaPos[vrt], 17.0/48, p);
				else if(bndNbrCnt == 1)
					VecScaleAppend(aaSmoothVolPos_toc[vrt], -1.0/16, aaPos[vrt], 17.0/48 + 17.0/(48*2), p);
				else
					VecScaleAppend(aaSmoothVolPos_toc[vrt], -1.0/16, aaPos[vrt], 51.0/48, p);
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
						UG_THROW("ERROR in CalculateConstrainedSmoothVolumePosInTopLevel: identified opposing vertex actually not included in current volume.");
					}

					if(j != i && j != (size_t)GetVertexIndex(vol, oppVrt))
					{
					//	Only consider non-manifold direct neighbors
						if(markSH.get_subset_index(vol->vertex(j)) == -1)
						{
							VecAdd(p, p, aaPos[vol->vertex(j)]);
						}
						else
						{
							bndNbrCnt++;
						}
					}
				}

			//	Smooth vertex position

			//	if opposing vertex is a non-manifold vertex
				if(markSH.get_subset_index(oppVrt) == -1)
				{
					//VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8 + bndNbrCnt*1.0/12, aaPos[vrt], 1.0/12, p, 7.0/24, aaPos[oppVrt]);

					if(bndNbrCnt == 4)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 15.0/24, aaPos[oppVrt]);
					else if(bndNbrCnt == 0)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 1.0/12, p, 7.0/24, aaPos[oppVrt]);
					else if(bndNbrCnt == 1)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 1.0/12 + 1.0/(12*3), p, 7.0/24, aaPos[oppVrt]);
					else if(bndNbrCnt == 2)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 2.0/12, p, 7.0/24, aaPos[oppVrt]);
					else
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 4.0/12, p, 7.0/24, aaPos[oppVrt]);
				}
				else
				{
					//VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8 + bndNbrCnt*1.0/12 + 7.0/24, aaPos[vrt], 1.0/12, p);

					if(bndNbrCnt == 4)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 1.0, aaPos[vrt]);
					else if(bndNbrCnt == 0)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 1.0/12 + 7.0/(24*4), p);
					else if(bndNbrCnt == 1)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 1.0/12 + 1.0/(12*3) + 7.0/(24*3), p);
					else if(bndNbrCnt == 2)
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 2.0/12 + 7.0/(24*2), p);
					else
						VecScaleAppend(aaSmoothVolPos_toc[vrt], 3.0/8, aaPos[vrt], 4.0/12 + 7.0/24, p);
				}
			}

		//	UNSUPPORTED VOLUME ELEMENT CASE
			else
			{
				UG_THROW("ERROR in CalculateConstrainedSmoothVolumePosInTopLevel: Volume type not supported for subdivision volumes refinement, only tetrahedra and octahedra.");
			}
		}
	}

//	Manage vertex attachment communication in parallel case -> COMMUNICATE aSmoothVolPos
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aSmoothVolPos_toc, PCL_RO_SUM);
	#endif
}


////////////////////////////////////////////////////////////////////////////////
void CalculateNumElemsVertexAttachmentInTopLevel(MultiGrid& mg, AInt& aNumElems_toc, AInt& aNumElems_prism, AInt& aNumElems_hex)
{
//	Define attachment accessor
	Grid::VertexAttachmentAccessor<AInt> aaNumElems_toc(mg, aNumElems_toc);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems_prism(mg, aNumElems_prism);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems_hex(mg, aNumElems_hex);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, AInt> comPolCopyNumElems_toc(mg, aNumElems_toc);
		ComPol_CopyAttachment<VertexLayout, AInt> comPolCopyNumElems_prism(mg, aNumElems_prism);
		ComPol_CopyAttachment<VertexLayout, AInt> comPolCopyNumElems_hex(mg, aNumElems_hex);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif

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
			if(vol->reference_object_id() == ROID_TETRAHEDRON || vol->reference_object_id() == ROID_OCTAHEDRON)
				++aaNumElems_toc[vol->vertex(i)];
			else if(vol->reference_object_id() == ROID_PRISM)
				++aaNumElems_prism[vol->vertex(i)];
			else if(vol->reference_object_id() == ROID_HEXAHEDRON)
				++aaNumElems_hex[vol->vertex(i)];
			else
				UG_THROW("ERROR in CalculateNumElemsVertexAttachmentInTopLevel: not supported element type included in grid.");
		}
	}

//	Manage vertex attachment communication in parallel case -> COMMUNICATE aNumElems
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aNumElems_toc, PCL_RO_SUM);
		AttachmentAllReduce<Vertex>(mg, aNumElems_prism, PCL_RO_SUM);
		AttachmentAllReduce<Vertex>(mg, aNumElems_hex, PCL_RO_SUM);

	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyNumElems_toc);
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyNumElems_prism);
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyNumElems_hex);
		com.communicate();
	#endif
}


////////////////////////////////////////////////////////////////////////////////
void CalculateNumManifoldEdgesVertexAttachmentInParentLevel(MultiGrid& mg, MGSubsetHandler& markSH,
															AInt& aNumManifoldEdges)
{
//	Define attachment accessor
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldEdges(mg, aNumManifoldEdges);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, AInt> comPolCopyNumManifoldEdges(mg, aNumManifoldEdges);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif

//	Catch use of CalculateNumManifoldEdgesVertexAttachmentInParentLevel on MultiGrids with just one level
	if(mg.num_levels() == 1)
		UG_THROW("CalculateNumManifoldEdgesVertexAttachmentInParentLevel: method may not be used in base level 0.");

//	Loop all edges of parent level and calculate number of associated manifold edges of each vertex
	for(EdgeIterator eIter = mg.begin<Edge>(mg.top_level()-1); eIter != mg.end<Edge>(mg.top_level()-1); ++eIter)
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

//	Manage vertex attachment communication in parallel case -> COMMUNICATE aNumManifoldEdges
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aNumManifoldEdges, PCL_RO_SUM);

	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyNumManifoldEdges);
		com.communicate();
	#endif
}


////////////////////////////////////////////////////////////////////////////////
void CalculateNumManifoldFacesVertexAttachmentInTopLevel(MultiGrid& mg, MGSubsetHandler& markSH, AInt& aNumManifoldFaces_tri, AInt& aNumManifoldFaces_quad)
{
//	Define attachment accessor
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldFaces_tri(mg, aNumManifoldFaces_tri);
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldFaces_quad(mg, aNumManifoldFaces_quad);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, AInt> comPolCopyNumManifoldFaces_tri(mg, aNumManifoldFaces_tri);
		ComPol_CopyAttachment<VertexLayout, AInt> comPolCopyNumManifoldFaces_quad(mg, aNumManifoldFaces_quad);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif

//	Loop all manifold faces of top level and calculate number of faces each vertex is contained by
	for(FaceIterator fIter = mg.begin<Face>(mg.top_level()); fIter != mg.end<Face>(mg.top_level()); ++fIter)
	{
		Face* f = *fIter;

	//	Only consider boundary manifold faces
		if(markSH.get_subset_index(f) != -1)
		{
		//	Skip ghosts
			#ifdef UG_PARALLEL
				if(dgm.is_ghost(f))
					continue;
			#endif

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				if(f->reference_object_id() == ROID_TRIANGLE)
					++aaNumManifoldFaces_tri[f->vertex(i)];
				else if(f->reference_object_id() == ROID_QUADRILATERAL)
					++aaNumManifoldFaces_quad[f->vertex(i)];
				else
					UG_THROW("ERROR in CalculateNumManifoldFacesVertexAttachmentInTopLevel: Non triangular/quadrilateral faces included in grid.");
			}
		}
	}

//	Manage vertex attachment communication in parallel case -> COMMUNICATE aNumElems
	#ifdef UG_PARALLEL
	//	Reduce add operations:
	//	sum up h_slaves into h_masters

	//	Copy operations:
	//	copy h_masters to h_slaves for consistency
		AttachmentAllReduce<Vertex>(mg, aNumManifoldFaces_tri, PCL_RO_SUM);
		AttachmentAllReduce<Vertex>(mg, aNumManifoldFaces_quad, PCL_RO_SUM);

	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyNumManifoldFaces_tri);
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyNumManifoldFaces_quad);
		com.communicate();
	#endif
}


////////////////////////////////////////////////////////////////////////////////
void InitLinearManifoldSubsetHandler(MultiGrid& mg, MGSubsetHandler& sh,
											   MGSubsetHandler& linearManifoldSH,
											   const char* linearManifoldSubsets)
{
	if(linearManifoldSubsets[0] == '\0')
		return;

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in InitLinearManifoldSubsetHandler: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}

//	tokenize user input
	std::vector<std::string> linearManifoldSubsetsString = TokenizeString(linearManifoldSubsets);

//	remove white space
	for(size_t i = 0; i < linearManifoldSubsetsString.size(); ++i)
	{
		RemoveWhitespaceFromString(linearManifoldSubsetsString[i]);
	}

//	if no subset passed, clear subsets
	if(linearManifoldSubsetsString.size() == 1 && linearManifoldSubsetsString[0].empty())
	{
		linearManifoldSubsetsString.clear();
	}

//	if subsets passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < linearManifoldSubsetsString.size(); ++i)
	{
		if(linearManifoldSubsetsString.empty())
		{
			UG_THROW(	"ERROR in InitLinearManifoldSubsetHandler: "
						"linear boundary manifold subsets string passed lacks a "
						"subset specification at position "<<i<<"(of "
						<<linearManifoldSubsetsString.size()-1<<")");
		}
	}

// 	assign all user specified vertices to linear boundary manifold SubsetHandler
	for(size_t i = 0; i < linearManifoldSubsetsString.size(); ++i)
	{
		int j = sh.get_subset_index(linearManifoldSubsetsString[i].c_str());
		UG_COND_THROW(j < 0, "Linear manifold subset named '"
			<< linearManifoldSubsetsString[i] << "' could not be identified.");

		for(VertexIterator vrtIter = sh.begin<Vertex>(j, mg.top_level()); vrtIter != sh.end<Vertex>(j, mg.top_level()); ++vrtIter)
		{
			Vertex* vrt = *vrtIter;
			linearManifoldSH.assign_subset(vrt, 0);
		}

		for(VertexIterator vrtIter = sh.begin<Vertex>(j, mg.top_level()-1); vrtIter != sh.end<Vertex>(j, mg.top_level()-1); ++vrtIter)
		{
			Vertex* vrt = *vrtIter;
			linearManifoldSH.assign_subset(vrt, 0);
		}
	}

// 	assign all user specified edges to linear boundary manifold SubsetHandler
	for(size_t i = 0; i < linearManifoldSubsetsString.size(); ++i)
	{
		int j = sh.get_subset_index(linearManifoldSubsetsString[i].c_str());

		for(EdgeIterator eIter = sh.begin<Edge>(j, mg.top_level()); eIter != sh.end<Edge>(j, mg.top_level()); ++eIter)
		{
			Edge* e = *eIter;
			linearManifoldSH.assign_subset(e, 0);
		}

		for(EdgeIterator eIter = sh.begin<Edge>(j, mg.top_level()-1); eIter != sh.end<Edge>(j, mg.top_level()-1); ++eIter)
		{
			Edge* e = *eIter;
			linearManifoldSH.assign_subset(e, 0);
		}
	}

// 	assign all user specified faces to linear boundary manifold SubsetHandler
	for(size_t i = 0; i < linearManifoldSubsetsString.size(); ++i)
	{
		int j = sh.get_subset_index(linearManifoldSubsetsString[i].c_str());

		for(FaceIterator fIter = sh.begin<Face>(j, mg.top_level()); fIter != sh.end<Face>(j, mg.top_level()); ++fIter)
		{
			Face* f = *fIter;
			linearManifoldSH.assign_subset(f, 0);
		}

		for(FaceIterator fIter = sh.begin<Face>(j, mg.top_level()-1); fIter != sh.end<Face>(j, mg.top_level()-1); ++fIter)
		{
			Face* f = *fIter;
			linearManifoldSH.assign_subset(f, 0);
		}
	}

//	Debug log
//	UG_LOG("InitLinearManifoldSubsetHandler:" << std::endl);
//	UG_LOG(">> Applying linear subdivision on the following boundary manifold subsets:" << std::endl);
//
//	for(size_t i = 0; i < linearManifoldSubsetsString.size(); ++i)
//	{
//		UG_LOG("Subset # " << sh.get_subset_index(linearManifoldSubsetsString[i].c_str()) << ": " << linearManifoldSubsetsString[i] << std::endl);
//	}
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void ApplySmoothManifoldPosToTopLevelLoopScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
												MGSubsetHandler& linearManifoldSH)
{
	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in ApplySmoothManifoldPosToTopLevelLoopScheme:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in ApplySmoothManifoldPosToTopLevelLoopScheme: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}


/*****************************************
 *
 *	(1) SETUP
 *
 *****************************************/

//	Position attachment value type
	typedef typename TAPosition::ValueType position_type;

//	Vertex attachments for associated number of manifold edges and smooth position
//	(distinguish between volume and boundary smooth vertex positions
//	 and in case of boundary between EVEN and ODD smooth vertex positions)
	AInt aNumManifoldEdges;
	TAPosition aSmoothBndPosEvenVrt;
	TAPosition aSmoothBndPosOddVrt;

//	attach previously declared vertex attachments with initial value 0
	mg.attach_to_vertices_dv(aNumManifoldEdges, 0);
	mg.attach_to_vertices_dv(aSmoothBndPosEvenVrt, position_type(0));
	mg.attach_to_edges_dv(aSmoothBndPosOddVrt, position_type(0));

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);
	Grid::VertexAttachmentAccessor<TAPosition> aaSmoothBndPosEvenVrt(mg, aSmoothBndPosEvenVrt);
	Grid::EdgeAttachmentAccessor<TAPosition> aaSmoothBndPosOddVrt(mg, aSmoothBndPosOddVrt);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment aPosition
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, TAPosition> comPolCopyAPosition(mg, aPos);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif


/*****************************************
 *
 *	(2) DETERMINE aNumManifoldEdges
 *
 *****************************************/

	CalculateNumManifoldEdgesVertexAttachmentInParentLevel(mg, markSH, aNumManifoldEdges);


/*****************************************
 *
 *	(3) CALCULATE aSmoothBndPosEvenVrt,
 *				  aSmoothBndPosOddVrt
 *
 *****************************************/

//	Calculate aSmoothBndPosEvenVrt, aSmoothBndPosOddVrt
	CalculateSmoothManifoldPosInParentLevelLoopScheme(mg, aPos, markSH, linearManifoldSH,
												aSmoothBndPosEvenVrt, aSmoothBndPosOddVrt, aNumManifoldEdges);


/*****************************************
 *
 *	(4) APPLY
 *
 *****************************************/

//	Loop all vertices of top_level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

	//	Catch vertices without parent
		if(mg.get_parent(vrt) == NULL)
			continue;

	//	In case of marked manifold vertices, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated Loop scheme refinement apply subdivision surfaces smoothing, else linear refinement
		if(markSH.get_subset_index(vrt) != -1 && linearManifoldSH.get_subset_index(vrt) == -1)
		{
		//	EVEN VERTEX
			if(mg.get_parent(vrt)->reference_object_id() == ROID_VERTEX)
			{
			//	Get parent vertex
				Vertex* parentVrt = static_cast<Vertex*>(mg.get_parent(vrt));

				aaPos[vrt] = aaSmoothBndPosEvenVrt[parentVrt];
			}

		//	ODD VERTEX
			else if(mg.get_parent(vrt)->reference_object_id() == ROID_EDGE)
			{
			//	Get parent edge
				Edge* parentEdge = static_cast<Edge*>(mg.get_parent(vrt));

				aaPos[vrt] =  aaSmoothBndPosOddVrt[parentEdge];
			}
		}
	}


/*****************************************
 *
 *	(5) COMMUNICATE VERTICALLY
 *	    AFTER SUBDIVISION SURFACES
 *
 *****************************************/

//	Communicate aPosition in parallel case
	#ifdef UG_PARALLEL
	//	copy ghosts = VMASTER to v_slaves
		com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, comPolCopyAPosition);
		com.communicate();
	#endif


/*****************************************
 *
 *	(6) CLEAN UP
 *
 *****************************************/

//	detach vertex attachments
	mg.detach_from_vertices(aNumManifoldEdges);
	mg.detach_from_vertices(aSmoothBndPosEvenVrt);
	mg.detach_from_edges(aSmoothBndPosOddVrt);
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void ApplySmoothManifoldPosToTopLevelButterflyScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
												MGSubsetHandler& linearManifoldSH)
{
	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in ApplySmoothManifoldPosToTopLevelButterflyScheme:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in ApplySmoothManifoldPosToTopLevelButterflyScheme: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}


/*****************************************
 *
 *	(1) SETUP
 *
 *****************************************/

//	Position attachment value type
	typedef typename TAPosition::ValueType position_type;

//	Vertex attachments for associated number of manifold edges and smooth position
//	(distinguish between volume and boundary smooth vertex positions
//	 and in case of boundary between EVEN and ODD smooth vertex positions)
	AInt aNumManifoldEdges;
	TAPosition aSmoothBndPosOddVrt;

//	attach previously declared vertex attachments with initial value 0
	mg.attach_to_vertices_dv(aNumManifoldEdges, 0);
	mg.attach_to_edges_dv(aSmoothBndPosOddVrt, position_type(0));

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);
	Grid::EdgeAttachmentAccessor<TAPosition> aaSmoothBndPosOddVrt(mg, aSmoothBndPosOddVrt);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment aPosition
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, TAPosition> comPolCopyAPosition(mg, aPos);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif


/*****************************************
 *
 *	(2) DETERMINE aNumManifoldEdges
 *
 *****************************************/

	CalculateNumManifoldEdgesVertexAttachmentInParentLevel(mg, markSH, aNumManifoldEdges);


/*****************************************
 *
 *	(3) CALCULATE aSmoothBndPosEvenVrt,
 *				  aSmoothBndPosOddVrt
 *
 *****************************************/

//	Calculate aSmoothBndPosOddVrt
	CalculateSmoothManifoldPosInParentLevelButterflyScheme(mg, aPos, markSH, linearManifoldSH, aSmoothBndPosOddVrt, aNumManifoldEdges);


/*****************************************
 *
 *	(4) APPLY
 *
 *****************************************/

//	Loop all vertices of top_level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

	//	Catch vertices without parent
		if(mg.get_parent(vrt) == NULL)
			continue;

	//	In case of marked manifold vertices, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated Loop scheme refinement apply subdivision surfaces smoothing, else linear refinement
		if(markSH.get_subset_index(vrt) != -1 && linearManifoldSH.get_subset_index(vrt) == -1)
		{
		//	ODD VERTEX
			if(mg.get_parent(vrt)->reference_object_id() == ROID_EDGE)
			{
			//	Get parent edge
				Edge* parentEdge = static_cast<Edge*>(mg.get_parent(vrt));

				aaPos[vrt] =  aaSmoothBndPosOddVrt[parentEdge];
			}
		}
	}


/*****************************************
 *
 *	(5) COMMUNICATE VERTICALLY
 *	    AFTER SUBDIVISION SURFACES
 *
 *****************************************/

//	Communicate aPosition in parallel case
	#ifdef UG_PARALLEL
	//	copy ghosts = VMASTER to v_slaves
		com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, comPolCopyAPosition);
		com.communicate();
	#endif


/*****************************************
 *
 *	(6) CLEAN UP
 *
 *****************************************/

//	detach vertex attachments
	mg.detach_from_vertices(aNumManifoldEdges);
	mg.detach_from_edges(aSmoothBndPosOddVrt);
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void ApplySmoothManifoldPosToTopLevelAveragingScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
													 MGSubsetHandler& linearManifoldSH)
{
/*****************************************
 *
 *	(1) SETUP
 *
 *****************************************/

	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in ApplySmoothManifoldPosToTopLevelAveragingScheme:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Position attachment value type
	typedef typename TAPosition::ValueType position_type;

//	Vertex attachments for associated number of manifold faces and smooth position
	AInt aNumManifoldFaces_tri;
	AInt aNumManifoldFaces_quad;
	TAPosition aSmoothBndPos_tri;
	TAPosition aSmoothBndPos_quad;

//	attach previously declared vertex attachments with initial value 0
	mg.attach_to_vertices_dv(aNumManifoldFaces_tri, 0);
	mg.attach_to_vertices_dv(aNumManifoldFaces_quad, 0);
	mg.attach_to_vertices_dv(aSmoothBndPos_tri, position_type(0));
	mg.attach_to_vertices_dv(aSmoothBndPos_quad, position_type(0));

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aPos);
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldFaces_tri(mg, aNumManifoldFaces_tri);
	Grid::VertexAttachmentAccessor<AInt> aaNumManifoldFaces_quad(mg, aNumManifoldFaces_quad);
	Grid::VertexAttachmentAccessor<TAPosition> aaSmoothBndPos_tri(mg, aSmoothBndPos_tri);
	Grid::VertexAttachmentAccessor<TAPosition> aaSmoothBndPos_quad(mg, aSmoothBndPos_quad);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment aPosition
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, TAPosition> comPolCopyAPosition(mg, aPos);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif


/*****************************************
 *
 *	(2) DETERMINE aNumManifoldEdges
 *
 *****************************************/

	CalculateNumManifoldFacesVertexAttachmentInTopLevel(mg, markSH, aNumManifoldFaces_tri, aNumManifoldFaces_quad);


/*****************************************
 *
 *	(3) CALCULATE
 *
 *****************************************/

//	Calculate aSmoothBndPosEvenVrt, aSmoothBndPosOddVrt
	CalculateSmoothManifoldPosInTopLevelAveragingScheme(mg, aPos, markSH, linearManifoldSH, aSmoothBndPos_tri, aSmoothBndPos_quad);


/*****************************************
 *
 *	(4) APPLY
 *
 *****************************************/

//	Move manifold vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

	//	In case of marked manifold vertices, which do not belong to the user-specified linear boundary manifold subsets,
	//	and activated averaging scheme apply subdivision surfaces smoothing, else linear refinement
		if(markSH.get_subset_index(vrt) != -1 && linearManifoldSH.get_subset_index(vrt) == -1)
		{
			if(aaNumManifoldFaces_tri[vrt] == 0 && aaNumManifoldFaces_quad[vrt] == 0)
				UG_THROW("ERROR in ApplySmoothManifoldPosToTopLevelAveragingScheme: grid contains manifold vertex not contained in any manifold face.");

		//	Scale smooth vertex position by the number of associated volume elements (SubdivisionVolumes smoothing)
			VecScale(aaSmoothBndPos_tri[vrt],  aaSmoothBndPos_tri[vrt], 1.0/6.0/(aaNumManifoldFaces_tri[vrt]/6.0 + aaNumManifoldFaces_quad[vrt]/4.0));
			VecScale(aaSmoothBndPos_quad[vrt],  aaSmoothBndPos_quad[vrt], 1.0/4.0/(aaNumManifoldFaces_tri[vrt]/6.0 + aaNumManifoldFaces_quad[vrt]/4.0));
			VecScaleAdd(aaPos[vrt], 1.0, aaSmoothBndPos_tri[vrt], 1.0, aaSmoothBndPos_quad[vrt]);
		}
	}


/*****************************************
 *
 *	(5) COMMUNICATE VERTICALLY
 *	    AFTER SUBDIVISION SURFACES
 *
 *****************************************/

//	Communicate aPosition in parallel case
	#ifdef UG_PARALLEL
	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyAPosition);
		com.communicate();
	#endif


/*****************************************
 *
 *	(6) CLEAN UP
 *
 *****************************************/

//	detach vertex attachments
	mg.detach_from_vertices(aNumManifoldFaces_tri);
	mg.detach_from_vertices(aNumManifoldFaces_quad);
	mg.detach_from_vertices(aSmoothBndPos_tri);
	mg.detach_from_vertices(aSmoothBndPos_quad);
}


////////////////////////////////////////////////////////////////////////////////
void ApplySmoothVolumePosToTopLevel(MultiGrid& mg, MGSubsetHandler& markSH,
									MGSubsetHandler& linearManifoldSH, bool bConstrained)
{
/*****************************************
 *
 *	(1) SETUP
 *
 *****************************************/

//	Vertex attachments for associated number of elements and smooth position
	AInt aNumElems_toc;
	AInt aNumElems_prism;
	AInt aNumElems_hex;
	APosition aSmoothVolPos_toc;
	APosition aSmoothVolPos_prism;
	APosition aSmoothVolPos_hex;

//	attach previously declared vertex attachments with initial value 0
	mg.attach_to_vertices_dv(aNumElems_toc, 0);
	mg.attach_to_vertices_dv(aNumElems_prism, 0);
	mg.attach_to_vertices_dv(aNumElems_hex, 0);
	mg.attach_to_vertices_dv(aSmoothVolPos_toc, vector3(0, 0, 0));
	mg.attach_to_vertices_dv(aSmoothVolPos_prism, vector3(0, 0, 0));
	mg.attach_to_vertices_dv(aSmoothVolPos_hex, vector3(0, 0, 0));

//	Define attachment accessors
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems_toc(mg, aNumElems_toc);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems_prism(mg, aNumElems_prism);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems_hex(mg, aNumElems_hex);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_toc(mg, aSmoothVolPos_toc);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_prism(mg, aSmoothVolPos_prism);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothVolPos_hex(mg, aSmoothVolPos_hex);

//	Manage vertex attachment communication in parallel case:
//	- Setup communication policy for the above attachment aPosition
//	- Setup interface communicator
//	- Setup distributed grid manager
//	- Setup grid layout map
	#ifdef UG_PARALLEL
	//	Attachment communication policies COPY
		ComPol_CopyAttachment<VertexLayout, AVector3> comPolCopyAPosition(mg, aPosition);

	//	Interface communicators and distributed domain manager
		pcl::InterfaceCommunicator<VertexLayout> com;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
	#endif


/*****************************************
 *
 *	(2) DETERMINE aNumElems
 *
 *****************************************/

	CalculateNumElemsVertexAttachmentInTopLevel(mg, aNumElems_toc, aNumElems_prism, aNumElems_hex);


/*****************************************
 *
 *	(3) CALCULATE
 *
 *****************************************/

//	Calculate aSmoothVolPos
	if(bConstrained)
		CalculateConstrainedSmoothVolumePosInTopLevel(mg, markSH, aSmoothVolPos_toc);
	else
		CalculateSmoothVolumePosInTopLevel(mg, markSH, aSmoothVolPos_toc, aSmoothVolPos_prism, aSmoothVolPos_hex);


/*****************************************
 *
 *	(4) APPLY
 *
 *****************************************/

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		if(aaNumElems_toc[vrt] == 0 &&  aaNumElems_prism[vrt] == 0 && aaNumElems_hex[vrt] == 0)
			UG_THROW("ERROR in ApplySmoothVolumePosToTopLevel: grid contains vertex not contained in any volume.");

		if(g_boundaryRefinementRule == SUBDIV_VOL)
		{
		//	Scale smooth vertex position by the number of associated volume elements (SubdivisionVolumes smoothing)
			VecScale(aaSmoothVolPos_toc[vrt],  aaSmoothVolPos_toc[vrt], 1.0/14.0/(aaNumElems_toc[vrt]/14.0 + aaNumElems_prism[vrt]/12.0 + aaNumElems_hex[vrt]/8.0));
			VecScale(aaSmoothVolPos_prism[vrt],  aaSmoothVolPos_prism[vrt], 1.0/12.0/(aaNumElems_toc[vrt]/14.0 + aaNumElems_prism[vrt]/12.0 + aaNumElems_hex[vrt]/8.0));
			VecScale(aaSmoothVolPos_hex[vrt],  aaSmoothVolPos_hex[vrt], 1.0/8.0/(aaNumElems_toc[vrt]/14.0 + aaNumElems_prism[vrt]/12.0 + aaNumElems_hex[vrt]/8.0));
			VecScaleAdd(aaPos[vrt], 1.0, aaSmoothVolPos_toc[vrt], 1.0, aaSmoothVolPos_prism[vrt], 1.0, aaSmoothVolPos_hex[vrt]);
		}
		else
		{
		//	Only in case of inner vertices
			if(markSH.get_subset_index(vrt) == -1)
			{
			//	Scale smooth vertex position by the number of associated volume elements
				VecScale(aaSmoothVolPos_toc[vrt],  aaSmoothVolPos_toc[vrt], 1.0/14.0/(aaNumElems_toc[vrt]/14.0 + aaNumElems_prism[vrt]/12.0 + aaNumElems_hex[vrt]/8.0));
				VecScale(aaSmoothVolPos_prism[vrt],  aaSmoothVolPos_prism[vrt], 1.0/12.0/(aaNumElems_toc[vrt]/14.0 + aaNumElems_prism[vrt]/12.0 + aaNumElems_hex[vrt]/8.0));
				VecScale(aaSmoothVolPos_hex[vrt],  aaSmoothVolPos_hex[vrt], 1.0/8.0/(aaNumElems_toc[vrt]/14.0 + aaNumElems_prism[vrt]/12.0 + aaNumElems_hex[vrt]/8.0));
				VecScaleAdd(aaPos[vrt], 1.0, aaSmoothVolPos_toc[vrt], 1.0, aaSmoothVolPos_prism[vrt], 1.0, aaSmoothVolPos_hex[vrt]);
			}
		}
	}


/*****************************************
 *
 *	(5) COMMUNICATE VERTICALLY
 *	    AFTER SUBDIVISION VOLUMES
 *
 *****************************************/

//	Communicate aPosition in parallel case
	#ifdef UG_PARALLEL
	//	copy v_slaves to ghosts = VMASTER
		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, comPolCopyAPosition);
		com.communicate();
	#endif


/*****************************************
 *
 *	(6) CLEAN UP
 *
 *****************************************/

//	detach vertex attachments
	mg.detach_from_vertices(aNumElems_toc);
	mg.detach_from_vertices(aNumElems_prism);
	mg.detach_from_vertices(aNumElems_hex);
	mg.detach_from_vertices(aSmoothVolPos_toc);
	mg.detach_from_vertices(aSmoothVolPos_prism);
	mg.detach_from_vertices(aSmoothVolPos_hex);
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void ApplySmoothSubdivisionSurfacesToTopLevel(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& sh,
											  MGSubsetHandler& markSH, MGSubsetHandler& linearManifoldSH)
{
/*****************************************
 *
 *	(1) SETUP
 *
 *****************************************/

	PROFILE_FUNC_GROUP("subdivision_volumes");

	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in ApplySmoothSubdivisionSurfacesToTopLevel:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in ApplySmoothSubdivisionSurfacesToTopLevel: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}


/*****************************************
 *
 *	(2) SUBDIVISION SURFACES
 *
 *****************************************/

	if(g_boundaryRefinementRule == SUBDIV_SURF_LOOP_SCHEME)
		ApplySmoothManifoldPosToTopLevelLoopScheme(mg, aPos, markSH, linearManifoldSH);
	else if(g_boundaryRefinementRule == SUBDIV_SURF_AVERAGING_SCHEME)
		ApplySmoothManifoldPosToTopLevelAveragingScheme(mg, aPos, markSH, linearManifoldSH);
	else if(g_boundaryRefinementRule == SUBDIV_SURF_BUTTERFLY_SCHEME)
		ApplySmoothManifoldPosToTopLevelButterflyScheme(mg, aPos, markSH, linearManifoldSH);
	else if(g_boundaryRefinementRule == SUBDIV_VOL){}
	else if(g_boundaryRefinementRule == LINEAR){}
	else
		UG_THROW("ERROR in ApplySmoothSubdivisionSurfacesToTopLevel: Unknown boundary refinement rule. Known rules are 'subdiv_surf_loop_scheme', 'subdiv_surf_averaging_scheme', 'subdiv_surf_butterfly_scheme' or 'linear'.");
}


////////////////////////////////////////////////////////////////////////////////
void ApplySmoothSubdivisionVolumesToTopLevel(MultiGrid& mg, MGSubsetHandler& sh, MGSubsetHandler& markSH,
											 MGSubsetHandler& linearManifoldSH, bool bConstrained)
{
/*****************************************
 *
 *	(1) SETUP
 *
 *****************************************/

	PROFILE_FUNC_GROUP("subdivision_volumes");

//	Ensure, that hybrid tet-/oct refinement is used as refinement rule for tetrahedrons
	if(tet_rules::GetRefinementRule() != tet_rules::HYBRID_TET_OCT)
		UG_THROW("ERROR in ApplySubdivisionVolumesToTopLevel: Set necessary refinement rule by SetTetRefinementRule('hybrid_tet_oct').");

//	Catch use of procedure for MultiGrids with just one level
	if(mg.num_levels() == 1)
	{
		UG_THROW("Error in ApplySmoothSubdivisionToTopLevel: "
				 "Procedure only to be used for MultiGrids with more than one level.");
	}


/*****************************************
 *
 *	(2) SUBDIVISION SURFACES
 *
 *****************************************/

	if(g_boundaryRefinementRule == SUBDIV_SURF_LOOP_SCHEME)
		ApplySmoothManifoldPosToTopLevelLoopScheme(mg, aPosition, markSH, linearManifoldSH);
	else if(g_boundaryRefinementRule == SUBDIV_SURF_AVERAGING_SCHEME)
		ApplySmoothManifoldPosToTopLevelAveragingScheme(mg, aPosition, markSH, linearManifoldSH);
	else if(g_boundaryRefinementRule == SUBDIV_SURF_BUTTERFLY_SCHEME)
		ApplySmoothManifoldPosToTopLevelButterflyScheme(mg, aPosition, markSH, linearManifoldSH);
	else if(g_boundaryRefinementRule == SUBDIV_VOL){}
	else if(g_boundaryRefinementRule == LINEAR){}
	else
		UG_THROW("ERROR in ApplySubdivisionVolumesToTopLevel: Unknown boundary refinement rule. Known rules are 'subdiv_surf_loop_scheme', 'subdiv_surf_averaging_scheme', 'subdiv_surf_butterfly_scheme', 'linear' or 'subdiv_vol'.");


/*****************************************
 *
 *	(3) SUBDIVISION VOLUMES
 *
 *****************************************/

	ApplySmoothVolumePosToTopLevel(mg, markSH, linearManifoldSH, bConstrained);
}


//////////////////////////////////////////////////////////////////////////////
//	Wrapper procedures
template <class TAPosition>
void ApplySmoothSubdivisionSurfacesToTopLevel(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& sh,
											  MGSubsetHandler& markSH, const char* linearManifoldSubsets)
{
	MGSubsetHandler linearManifoldSH(mg);
	InitLinearManifoldSubsetHandler(mg, sh, linearManifoldSH, linearManifoldSubsets);

	ApplySmoothSubdivisionSurfacesToTopLevel(mg, aPos, sh, markSH, linearManifoldSH);
}

void ApplySmoothSubdivisionVolumesToTopLevel(MultiGrid& mg, MGSubsetHandler& sh, MGSubsetHandler& markSH,
											 const char* linearManifoldSubsets)
{
	MGSubsetHandler linearManifoldSH(mg);
	InitLinearManifoldSubsetHandler(mg, sh, linearManifoldSH, linearManifoldSubsets);

	ApplySmoothSubdivisionVolumesToTopLevel(mg, sh, markSH, linearManifoldSH, false);
}

void ApplyConstrainedSmoothSubdivisionVolumesToTopLevel(MultiGrid& mg, MGSubsetHandler& sh, MGSubsetHandler& markSH,
														const char* linearManifoldSubsets)
{
	MGSubsetHandler linearManifoldSH(mg);
	InitLinearManifoldSubsetHandler(mg, sh, linearManifoldSH, linearManifoldSubsets);

	ApplySmoothSubdivisionVolumesToTopLevel(mg, sh, markSH, linearManifoldSH, true);
}


//////////////////////////////////////////////////////////////////////////////
//	Explicit instantiations
template void ApplySmoothSubdivisionSurfacesToTopLevel<APosition1>(MultiGrid& mg, APosition1& aPos, MGSubsetHandler& sh,
		  	  	  	  	  	  	  	  	  	  	  	  	  	  MGSubsetHandler& markSH, const char* linearManifoldSubsets);
template void ApplySmoothSubdivisionSurfacesToTopLevel<APosition2>(MultiGrid& mg, APosition2& aPos, MGSubsetHandler& sh,
		  	  	  	  	  	  	  	  	  	  	  	  	  	  MGSubsetHandler& markSH, const char* linearManifoldSubsets);
template void ApplySmoothSubdivisionSurfacesToTopLevel<APosition>(MultiGrid& mg, APosition& aPos, MGSubsetHandler& sh,
		  	  	  	  	  	  	  	  	  	  	  	  	  	  MGSubsetHandler& markSH, const char* linearManifoldSubsets);

template void ProjectHierarchyToSubdivisionLimit(MultiGrid& mg, APosition1& aPos);
template void ProjectHierarchyToSubdivisionLimit(MultiGrid& mg, APosition2& aPos);
template void ProjectHierarchyToSubdivisionLimit(MultiGrid& mg, APosition& aPos);



}//	end of namespace
