/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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
#include <stack>
#include <queue>
#include <map>
#include "subset_util.h"
#include "geom_obj_util/geom_obj_util.h"
#include "polychain_util.h"
#include "orientation_util.h"
#include "subset_color_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
int GetFirstFreeSubset(const ISubsetHandler& sh)
{
	for(int i = 0; i < sh.num_subsets(); ++i){
		if(!(sh.contains_vertices(i) || sh.contains_edges(i)
			 || sh.contains_faces(i) || sh.contains_volumes(i)))
		 {
			return i;
		 }
	}

	return sh.num_subsets();
}


////////////////////////////////////////////////////////////////////////
//	EraseEmptySubsets
///	Erases all subsets which do not contain any geometric objects
void EraseEmptySubsets(ISubsetHandler& sh)
{
	int si = 0;
	while(si < sh.num_subsets()){
		if(!(	sh.contains_vertices(si)
			||	sh.contains_edges(si)
			||	sh.contains_faces(si)
			||	sh.contains_volumes(si)))
		{
		//	the subset is empty
			sh.erase_subset(si);
		}
		else
			++si;
	}
}

////////////////////////////////////////////////////////////////////////
//	AssignGridToSubset
void AssignGridToSubset(Grid& g, ISubsetHandler& sh, int subsetInd)
{
	UG_ASSERT(&g == sh.grid(), "Specified subset-handler has to operate on the specified grid!");

	sh.assign_subset(g.begin<Vertex>(), g.end<Vertex>(), subsetInd);
	sh.assign_subset(g.begin<Edge>(), g.end<Edge>(), subsetInd);
	sh.assign_subset(g.begin<Face>(), g.end<Face>(), subsetInd);
	sh.assign_subset(g.begin<Volume>(), g.end<Volume>(), subsetInd);
}

////////////////////////////////////////////////////////////////////////
//	AssignSelectionToSubset
void AssignSelectionToSubset(ISelector& sel, ISubsetHandler& sh, int subsetInd)
{
	UG_ASSERT(sel.grid() == sh.grid(), "Specified selector and subset-handler "
									   "have to operate on the same grid!");

	GridObjectCollection selGoc = sel.get_grid_objects();

	for(size_t i = 0; i < selGoc.num_levels(); ++i){
		sh.assign_subset(selGoc.begin<Vertex>(i),
						 selGoc.end<Vertex>(i), subsetInd);
		sh.assign_subset(selGoc.begin<Edge>(i),
						 selGoc.end<Edge>(i), subsetInd);
		sh.assign_subset(selGoc.begin<Face>(i),
						 selGoc.end<Face>(i), subsetInd);
		sh.assign_subset(selGoc.begin<Volume>(i),
						 selGoc.end<Volume>(i), subsetInd);
	}
}

////////////////////////////////////////////////////////////////////////
//	AssignInterfaceEdgesToSubsets
void AssignFaceInterfaceEdgesToSubsets(Grid& grid, SubsetHandler& sh)
{
//TODO:	This algorithm can be improved.
//	In the moment it runs with arbitrarily O(numSubs^2*averageNumFacesPerSub)

//	get the first free edge-subset
	int edgeSubsetIndex = GetMaxSubsetIndex<Edge>(sh) + 1;

//	use this vector to collect adjacent faces
	vector<Edge*> vEdges;
	vector<Face*> vFaces;

//	iterate through all subsets
	for(int j = -1; j < sh.num_subsets(); ++j)
	{
	//	find interface-edges to the other subsets.
		for(int i = j + 1; i < sh.num_subsets(); ++i)
		{
		//	iterate through all faces of subset i and check
		//	whether some edges are adjacent to subset j.
			for(FaceIterator iter = sh.begin<Face>(i);
				iter != sh.end<Face>(i); ++iter)
			{
				Face* f = *iter;
				CollectEdges(vEdges, grid, f);
				for(uint k = 0; k < vEdges.size(); ++k)
				{
					Edge* e = vEdges[k];

					if(sh.get_subset_index(e) == -1)
					{
						CollectFaces(vFaces, grid, e);

						if(vFaces.size() == 1)
							sh.assign_subset(e, edgeSubsetIndex);
						else
						{
							for(uint l = 0; l < vFaces.size(); ++l)
							{
								if(sh.get_subset_index(vFaces[l]) == j)
								{
									sh.assign_subset(e, edgeSubsetIndex);
									break;
								}
							}
						}
					}
				}
			}

		//	if we assigned edges, we will proceed with the next subset
			if(edgeSubsetIndex < (int)sh.num_subsets())
				if(sh.num_elements<Edge>(edgeSubsetIndex) > 0)
					edgeSubsetIndex++;
		}
	}
}

////////////////////////////////////////////////////////////////////////
// AssignVolumeInterfaceFacesToSubsets
void AssignVolumeInterfaceFacesToSubsets(Grid& grid, SubsetHandler& sh)
{
//TODO:	This algorithm can be improved.
//	In the moment it runs with arbitrarily O(numSubs^2*averageNumVolumesPerSub)

//	get the first free face-subset
	int faceSubsetIndex = GetMaxSubsetIndex<Face>(sh) + 1;

//	use this vector to collect adjacent faces and volumes
	vector<Face*> vFaces;
	vector<Volume*> vVolumes;

//	iterate through all subsets
	for(int j = -1; j < (int)sh.num_subsets(); ++j)
	{
	//	find interface-edges to the other subsets.
		for(int i = j + 1; i < (int)sh.num_subsets(); ++i)
		{
		//	iterate through all volumes of subset i and check
		//	whether some faces are adjacent to subset j.
			for(VolumeIterator iter = sh.begin<Volume>(i);
				iter != sh.end<Volume>(i); ++iter)
			{
				Volume* v = *iter;
				CollectFaces(vFaces, grid, v);
				for(uint k = 0; k < vFaces.size(); ++k)
				{
					Face* f = vFaces[k];

					if(sh.get_subset_index(f) == -1)
					{
						CollectVolumes(vVolumes, grid, f);

						if(vVolumes.size() == 1)
							sh.assign_subset(f, faceSubsetIndex);
						else
						{
							for(uint l = 0; l < vVolumes.size(); ++l)
							{
								if(sh.get_subset_index(vVolumes[l]) == j)
								{
									sh.assign_subset(f, faceSubsetIndex);
									break;
								}
							}
						}
					}
				}
			}

		//	if we assigned edges, we will proceed with the next subset
			if(faceSubsetIndex < (int)sh.num_subsets())
				if(sh.num_elements<Face>(faceSubsetIndex) > 0)
					faceSubsetIndex++;
		}
	}
}


////////////////////////////////////////////////////////////////////////
/**	Helper method that copies subset-indices based on the subset-dimension
 * property of each element. Subsets are only assigned to previously unassigned elements.*/
template <class TIterator>
static void CopySubsetFromHigherDimNbr(
				ISubsetHandler& sh,
				TIterator elemsBegin,
				TIterator elemsEnd,
				AChar aDimension)
{
	typedef typename PtrToValueType<typename TIterator::value_type>::base_type TElem;

	UG_COND_THROW(!sh.grid(), "The subset-handler has to operate on a grid!");
	Grid& grid = *sh.grid();
	MultiElementAttachmentAccessor<AChar> aaDim(grid, aDimension, true, true, true, true);

	Grid::edge_traits::secure_container		edges;
	Grid::face_traits::secure_container		faces;
	Grid::volume_traits::secure_container	vols;

	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		TElem* elem = *iter;
		if(sh.get_subset_index(elem) != -1)
			continue;

		const int dim = aaDim[elem];
		switch(dim){
			case 1:{
			//	get the edge with dim 1 and with the lowest subset index
				int si = sh.num_subsets();
				grid.associated_elements(edges, elem);
				for(size_t i = 0; i < edges.size(); ++i){
					Edge* e = edges[i];
					int esi = sh.get_subset_index(e);
					if((aaDim[e] == 1) && (esi < si))
						si = esi;
				}
				sh.assign_subset(elem, si);
			}break;

			case 2:{
			//	get the face with dim 2 and with the lowest subset index
				int si = sh.num_subsets();
				grid.associated_elements(faces, elem);
				for(size_t i = 0; i < faces.size(); ++i){
					Face* e = faces[i];
					int esi = sh.get_subset_index(e);
					if((aaDim[e] == 2) && (esi < si))
						si = esi;
				}
				sh.assign_subset(elem, si);
			}break;

			case 3:{
			//	get the volume with dim 3 and with the lowest subset index
				int si = sh.num_subsets();
				grid.associated_elements(vols, elem);
				for(size_t i = 0; i < vols.size(); ++i){
					Volume* e = vols[i];
					int esi = sh.get_subset_index(e);
					if((aaDim[e] == 3) && (esi < si))
						si = esi;
				}
				sh.assign_subset(elem, si);
			}break;
		}
	}
}

void CopySubsetIndicesToSides(
		ISubsetHandler& sh,
		GridObjectCollection goc,
		bool toUnassignedOnly)
{
	UG_COND_THROW(!sh.grid(), "The subset-handler has to operate on a grid!");
	Grid& grid = *sh.grid();

	if(toUnassignedOnly){
		AChar aDimension;
		grid.attach_to_all(aDimension);
		ComputeLocalSubsetDimensions(sh, aDimension, true);
		for(int i = 0; i < (int)goc.num_levels(); ++i)
			CopySubsetFromHigherDimNbr(sh, goc.begin<Vertex>(i), goc.end<Vertex>(i), aDimension);
		for(int i = 0; i < (int)goc.num_levels(); ++i)
			CopySubsetFromHigherDimNbr(sh, goc.begin<Edge>(i), goc.end<Edge>(i), aDimension);
		for(int i = 0; i < (int)goc.num_levels(); ++i)
			CopySubsetFromHigherDimNbr(sh, goc.begin<Face>(i), goc.end<Face>(i), aDimension);
		grid.detach_from_all(aDimension);
	}
	else{
		for(int i = 0; i < (int)goc.num_levels(); ++i){
			CopySubsetIndicesToSides(sh, goc.begin<Volume>(i), goc.end<Volume>(i), false);
			CopySubsetIndicesToSides(sh, goc.begin<Face>(i), goc.end<Face>(i), false);
			CopySubsetIndicesToSides(sh, goc.begin<Edge>(i), goc.end<Edge>(i), false);
			CopySubsetIndicesToSides(sh, goc.begin<Vertex>(i), goc.end<Vertex>(i), false);
		}
	}
}

void CopySubsetIndicesToSides(ISubsetHandler& sh, ISelector& sel,
							  bool toUnassignedOnly)
{
	CopySubsetIndicesToSides (sh, sel.get_grid_objects(), toUnassignedOnly);
}


////////////////////////////////////////////////////////////////////////
void CopySubsetIndicesToSides(ISubsetHandler& sh, bool toUnassignedOnly)
{
	UG_COND_THROW(!sh.grid(), "The subset-handler has to operate on a grid!");
	Grid& grid = *sh.grid();
	CopySubsetIndicesToSides(sh, grid.get_grid_objects(), toUnassignedOnly);
}


////////////////////////////////////////////////////////////////////////
// AdjustSubsetsForLgmNg
void AdjustSubsetsForLgmNg(Grid& grid, SubsetHandler& sh,
							bool keepExistingInterfaceSubsets)
{
	if(grid.num_volumes() > 0)
	{
	//	the 3d case
		vector<Volume*> vVolumes;

	//	first make sure that all volumes are assigned to subsets,
	//	and that no empty volume-subset is between filled ones.
		MakeSubsetsConsecutive<Volume>(sh);

	//	the first index for new subsets.
		int newVolSubsetIndex = GetMaxSubsetIndex<Volume>(sh) + 1;
	//	iterate through all volumes and assign each, that is not assigned
	//	to a subset, to the new one.
		for(VolumeIterator iter = grid.volumes_begin();
			iter != grid.volumes_end(); ++iter)
		{
			if(sh.get_subset_index(*iter) == -1)
				sh.assign_subset(*iter, newVolSubsetIndex);
		}

	//	assign all edges to subset -1
		sh.assign_subset(grid.edges_begin(), grid.edges_end(), -1);

		if(keepExistingInterfaceSubsets){
		//	now we'll assign all faces which are no interface
		//	elements, to subset -1.
			for(FaceIterator iter = grid.faces_begin();
				iter != grid.faces_end(); ++iter)
			{
				Face* f = *iter;
				if(sh.get_subset_index(f) != -1)
				{
					CollectVolumes(vVolumes, grid, f);

					if(vVolumes.size() > 1)
					{
					//	check if there are different adjacent volumes.
						int si = sh.get_subset_index(vVolumes[0]);
						bool gotOne = false;
						for(uint i = 0; i < vVolumes.size(); ++i)
						{
							if(sh.get_subset_index(vVolumes[i]) != si)
							{
								gotOne = true;
								break;
							}
						}

						if(!gotOne)
						{
						//	no, they are all from the same subset.
							sh.assign_subset(f, -1);
						}
					}
				}
			}
		}
		else{
		//	since we don't have to keep existing interface subsets we
		//	can assign all faces to subset -1
			sh.assign_subset(grid.faces_begin(), grid.faces_end(), -1);
		}

	//	make sure that there are no empty face-subsets between filled ones.
	//	since we may not swap subsets (otherwise the volume-sequence would be destroyed)
	//	we have to copy the elements.
		for(int i = 0; i < sh.num_subsets(); ++i){
			if(sh.num<Face>(i) == 0){
			//	find the next that has faces
				int next = sh.num_subsets();
				for(int j = i+1; j < sh.num_subsets(); ++j){
					if(sh.num<Face>(j) > 0){
						next = j;
						break;
					}
				}

			//	if a filled one has been found, we'll copy the elements
				if(next < sh.num_subsets()){
				//	assign all faces in next to subset i
					sh.assign_subset(sh.begin<Face>(next), sh.end<Face>(next), i);
				}
				else{
				//	we're done
					break;
				}
			}
		}

	//	now we have to assign the interface faces to subsets.
		AssignVolumeInterfaceFacesToSubsets(grid, sh);

	//	we now have to make sure that all face-subsets are regular manifolds.
		for(int i = 0; i < sh.num_subsets(); ++i){
			int firstFree = GetMaxSubsetIndex<Face>(sh) + 1;
			SplitIrregularManifoldSubset(sh, i, firstFree, true);
		}

	//	fix orientation of all face subsets
		for(int i = 0; i < sh.num_subsets(); ++i){
			if(sh.num<Face>(i) != 0){
				FixFaceOrientation(grid, sh.begin<Face>(i), sh.end<Face>(i));
			}
		}
	}
	else
	{
	//	the 2d - case.
		vector<Face*> vFaces;

	//	first make sure that all faces are assigned to subsets,
	//	and that no empty face-subset is between filled ones.
		MakeSubsetsConsecutive<Face>(sh);

	//	the first index for new subsets.
		int newFaceSubsetIndex = GetMaxSubsetIndex<Face>(sh) + 1;

	//	iterate through all faces and assign each, that is not assigned
	//	to a subset, to the new one.
		for(FaceIterator iter = grid.faces_begin();
			iter != grid.faces_end(); ++iter)
		{
			if(sh.get_subset_index(*iter) == -1)
				sh.assign_subset(*iter, newFaceSubsetIndex);
		}

	//	now we'll assign all edges which are no interface
	//	elements, to subset -1.
		if(keepExistingInterfaceSubsets){
			for(EdgeIterator iter = grid.edges_begin();
				iter != grid.edges_end(); ++iter)
			{
				Edge* e = *iter;
				if(sh.get_subset_index(e) != -1)
				{
					CollectFaces(vFaces, grid, e);

					if(vFaces.size() > 1)
					{
					//	check if there are different adjacent volumes.
						int si = sh.get_subset_index(vFaces[0]);
						bool gotOne = false;
						for(uint i = 0; i < vFaces.size(); ++i)
						{
							if(sh.get_subset_index(vFaces[i]) != si)
							{
								gotOne = true;
								break;
							}
						}

						if(!gotOne)
						{
						//	no, they are all from the same subset.
							sh.assign_subset(e, -1);
						}
					}
				}
			}
		}
		else{
		//	since we don't have to keep existing interface subsets we can
		//	assign all edges to subset -1.
			sh.assign_subset(grid.edges_begin(), grid.edges_end(), -1);
		}

	//	make sure that there are no empty edge-subsets between filled ones.
	//	since we may not swap subsets (otherwise the face-sequence would be destroyed)
	//	we have to copy the elements.
		for(int i = 0; i < sh.num_subsets(); ++i){
			if(sh.num<Edge>(i) == 0){
			//	find the next that has faces
				int next = sh.num_subsets();
				for(int j = i+1; j < sh.num_subsets(); ++j){
					if(sh.num<Edge>(j) > 0){
						next = j;
						break;
					}
				}

			//	if a filled one has been found, we'll copy the elements
				if(next < sh.num_subsets()){
				//	assign all edges in next to subset i
					sh.assign_subset(sh.begin<Edge>(next), sh.end<Edge>(next), i);
				}
				else{
				//	we're done
					break;
				}
			}
		}

	//	now we have to assign the interface edges to subsets.
		AssignFaceInterfaceEdgesToSubsets(grid, sh);
		
	//	make sure that all edge-subsets are regular poly-chains
		int firstFree = GetMaxSubsetIndex<Edge>(sh) + 1;
		for(int i = 0; i < sh.num_subsets(); ++i){
			if(SplitIrregularPolyChain(sh, i, firstFree))
				++firstFree;
		}

	//	Since ug3 does not support loops, we have to remove those.
		for(int i = 0; i < sh.num_subsets(); ++i){
			size_t chainType = GetPolyChainType(grid, sh.begin<Edge>(i),
												sh.end<Edge>(i),
												IsInSubset(sh, i));
			if(chainType == PCT_CLOSED){
			//	since the chain is regular (see loop above) and since it is
			//	closed, it is enough to simply move the first edge of the
			//	chain to a new subset.
				sh.assign_subset(*sh.begin<Edge>(i), firstFree++);
			}
		}
		
	//	fix orientation
		for(int i = 0; i < sh.num_subsets(); ++i){
			if(sh.num<Edge>(i) != 0){
				FixEdgeOrientation(grid, sh.begin<Edge>(i), sh.end<Edge>(i));
			}
		}

	}
}

////////////////////////////////////////////////////////////////////////
bool SplitIrregularManifoldSubset(SubsetHandler& sh, int srcIndex,
								  int targetIndex, bool strictSplitting)
{
//	Begin with the first face in the subset-handler and find
//	associated faces which build a regular manifold.
//	if the subset is empty, there's nothing to do.
	if(sh.empty<Face>(srcIndex))
		return false;

//	get the grid behind the subset-handler
	if(!sh.grid()){
		UG_LOG("ERROR in SplitIrregularManifoldSubset: No Grid associated");
		UG_LOG(" with the given SubsetHandler.\n");
		return false;
	}

	Grid& grid = *sh.grid();

//	edges are required
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in SplitIrregularManifoldSubset: Autoenabling FACEOPT_AUTOGENERATE_EDGES.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	here we'll store some temporary results
	vector<Edge*> edges;
	vector<Face*> faces;

//	the queue is used to keep a list of elements which have to be checked.
//	begin with the first face.
	queue<Face*> queFaces;
	queFaces.push(*sh.begin<Face>(srcIndex));

//	We have to mark all faces which have once been pushed to the stack
	grid.begin_marking();
	grid.mark(queFaces.front());

//	If we notice that a face would cause a manifold to be irregular,
//	we'll mark it and push it to this vector.
	vector<Face*> vIrregFaces;

//	This counter tells us whether all faces have been processed.
//	If so, the surface is already a regular manifold.
	size_t numProcessed = 0;

//	now while there are unprocessed manifold faces
	while(!queFaces.empty())
	{
		Face* curFace = queFaces.front();
		queFaces.pop();
		++numProcessed;

	//	collect associated faces of associated edges
		CollectAssociated(edges, grid, curFace);
		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge)
		{
		//	collect associated faces
			CollectAssociated(faces, grid, edges[i_edge]);

		//	we only have to continue if we found exactly one face
		//	in the same subset other than curFace.
		//	(Do not check marks here, since we otherwise may find
		//	unwanted neighbors).
			Face* nbr = NULL;
			size_t numFound = 0;
			for(size_t i_face = 0; i_face < faces.size(); ++i_face){
				Face* f = faces[i_face];
				if((f != curFace) &&
					(	(strictSplitting && (sh.get_subset_index(f) != -1))
					 || (!strictSplitting && (sh.get_subset_index(f) == srcIndex))))
				{
				//	we found an associated face in the same subset
					nbr = f;
					numFound++;
				}
			}

		//	If we found exactly one valid neighbor
			if(numFound == 1){
			//	if the neighbor was not already processed we'll schedule it
				if(!grid.is_marked(nbr)){
					grid.mark(nbr);
					queFaces.push(nbr);
				}
			}
			else if(numFound > 1){
			//	here we have to make sure that no face which builds an
			//	irregular manifold with another face can be identified
			//	as a valid neighbor through another element.
			//	Thus we'll mark them. Since they wouldn't be added to the
			//	subset at targetInd later on (since they are marked) we have
			//	to additionally push them into a vector.
				for(size_t i_face = 0; i_face < faces.size(); ++i_face){
					Face* f = faces[i_face];
					if((f != curFace) && (sh.get_subset_index(f) == srcIndex)){
						if(!grid.is_marked(f)){
							grid.mark(f);
							vIrregFaces.push_back(f);
						}
					}
				}
			}
		}
	}

//	if all faces of the subset have been processed, then the whole subset is
//	a regular manifold.
	if(numProcessed == sh.num<Face>(srcIndex)){
		grid.end_marking();
		return false;
	}

//	if not, we'll iterate over all faces of the subset and assign unmarked ones
//	to the subset at targetIndex.
	for(FaceIterator iter = sh.begin<Face>(srcIndex);
		iter != sh.end<Face>(srcIndex);)
	{
	//	Take care with the iterator, since the element may be assigned
	//	to another subset.
		Face* f = *iter;
		++iter;

		if(!grid.is_marked(f))
			sh.assign_subset(f, targetIndex);
	}

//	additionally we have to assign all faces in vIrregFaces to the new subset
	for(size_t i = 0; i < vIrregFaces.size(); ++i)
		sh.assign_subset(vIrregFaces[i], targetIndex);

//	The subset was split - return true.
	grid.end_marking();
	return true;
}

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
//	separates faces by orthogonal axis-aligned normals.
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								APosition aPos, ANormal* paNorm,
								int applyToSubset)
{
	vector<vector3> vNormals(6);
	vNormals[0] = vector3(1, 0, 0);
	vNormals[1] = vector3(0, 1, 0);
	vNormals[2] = vector3(0, 0, 1);
	vNormals[3] = vector3(-1, 0, 0);
	vNormals[4] = vector3(0, -1, 0);
	vNormals[5] = vector3(0, 0, -1);

	SeparateFaceSubsetsByNormal(grid, sh, vNormals, aPos, paNorm,
								applyToSubset);
}

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
//	separates subset by the given normals.
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								std::vector<vector3> vNormals,
								APosition aPos, ANormal* paNorm,
								int applyToSubset)
{
//	iterate through all existing face-subsets.
//		iterate through all faces of the subset
//			find the normal in vNormals that is the closest to the face-normal.
//				add the faces to the matching subset.

	if(!grid.has_vertex_attachment(aPos))
		return;

	if(vNormals.empty())
		return;

//	the position accessor
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	the normal accessor
	vector3 tmpNorm;//	this is only required if paNorm == NULL.
	Grid::FaceAttachmentAccessor<ANormal> aaNorm;
	if(paNorm)
	{
		if(grid.has_face_attachment(*paNorm))
			aaNorm.access(grid, *paNorm);
	}

//	only separate faces of subsets that already existed at the beginning.
	int numFaceSubsets = GetMaxSubsetIndex<Face>(sh) + 1;
	int nextFreeSubset = numFaceSubsets;

//	if a subset was specified, we'll only apply the separation to this subset
	int i = 0;
	int iMax = numFaceSubsets - 1;
	if(applyToSubset >= 0){
		i = applyToSubset;
		iMax = applyToSubset;
	}

	for(; i <= iMax; ++i)
	{
		bool firstFace = true;
	//	this vector holds associated subset indices for each normal.
		vector<int> vSubsetIndices(vNormals.size(), -1);

	//	iterate through the faces of the active subset
		FaceIterator iter = sh.begin<Face>(i);
		while(iter != sh.end<Face>(i))
		{
			Face* f = *iter;
			iter++;//	increment iterator here, since it would else be invalidated later on.

			vector3* pN;//	pointer to the normal of the face.
			if(aaNorm.valid())
				pN = &aaNorm[f];
			else
			{
			//	calculate normal on the fly
				CalculateNormal(tmpNorm, f, aaPos);
				pN = &tmpNorm;
			}

		//	find the index of the closest normal
			int closestInd = -1;
			number closestDist = 100;// 4 is the max-distance-square for normized normals.
			for(uint j = 0; j < (uint)vNormals.size(); ++j)
			{
				number nDist = VecDistanceSq(*pN, vNormals[j]);
				if(nDist < closestDist)
				{
					closestDist = nDist;
					closestInd = j;
				}
			}

		//	get the index of the matching subset.
			if(vSubsetIndices[closestInd] == -1)
			{
			//	all faces with the same normal as the first face
			//	will stay in the subset
				if(firstFace)
				{
					firstFace = false;
					vSubsetIndices[closestInd] = i;
				}
				else
				{
				//	choose an empty subset
					vSubsetIndices[closestInd] = nextFreeSubset++;
				}
			}

		//	assign the face
			if(vSubsetIndices[closestInd] != i)
				sh.assign_subset(f, vSubsetIndices[closestInd]);
		}
	}

}

/*
void SeparateVolumesByFaceSubsets(Grid& grid, SubsetHandler& sh)
{
//	first we'll assign all volumes to subset -1.
	sh.assign_subset(grid.volumes_begin(), grid.volumes_end(), -1);

//	we'll keep all unassigned volumes in a selector.
	Selector sel(grid);
	sel.select(grid.volumes_begin(), grid.volumes_end());

//	those vectors will be used to gather element neighbours.
	vector<Face*> vFaces;
	vector<Volume*> vVolumes;

//	this stack contains all volumes that we still have to check for neighbours.
	stack<Volume*> stkVols;

//	will be used to store sides of volumes
	FaceDescriptor fd;

//	now - while there are unassigned volumes.
	int subsetIndex = 0;
	while(!sel.empty())
	{
	//	choose the volume with which we want to start
	//	TODO: if material-points are supplied, this should be the
	//		the volume that contains the i-th material point.
		stkVols.push(*sel.begin<Volume>());
		while(!stkVols.empty())
		{
			Volume* v = stkVols.top();
			stkVols.pop();
		//	if the volume is unselected it has already been processed.
			if(!sel.is_selected(v))
				continue;
			sel.deselect(v);

		//	assign v to its new subset
			sh.assign_subset(v, subsetIndex);

		//	check neighbour-volumes, whether they belong to the same subset.
		//	iterate through the sides of the volume
			for(uint i = 0; i < v->num_faces(); ++i)
			{
				v->face(i, fd);

			//	get the corresponding face
				Face* f = grid.get_face(fd);

				bool bSeparator = false;
			//	if it belongs to a subset other that -1, it is a separator.
				if(f)
					bSeparator = (sh.get_subset_index(f) != -1);

			//	if fd is not corresponding to a separator, we'll add all connected volumes
				if(!bSeparator)
				{
					CollectVolumes(vVolumes, grid, fd);

				//	add all volumes that are still selected (v is not selected anymore).
					for(uint j = 0; j < vVolumes.size(); ++j)
					{
						if(sel.is_selected(vVolumes[j]))
							stkVols.push(vVolumes[j]);
					}
				}
			}
		}
	//	the stack is empty. increase subset index.
		subsetIndex++;
	}
}
*/
void AssignRegionToSubset(Grid& grid, ISubsetHandler& shVolsOut,
						  const ISubsetHandler& shFaces,
						  Volume* proxyVol, int newSubsetIndex)
{
//	we'll utilize a stack to gather all volumes that have to
//	be examined.
	stack<Volume*> stkVols;
	stkVols.push(proxyVol);

//	this vector will be used to collect neighbours
	vector<Volume*> vVols;

//	while there are volumes on the stack we'll go on
	while(!stkVols.empty()){
		Volume* v = stkVols.top();
		stkVols.pop();

		shVolsOut.assign_subset(v, newSubsetIndex);

	//	check all neighbours and decide whether they have to be
	//	pushed to the stack.
	//	First we'll have to check for each side of v whether we may traverse it
		for(uint i = 0; i < v->num_faces(); ++i){
			Face* f = grid.get_face(v, i);
			if(f){
			//	check whether f lies in a subset
			//	if so we may not traverse it and we'll continue with the next face.
				if(shFaces.get_subset_index(f) != -1)
					continue;

			//	we may traverse it. get the neighbour-volumes
				CollectVolumes(vVols, grid, f);
			}
			else{
			//	no associated face existed. We'll use GetNeighbours to find
			//	neighbouring volumes. Since no face can separate them, they
			//	belong to the region, too.
				GetNeighbours(vVols, grid, v, i);
			}

		//	add the volumes in vVols to the stack - but only if they don't
		//	already belong to the region.
			for(size_t j = 0; j < vVols.size(); ++j){
				if(vVols[j] != v){
					if(shVolsOut.get_subset_index(vVols[j]) != newSubsetIndex)
						stkVols.push(vVols[j]);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
bool SeparateRegions(Grid& grid, ISubsetHandler& shVolsOut,
					 const ISubsetHandler& shFaces,
					 const MarkerPointManager& mpm,
					 int firstSubsetIndex)
{
	if(!grid.has_vertex_attachment(aPosition))
		return false;

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	for(size_t i = 0; i < mpm.num_markers(); ++i)
	{
		const vector3& pos = mpm.get_marker(i).pos;

//TODO: extend to volumes in general. Add a PointIsInsideVolume method.
	//	find the tetrahedron that contains the marker point
		for(TetrahedronIterator iter = grid.begin<Tetrahedron>();
			iter != grid.end<Tetrahedron>(); ++iter)
		{
			if(PointIsInsideTetrahedron(pos, *iter, aaPos)){
			//	assign region to subset
				int si = firstSubsetIndex + i;
				AssignRegionToSubset(grid, shVolsOut, shFaces,
									*iter, si);
			//	set subset name
				shVolsOut.subset_info(si).name = mpm.get_marker(i).name;
			//	we don't have to check other volumes.
				break;
			}
		}
	}

//	done
	return true;
}

void AssignInnerAndBoundarySubsets(Grid& grid, ISubsetHandler& shOut,
									int inSubset, int bndSubset)
{
	if(grid.num<Volume>() > 0){
	//	assign volumes to inSubset
		for(VolumeIterator iter = grid.begin<Volume>();
			iter != grid.end<Volume>(); ++iter)
		{
			if(shOut.get_subset_index(*iter) == -1)
				shOut.assign_subset(*iter, inSubset);
		}

		vector<Edge*> vEdges;
	//	assign volume-boundary elements to bndSubset
		for(FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>(); ++iter)
		{
			Face* f = *iter;
			if(shOut.get_subset_index(f) == -1){
				if(IsVolumeBoundaryFace(grid, f)){
				//	assign the face to the boundary subset.
					shOut.assign_subset(f, bndSubset);

				//	assign associated vertices and edges to the boundary subset, too.
					for(size_t i = 0; i < f->num_vertices(); ++i){
						if(shOut.get_subset_index(f->vertex(i)) == -1)
							shOut.assign_subset(f->vertex(i), bndSubset);
					}

					CollectEdges(vEdges, grid, f);
					for(size_t i = 0; i < vEdges.size(); ++i){
						if(shOut.get_subset_index(vEdges[i]) == -1)
							shOut.assign_subset(vEdges[i], bndSubset);
					}
				}
				else{
				//	assing the face to the inner subset
					shOut.assign_subset(f, inSubset);
				}
			}
		}
	}
	else{
		for(FaceIterator iter = grid.begin<Face>();
			iter != grid.end<Face>(); ++iter)
		{
			if(shOut.get_subset_index(*iter) == -1)
				shOut.assign_subset(*iter, inSubset);
		}
	}

//	assign edges and vertices
	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(shOut.get_subset_index(e) == -1){
			if(IsBoundaryEdge2D(grid, e))
			{
				shOut.assign_subset(e, bndSubset);
				for(size_t i = 0; i < 2; ++i)
				{
					if(shOut.get_subset_index(e->vertex(i)) == -1)
						shOut.assign_subset(e->vertex(i), bndSubset);
				}
			}
			else
				shOut.assign_subset(e, inSubset);
		}
	}

//	assign unassigned vertices
	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter)
	{
		if(shOut.get_subset_index(*iter) == -1)
			shOut.assign_subset(*iter, inSubset);
	}
}


////////////////////////////////////////////////////////////////////////
vector3 GetColorFromStandardPalette(int index)
{
	return GetColorFromDefaultPalette (index);
}

////////////////////////////////////////////////////////////////////////
void AssignSubsetColors(ISubsetHandler& sh)
{
	AssignDefaultSubsetColors (sh);
}


////////////////////////////////////////////////////////////////////////
template <class TElem>
void AssignSidesToSubsets(ISubsetHandler& sh, ISelector* psel)
{
	typedef typename TElem::lower_dim_base_object 		Side;
	typedef typename geometry_traits<Side>::iterator 	SideIter;
	typedef basic_string<int> IntString;

//	access the grid on which sh operates.
	if(!sh.grid())
		return;

	Grid& grid = *sh.grid();

//	we'll use those marks to check whether a subset has already been
//	processed in an iteration. A subset is considered to be marked, if
//	marks[subInd] == curMarker.
//	Note that curMarker is increased in each iteration.
//	Note also, that marks has to be resized whenever a new subset is added.
	size_t curMark = 1;
	vector<size_t> marks(sh.num_subsets(), 0);

//	here we'll check whether a subset already contains elements with
//	a given neighborhood subset-constellation.
//	Note that the key will be calculated from the subset-constellation.
//	we're using double to allow more values
	map<IntString, int> subsetMap;

//	this vector will be used to collect all associated subsets of an element.
//	we'll use the subset-mark mechanism to avoid duplicate entries.
	IntString skey;

//	used to collect neighbors of type TElem
	vector<TElem*> elems;

//	Now iterate over all elements of type Side
	for(SideIter iterSide = grid.begin<Side>();
		iterSide != grid.end<Side>(); ++iterSide)
	{
		Side* side = *iterSide;

		if(sh.get_subset_index(side) != -1)
			continue;

	//	check whether we are working on a selected side
		bool isSelected = false;
		if(psel && psel->is_selected(side))
			isSelected = true;

	//	increase curMark
		++curMark;

	//	collect all associated elements of side
		CollectAssociated(elems, grid, side);

	//	collect the subsets in which they lie. Note that we won't push the
	//	same subset index twice to skey.
		skey.clear();

	//	this subset index shall be used if the key is new
		int useSubsetInd = GetFirstFreeSubset(sh);

		for(size_t i_elems = 0; i_elems < elems.size(); ++i_elems){
			bool nbrIsSelected = false;
			if(psel && psel->is_selected(elems[i_elems]))
				nbrIsSelected = true;

		//	either the element is not selected or it is selected and its neighbor
		//	is selected, too.
			if(!isSelected || nbrIsSelected){
				int si = sh.get_subset_index(elems[i_elems]);
			//	check if the subset is marked
				if(marks[si] != curMark){
				//	no - mark it and add si to subsets
					marks[si] = curMark;
					skey += si;
				}
			}
		}

		if(skey.size() == 1){
			if(elems.size() > 1){
			//	we've got an inner element. Assign the subset-index
			//	of those inner neighbor elements.
				useSubsetInd = skey[0];
			}
			else{
			//	An outer interface element.
			//	add -2 to skey to make sure that they go into a separate subset
			//	(-1 might occur as a normal subset index)
				skey += -2;
			}
		}

	//	The subset indices in skey have to be sorted before they can be
	//	used as unique key
		sort(skey.begin(), skey.end());

	//	get the subset index from the subsetMap, if one exists.
		map<IntString, int>::iterator tMapIter = subsetMap.find(skey);
		if(tMapIter != subsetMap.end()){
		//	got one
			sh.assign_subset(side, tMapIter->second);
		}
		else{
		//	create a new entry and resize the marks array
			subsetMap[skey] = useSubsetInd;
			sh.assign_subset(side, useSubsetInd);
			if(sh.num_subsets() > (int)marks.size())
				marks.resize(sh.num_subsets(), 0);
		}
	}
}

//	template specialization
template void AssignSidesToSubsets<Edge>(ISubsetHandler&, ISelector*);
template void AssignSidesToSubsets<Face>(ISubsetHandler&, ISelector*);
template void AssignSidesToSubsets<Volume>(ISubsetHandler&, ISelector*);


void AssignSubsetsByElementType(ISubsetHandler& sh)
{
	if(!sh.grid())
		return;

	AssignSubsetsByElementType(sh, sh.grid()->get_grid_objects());
}


void AssignSubsetsByElementType(ISubsetHandler& sh, GridObjectCollection g)
{
	int subsetInd = 0;

	if(g.num<RegularVertex>() > 0){
		sh.assign_subset(g.begin<RegularVertex>(), g.end<RegularVertex>(), subsetInd);
		sh.subset_info(subsetInd++).name = "RegularVertex";
	}

	if(g.num<ConstrainedVertex>() > 0){
		sh.assign_subset(g.begin<ConstrainedVertex>(), g.end<ConstrainedVertex>(), subsetInd);
		sh.subset_info(subsetInd++).name = "HangingVertex";
	}

	if(g.num<RegularEdge>() > 0){
		sh.assign_subset(g.begin<RegularEdge>(), g.end<RegularEdge>(), subsetInd);
		sh.subset_info(subsetInd++).name = "RegularEdge";
	}

	if(g.num<ConstrainingEdge>() > 0){
		sh.assign_subset(g.begin<ConstrainingEdge>(), g.end<ConstrainingEdge>(), subsetInd);
		sh.subset_info(subsetInd++).name = "ConstrainingEdge";
	}

	if(g.num<ConstrainedEdge>() > 0){
		sh.assign_subset(g.begin<ConstrainedEdge>(), g.end<ConstrainedEdge>(), subsetInd);
		sh.subset_info(subsetInd++).name = "ConstrainedEdge";
	}

	if(g.num<Triangle>() > 0){
		sh.assign_subset(g.begin<Triangle>(), g.end<Triangle>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Triangle";
	}

	if(g.num<ConstrainingTriangle>() > 0){
		sh.assign_subset(g.begin<ConstrainingTriangle>(), g.end<ConstrainingTriangle>(), subsetInd);
		sh.subset_info(subsetInd++).name = "ConstrainingTriangle";
	}

	if(g.num<ConstrainedTriangle>() > 0){
		sh.assign_subset(g.begin<ConstrainedTriangle>(), g.end<ConstrainedTriangle>(), subsetInd);
		sh.subset_info(subsetInd++).name = "ConstrainedTriangle";
	}

	if(g.num<Quadrilateral>() > 0){
		sh.assign_subset(g.begin<Quadrilateral>(), g.end<Quadrilateral>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Quadrilateral";
	}

	if(g.num<ConstrainingQuadrilateral>() > 0){
		sh.assign_subset(g.begin<ConstrainingQuadrilateral>(), g.end<ConstrainingQuadrilateral>(), subsetInd);
		sh.subset_info(subsetInd++).name = "ConstrainingQuadrilateral";
	}

	if(g.num<ConstrainedQuadrilateral>() > 0){
		sh.assign_subset(g.begin<ConstrainedQuadrilateral>(), g.end<ConstrainedQuadrilateral>(), subsetInd);
		sh.subset_info(subsetInd++).name = "ConstrainedQuadrilateral";
	}

	if(g.num<Tetrahedron>() > 0){
		sh.assign_subset(g.begin<Tetrahedron>(), g.end<Tetrahedron>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Tetrahedron";
	}

	if(g.num<Pyramid>() > 0){
		sh.assign_subset(g.begin<Pyramid>(), g.end<Pyramid>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Pyramid";
	}

	if(g.num<Prism>() > 0){
		sh.assign_subset(g.begin<Prism>(), g.end<Prism>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Prism";
	}

	if(g.num<Hexahedron>() > 0){
		sh.assign_subset(g.begin<Hexahedron>(), g.end<Hexahedron>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Hexahedron";
	}

	if(g.num<Octahedron>() > 0){
		sh.assign_subset(g.begin<Octahedron>(), g.end<Octahedron>(), subsetInd);
		sh.subset_info(subsetInd++).name = "Octahedron";
	}
}

template <class TElem>
static char GetSmallestLocalSubsetDimension(
				typename Grid::traits<TElem>::secure_container& nbrs,
				MultiElementAttachmentAccessor<AChar>& aaDim)
{
	int d = -1;
	for(size_t i = 0; i < nbrs.size(); ++i){
		TElem* e = nbrs[i];
		char ed = aaDim[e];
		if((ed != -1 ) && ((d == -1) || (ed < d)))
			d = ed;
	}
	return d;
}

template <class TElem>
static bool NbrIsInSubset(
				ISubsetHandler& sh,
				typename Grid::traits<TElem>::secure_container& nbrs,
				int si)
{
	for(size_t i = 0; i < nbrs.size(); ++i){
		if(sh.get_subset_index(nbrs[i]) == si)
			return true;
	}
	return false;
}

void ComputeLocalSubsetDimensions(
		ISubsetHandler& sh, 
		AChar aDimension,
		bool includeUnassigned)
{
	UG_COND_THROW(!sh.grid(), "The subset-handler has to operate on a grid!");
	Grid& grid = *sh.grid();
	grid.attach_to_all(aDimension);
	MultiElementAttachmentAccessor<AChar> aaDim(grid, aDimension, true, true, true, true);


	for(VolumeIterator ivol = grid.begin<Volume>(); ivol != grid.end<Volume>(); ++ivol){
		if(sh.get_subset_index(*ivol) == -1)
			aaDim[*ivol] = -1;
		else
			aaDim[*ivol] = 3;
	}

	Grid::volume_traits::secure_container vols;
	for(FaceIterator iface = grid.begin<Face>(); iface != grid.end<Face>(); ++iface){
		Face* f = *iface;
		int si = sh.get_subset_index(f);
		if(si == -1){
			if(includeUnassigned){
				grid.associated_elements(vols, f);
				aaDim[f] = GetSmallestLocalSubsetDimension<Volume>(vols, aaDim);
			}
			else
				aaDim[f] = -1;
		}
		else{
			aaDim[f] = 2;
			grid.associated_elements(vols, f);
			if(NbrIsInSubset<Volume>(sh, vols, si))
				aaDim[f] = 3;
		}
	}

	Grid::face_traits::secure_container faces;
	for(EdgeIterator iedge = grid.begin<Edge>(); iedge != grid.end<Edge>(); ++iedge){
		Edge* e = *iedge;
		int si = sh.get_subset_index(e);
		if(si == -1){
			if(includeUnassigned){
				grid.associated_elements(faces, e);
				aaDim[e] = GetSmallestLocalSubsetDimension<Face>(faces, aaDim);
			}
			else
				aaDim[e] = -1;
		}
		else{
			aaDim[e] = 1;
			grid.associated_elements(vols, e);
			if(NbrIsInSubset<Volume>(sh, vols, si))
				aaDim[e] = 3;
			else{
				grid.associated_elements(faces, e);
				if(NbrIsInSubset<Face>(sh, faces, si))
					aaDim[e] = 2;
			}
		}
	}

	Grid::edge_traits::secure_container edges;
	for(VertexIterator ivrt = grid.begin<Vertex>(); ivrt != grid.end<Vertex>(); ++ivrt){
		Vertex* v = *ivrt;
		int si = sh.get_subset_index(v);
		if(si == -1){
			if(includeUnassigned){
				grid.associated_elements(edges, v);
				aaDim[v] = GetSmallestLocalSubsetDimension<Edge>(edges, aaDim);
			}
			else
				aaDim[v] = -1;
		}
		else{
			aaDim[v] = 0;
			grid.associated_elements(vols, v);
			if(NbrIsInSubset<Volume>(sh, vols, si))
				aaDim[v] = 3;
			else{
				grid.associated_elements(faces, v);
				if(NbrIsInSubset<Face>(sh, faces, si))
					aaDim[v] = 2;
				else{
					grid.associated_elements(edges, v);
					if(NbrIsInSubset<Edge>(sh, edges, si))
						aaDim[v] = 1;
				}
			}
		}
	}
}

}//	end of namespace

