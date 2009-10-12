//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#include <vector>
#include "subset_util.h"
#include "geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	AssignInterfaceEdgesToSubsets
void AssignFaceInterfaceEdgesToSubsets(Grid& grid, SubsetHandler& sh)
{
//TODO:	This algorithm can be improved.
//	In the moment it runs with arbitrarily O(numSubs^2*averageNumFacesPerSub)

//	get the first free edge-subset
	int edgeSubsetIndex = GetMaxSubsetIndex<EdgeBase>(sh) + 1;

//	use this vector to collect adjacent faces
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;

//	iterate through all subsets
	for(int j = -1; j < (int)sh.num_subsets(); ++j)
	{
	//	find interface-edges to the other subsets.
		for(uint i = j + 1; i < sh.num_subsets(); ++i)
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
					EdgeBase* e = vEdges[k];

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
				if(sh.num_elements<EdgeBase>(edgeSubsetIndex) > 0)
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
// AdjustSubsetsForLgmNg
void AdjustSubsetsForLgmNg(Grid& grid, SubsetHandler& sh)
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

	//TODO: as soon as the subset-handler supports reordering of subsets
	//		for each element-type separately, we should bring the face
	//		subsets in consecutive order.
	//	now we have to assign the interface faces to subsets.
		AssignVolumeInterfaceFacesToSubsets(grid, sh);
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
		for(EdgeBaseIterator iter = grid.edges_begin();
			iter != grid.edges_end(); ++iter)
		{
			EdgeBase* e = *iter;
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

	//TODO: as soon as the subset-handler supports reordering of subsets
	//		for each element-type separately, we should bring the edge
	//		subsets in consecutive order.
	//	now we have to assign the interface edges to subsets.
		AssignFaceInterfaceEdgesToSubsets(grid, sh);
	}
}

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
//	separates faces by orthogonal axis-aligned normals.
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								APosition aPos, ANormal* paNorm)
{
	vector<vector3> vNormals(6);
	vNormals[0] = vector3(1, 0, 0);
	vNormals[1] = vector3(0, 1, 0);
	vNormals[2] = vector3(0, 0, 1);
	vNormals[3] = vector3(-1, 0, 0);
	vNormals[4] = vector3(0, -1, 0);
	vNormals[5] = vector3(0, 0, -1);

	SeparateFaceSubsetsByNormal(grid, sh, vNormals, aPos, paNorm);
}

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
//	separates subset by the given normals.
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								std::vector<vector3> vNormals,
								APosition aPos, ANormal* paNorm)
{
//	iterate through all existing face-subsets.
//		iterate through all faces of the subset
//			find the normal in vNormals that is the closest to the face-normal.
//				add the faces to the matching subset.

	if(!grid.has_vertex_attachment(aPos))
		return;

	if(vNormals.size() == 0)
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

	for(int i = 0; i < numFaceSubsets; ++i)
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

void SeparateVolumesByFaceSubsets(Grid& grid, SubsetHandler& sh,
									vector3* pMaterialPoints,
									int numMaterialPoints)
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
				Face* f = FindFace(grid, fd);

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

}//	end of namespace

