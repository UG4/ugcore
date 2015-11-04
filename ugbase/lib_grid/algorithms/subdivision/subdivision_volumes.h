// created by mstepnie
// martin.stepniewski@gcsc.uni-frankfurt.de
// Juli 14, 2014

#ifndef __H__UG__SUBDIVISION_VOLUMES__
#define __H__UG__SUBDIVISION_VOLUMES__

#include <vector>
#include <cassert>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"

namespace ug
{


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

		//	Tetrahedron pattern from tetrahedron_rules.cpp
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

		//	Tetrahedron pattern from tetrahedron_rules.cpp
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

		//	Tetrahedron pattern from tetrahedron_rules.cpp
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
 * 	@param vrt			pointer to Vertex
 * 	@param aaPos		reference to VertexAttachmentAccessor accessing aPos
 * 	@param aaSmoothPos	reference to VertexAttachmentAccessor accessing aSmoothPos
**/
void MoveVertexToSubdivisionVolumePosition(	MultiGrid& mg, Vertex* vrt,
											Grid::VertexAttachmentAccessor<APosition>& aaPos,
											Grid::VertexAttachmentAccessor<APosition>& aaSmoothPos)
{
//	Declare centroid coordinate vector
	typedef APosition::ValueType pos_type;
	pos_type p;

//	Declare vertex volume valence
	size_t valence = 0;

//	Collect associated volumes
	std::vector<Volume*> volumes;
	CollectAssociated(volumes, mg, vrt);

//	Iterate over all associated volumes
	for(Grid::AssociatedVolumeIterator vIter = mg.associated_volumes_begin(vrt); vIter != mg.associated_volumes_end(vrt); ++vIter)
	{
		VecSet(p, 0);
		Volume* vol = *vIter;
		++valence;

	//	TETRAHEDRON CASE
		if(vol->reference_object_id() == ROID_TETRAHEDRON)
		{
		//	Iterate over all associated vertices inside tetrahedron
			int vrtIndex = 0;
			for(size_t i = 0; i < vol->num_vertices(); ++i)
			{
				if(vrtIndex != GetVertexIndex(vol, vrt))
				{
					VecAdd(p, p, aaPos[vol->vertex(i)]);
				}

				++vrtIndex;
			}

		//	TODO: refer to subdivision rules object
			number centerWgt 	= -1.0/16;
			number nbrWgt 		= 17.0/48;

			VecScaleAppend(aaSmoothPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p);
		}

	//	OCTAHEDRON CASE
		else if(vol->reference_object_id() == ROID_OCTAHEDRON)
		{
		//	Summate positions of all other vertices inside octahedron, first associate ones and last the opposing one
			for(Grid::AssociatedEdgeIterator eIter = mg.associated_edges_begin(vol); eIter != mg.associated_edges_end(vol); ++eIter)
			{
				Edge* e = *eIter;
				if(GetConnectedVertex(e, vrt) != NULL)
				{
					VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
				}
			}

			Vertex* oppVrt = vol->vertex(vol->get_opposing_object(vrt).second);

		//	TODO: refer to subdivision rules object
			number centerWgt 	= 3.0/8;
			number nbrWgt 		= 1.0/12;
			number oppNbrWgt 	= 7.0/24;

			VecScaleAppend(aaSmoothPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p, oppNbrWgt, aaPos[oppVrt]);
		}

	//	UNSUPPORTED VOLUME ELEMENT CASE
		else
		{
			UG_THROW("Volume type not supported for subdivision volumes refinement.");
		}
	}

//	Scale vertex position by the number of associated volume elements
	VecScale(aaSmoothPos[vrt],  aaSmoothPos[vrt], 1.0/valence);
}


/// Vertex repositioning function for subdivision surfaces refinement (by C. Loop, 1987)
/** This function repositions a vertex to its smoothed position in the current
 * 	level determined by the subdivision surfaces refinement.
 *
 * 	@param mg			reference to MultiGrid
 * 	@param vrt			pointer to Vertex
 * 	@param aaPos		reference to VertexAttachmentAccessor accessing aPos
 * 	@param aaSmoothPos	reference to VertexAttachmentAccessor accessing aSmoothPos
**/
void MoveVertexToSubdivisionSurfacePosition(MultiGrid& mg, Vertex* vrt,
											Grid::VertexAttachmentAccessor<APosition>& aaPos,
											Grid::VertexAttachmentAccessor<APosition>& aaSmoothPos)
{
//	Check, if volumes are included in input grid
	bool volumesExist = mg.num<Volume>() > 0;

//	Declare centroid coordinate vector
	typedef APosition::ValueType pos_type;
	pos_type p;

//	Load subdivision surfaces rules
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	Even vertex
	if(mg.get_parent(vrt)->reference_object_id() == ROID_VERTEX)
	{
		//Vertex* parentVrt = dynamic_cast<Vertex*>(mg.get_parent(vrt));
		Vertex* parentVrt = static_cast<Vertex*>(mg.get_parent(vrt));

	//	perform loop subdivision on even surface vertices
	//	first get neighboured vertices
		size_t valence = 0;
		pos_type p;
		VecSet(p, 0);

		for(Grid::AssociatedEdgeIterator iter = mg.associated_edges_begin(parentVrt); iter != mg.associated_edges_end(parentVrt); ++iter)
		{
			if((!volumesExist) || IsBoundaryEdge3D(mg, *iter))
			{
				VecAdd(p, p, aaPos[GetConnectedVertex(*iter, parentVrt)]);
				++valence;
			}
		}

		number centerWgt 	= subdiv.ref_even_inner_center_weight(valence);
		number nbrWgt 		= subdiv.ref_even_inner_nbr_weight(valence);

		VecScaleAdd(aaSmoothPos[vrt], centerWgt, aaPos[parentVrt], nbrWgt, p);
	}

//	Odd vertex
	else if(mg.get_parent(vrt)->reference_object_id() == ROID_EDGE)
	{
	//	Get parent edge
		//Edge* parentEdge = dynamic_cast<Edge*>(mg.get_parent(vrt));
		Edge* parentEdge = static_cast<Edge*>(mg.get_parent(vrt));

	//	apply loop-subdivision on inner elements
	//	get the neighboured triangles
		Face* f[2];
		int numAssociatedBndFaces = 0;

		std::vector<Face*> faces;
		CollectAssociated(faces, mg, parentEdge);

		for(size_t i = 0; i < faces.size(); ++i)
		{
			if(IsBoundaryFace3D(mg, faces[i]))
			{
				if(numAssociatedBndFaces < 2)
				{
					f[numAssociatedBndFaces] = faces[i];
				}
				++numAssociatedBndFaces;
			}
		}

		if(numAssociatedBndFaces == 2)
		{
			if(f[0]->num_vertices() == 3 && f[1]->num_vertices() == 3)
			{
			//	the 4 vertices that are important for the calculation
				Vertex* v[4];
				v[0] = parentEdge->vertex(0); v[1] = parentEdge->vertex(1);
				v[2] = GetConnectedVertex(parentEdge, f[0]);
				v[3] = GetConnectedVertex(parentEdge, f[1]);

				vector4 wghts;

				wghts = subdiv.ref_odd_inner_weights();

				VecScaleAdd(aaSmoothPos[vrt],
							wghts.x(), aaPos[v[0]], wghts.y(), aaPos[v[1]],
							wghts.z(), aaPos[v[2]], wghts.w(), aaPos[v[3]]);
			}
			else
				UG_THROW("Non triangular faces included in grid.");
		}
		else
			UG_THROW("numAssociatedBndFaces != 2.");
	}
}


/// Function to create a smooth subdivision volumes hierarchy
/** This function transforms a linearly refined hybrid tetra-/octahedral volume
 * 	grid hierarchy into a smoothed subdivision volumes hierarchy
 * 	(s. Schaefer et al, "Smooth subdivision of tetrahedral meshes", 2004)
 *
 * 	@param mg					reference to MultiGrid
 * 	@param bPreserveBnd			bool flag to preserve boundary
 * 	@param bSubdivisionLoopBnd	bool flag to apply subdivision surfaces smoothing
 * 								on the preserved boundary
**/
void SubdivisionVolumes(MultiGrid& mg, bool bPreserveBnd, bool bSubdivisionLoopBnd)
{
//	Ensure, that hybrid tet-/oct refinement is used as refinement rule for tetrahedrons
	if(tet_rules::GetRefinementRule() != tet_rules::HYBRID_TET_OCT)
		UG_THROW("ERROR in SubdivisionVolumes: Set necessary refinement rule by SetTetRefinementRule(hybrid_tet_oct).");

//	Ensure correct user input
	if(!bPreserveBnd)
		bSubdivisionLoopBnd = false;

//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

//	Smooth position attachment plus (default) initialization with zero
	APosition aSmoothPos;
	mg.attach_to_vertices_dv(aSmoothPos, vector3(0, 0, 0));
	Grid::VertexAttachmentAccessor<APosition> aaSmoothPos(mg, aSmoothPos);

//	Load subdivision surfaces rules
	//SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	Check, if volumes are included in input grid
	bool volumesExist = mg.num<Volume>() > 0;
	if(!volumesExist)
		UG_THROW("SubdivisionVolumes: No volumes included in input grid for smooth TetGrid subdivision refinement.");

// 	Loop all vertices of top_level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		if(!bPreserveBnd)
		{
			MoveVertexToSubdivisionVolumePosition(mg, vrt, aaPos, aaSmoothPos);
		}

	//	Surface preservation
		else
		{
		//	BOUNDARY VERTEX
			if(IsBoundaryVertex3D(mg, vrt))
			{
				if(bSubdivisionLoopBnd)
					MoveVertexToSubdivisionSurfacePosition(mg, vrt, aaPos, aaSmoothPos);
				else
					aaSmoothPos[vrt] = aaPos[vrt];
			}

		//	INNER VERTEX
			else
			{
				MoveVertexToSubdivisionVolumePosition(mg, vrt, aaPos, aaSmoothPos);
			}
		}
	}

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		VecScale(aaPos[vrt], aaSmoothPos[vrt], 1.0);
	}

}


}//	end of namespace

#endif
