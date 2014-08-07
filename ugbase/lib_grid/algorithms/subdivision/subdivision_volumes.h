// created by mstepnie
// martin.stepniewski@gcsc.uni-frankfurt.de
// Juli 14, 2014

#ifndef __H__UG__SUBDIVISION_VOLUMES__
#define __H__UG__SUBDIVISION_VOLUMES__

#include <vector>
#include <cassert>
#include "lib_grid/lg_base.h"

namespace ug
{

void ProjectToLimitSmoothTetGrid(MultiGrid& mg)
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


////////////////////////////////////////////////////////////////////////////////////////////
//	MoveVertexToSmoothTetGridSubdivisionPosition
void MoveVertexToSmoothTetGridSubdivisionPosition(MultiGrid& mg, Vertex* vrt,
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

			//
			// Alternative iteration
			//
			/*
			for(Grid::AssociatedEdgeIterator eIter = mg.associated_edges_begin(vol); eIter != mg.associated_edges_end(vol); ++eIter)
			{
				Edge* e = *eIter;
				if(GetConnectedVertex(e, vrt) != NULL)
					VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
			}
			*/

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
		//	Iterate over all vertices inside octahedron, first associate ones and last the opposing one
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

////////////////////////////////////////////////////////////////////////////////////////////
//	SubdivisionTetGridSmooth
/// (see Schaefer et al, "Smooth subdivision of tetrahedral meshes")
//void SubdivisionTetGridSmooth(MultiGrid& mg, MGSubsetHandler& sh)
void SubdivisionTetGridSmooth(MultiGrid& mg)
{
	typedef APosition::ValueType pos_type;

//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

//	Smooth position attachment plus initialisation with zero
	APosition aSmoothPosition(0.0);
	mg.attach_to_vertices(aSmoothPosition);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothPos(mg, aSmoothPosition);

	vector3 zeroVector(0.0, 0.0, 0.0);
	SetAttachmentValues(aaSmoothPos, mg.begin<Vertex>(mg.top_level()), mg.end<Vertex>(mg.top_level()), zeroVector);

//	Load subdivision surfaces rules
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	Check, if volumes are included in input grid
	bool volumesExist = mg.num<Volume>() > 0;
	if(!volumesExist)
		UG_THROW("SubdivisionTetGridSmooth: No volumes included in input grid for smooth TetGrid subdivision refinement.");

// 	Loop all vertices of top_level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

	//	BOUNDARY VERTEX
		if(IsBoundaryVertex3D(mg, vrt))
		{
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

	//	INNER VERTEX
		else
		{
			MoveVertexToSmoothTetGridSubdivisionPosition(mg, vrt, aaPos, aaSmoothPos);
		}
	}

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		VecScale(aaPos[vrt], aaSmoothPos[vrt], 1.0);
	}

}

////////////////////////////////////////////////////////////////////////////////////////////
//	SubdivisionTetGridSmoothBasic
void SubdivisionTetGridSmoothBasic(MultiGrid& mg, MGSubsetHandler& sh)
{
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

	APosition aTmpPosition;
	mg.attach_to_vertices(aTmpPosition);
	Grid::VertexAttachmentAccessor<APosition> aaTmpPos(mg, aTmpPosition);

//	Select all elements for linear refinement
	//sel.select(mg.begin<Vertex>(0), mg.end<Vertex>(0));
	//sel.select(mg.begin<Face>(0), mg.end<Face>(0));
	//sel.select(mg.begin<Volume>(0), mg.end<Volume>(0));

	//tet_rules::SetRefinementRule(tet_rules::HYBRID_TET_OCT);
	//Refine(mg, sel);

//	Initialize new TmpPosition of vertices with 0
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	//for(VertexIterator vrtIter = mg.begin<Vertex>(i); vrtIter != mg.end<Vertex>(i); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		aaTmpPos[vrt].x() = 0.0;
		aaTmpPos[vrt].y() = 0.0;
		aaTmpPos[vrt].z() = 0.0;
	}

//	Smooth tetrahedral vertices (see Schaefer et al, "Smooth subdivision of tetrahedral meshes")
	for(VolumeIterator volIter = mg.begin<Tetrahedron>(mg.top_level()); volIter != mg.end<Tetrahedron>(mg.top_level()); ++volIter)
	//for(VolumeIterator volIter = mg.begin<Tetrahedron>(i); volIter != mg.end<Tetrahedron>(i); ++volIter)
	{
		Volume* vol = *volIter;

		VecScaleAppend(	aaTmpPos[vol->vertex(0)],
						-1.0/16, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(2)],
						17.0/48, aaPos[vol->vertex(3)]);

		VecScaleAppend(	aaTmpPos[vol->vertex(1)],
						-1.0/16, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(2)],
						17.0/48, aaPos[vol->vertex(3)]);

		VecScaleAppend(	aaTmpPos[vol->vertex(2)],
						-1.0/16, aaPos[vol->vertex(2)],
						17.0/48, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(3)]);

		VecScaleAppend(	aaTmpPos[vol->vertex(3)],
						-1.0/16, aaPos[vol->vertex(3)],
						17.0/48, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(2)]);
	}

//	Smooth octahedral vertices (see Schaefer et al, "Smooth subdivision of tetrahedral meshes")
	for(VolumeIterator volIter = mg.begin<Octahedron>(mg.top_level()); volIter != mg.end<Octahedron>(mg.top_level()); ++volIter)
	//for(VolumeIterator volIter = mg.begin<Octahedron>(i); volIter != mg.end<Octahedron>(i); ++volIter)
	{
		Volume* vol = *volIter;

		Vertex* vrt1 = vol->vertex(0);
		Vertex* vrt2 = vol->vertex(1);
		Vertex* vrt3 = vol->vertex(2);
		Vertex* vrt4 = vol->vertex(3);
		Vertex* vrt5 = vol->vertex(4);
		Vertex* vrt6 = vol->vertex(5);

	//	1
		VecScaleAppend(	aaTmpPos[vrt1],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt5]);

		VecScaleAppend(	aaTmpPos[vrt1],
						3.0/8,  aaPos[vrt1],
						7.0/24, aaPos[vrt6]);

	//	2
		VecScaleAppend(	aaTmpPos[vrt2],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt5],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt2],
						3.0/8,  aaPos[vrt2],
						7.0/24, aaPos[vrt4]);

	//	3
		VecScaleAppend(	aaTmpPos[vrt3],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt3],
						3.0/8,  aaPos[vrt3],
						7.0/24, aaPos[vrt5]);

	//	4
		VecScaleAppend(	aaTmpPos[vrt4],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt5],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt4],
						3.0/8,  aaPos[vrt4],
						7.0/24, aaPos[vrt2]);

	//	5
		VecScaleAppend(	aaTmpPos[vrt5],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt5],
						3.0/8,  aaPos[vrt5],
						7.0/24, aaPos[vrt3]);

	//	6
		VecScaleAppend(	aaTmpPos[vrt6],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt5]);

		VecScaleAppend(	aaTmpPos[vrt6],
						3.0/8,  aaPos[vrt6],
						7.0/24, aaPos[vrt1]);
	}

//	Calculate cell valencies of each vertex (:= of associated volumes)
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	//for(VertexIterator vrtIter = mg.begin<Vertex>(i); vrtIter != mg.end<Vertex>(i); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		int num_associated_volumes = 0;

	//	Calculate number of associated volumes
		for(Grid::AssociatedVolumeIterator vIter = mg.associated_volumes_begin(vrt); vIter != mg.associated_volumes_end(vrt); ++vIter)
			num_associated_volumes++;

		//UG_LOG(aaTmpPos[vrt].x() << "; " << aaTmpPos[vrt].y() << "; " << aaTmpPos[vrt].z() << "; " << num_associated_volumes << endl);

		VecScale(aaTmpPos[vrt], aaTmpPos[vrt], 1.0/num_associated_volumes);
	}

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	//for(VertexIterator vrtIter = mg.begin<Vertex>(i); vrtIter != mg.end<Vertex>(i); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		VecScale(aaPos[vrt], aaTmpPos[vrt], 1.0);
	}

//	Export grid
	//SaveGridToUGX(mg, sh, "test.ugx");
	//SaveGridToUGX(mg, sh, "test.ugx");
}

}//	end of namespace

#endif
