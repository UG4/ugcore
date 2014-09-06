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


////////////////////////////////////////////////////////////////////////////////////////////
//	SplitOctahedronToTetrahedrons
void SplitOctahedronToTetrahedrons(Grid& grid, Octahedron* oct, Volume* parentVol, std::vector<Tetrahedron*>& vTetsOut)
{
	UG_LOG(">SplitOctahedronToTetrahedrons" << std::endl);
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	Determine the shortest diagonal to split upon the octahedron
	int bestDiag = 2;

/*
 *	Recall the refinement of a tetrahedron (s. tetrahdron_rules.cpp). A tetrahedron is
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
 */

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
			UG_LOG(">case 0" << std::endl);
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
			UG_LOG(">case 1" << std::endl);
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
			UG_LOG(">case 2" << std::endl);
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
	//grid.erase(oct);
	UG_LOG(">Done." << std::endl);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	ConvertHybridTetOctGridToTetGrid
void ConvertHybridTetOctGridToTetGrid(MultiGrid& mg)
{
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

	if(mg.num_levels() == 1)
		UG_THROW("No octahedral elements to split in base level." << std::endl);

	std::vector<Tetrahedron*> vTetsOut;
	std::vector<Tetrahedron*> vTetChildrenOfOct;
	std::vector<Tetrahedron*> vTetChildrenOfOctParent;


//	>>>>>>>>>>>>>>>>>>>>>
//	Topological splitting
//	>>>>>>>>>>>>>>>>>>>>>

//	Loop over all levels
	for(size_t i = 1; i < mg.num_levels(); ++i)
	{
		UG_LOG(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Topological splitting " << " LEVEL " << i << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl);
		UG_LOG("#Octahedrons  in level " << i << ": " << mg.num<Octahedron>(i) << std::endl);
		UG_LOG("#Tetrahedrons in level " << i << ": " << mg.num<Tetrahedron>(i) << std::endl);

	//	Loop over all octahedrons on each level
		for(VolumeIterator octIter = mg.begin<Octahedron>(i); octIter != mg.end<Octahedron>(i); ++octIter)
		{
			Octahedron* oct 	= dynamic_cast<Octahedron*>(*octIter);
			Volume* parentVol 	= dynamic_cast<Volume*>(mg.get_parent(oct));

			SplitOctahedronToTetrahedrons(mg, oct, parentVol, vTetsOut);
		}
	}
	UG_LOG(std::endl);


//	>>>>>>>>>>>>>>>>>>>>>
//	Re-parenting VOLUMES
//	>>>>>>>>>>>>>>>>>>>>>

//	Loop over all levels starting with toplevel-1 decrementing till level 1
//	and get current oct o, its children tc and its parent's tetrahedral children tp
	for(size_t i = mg.top_level()-1; i > 0; --i)
	{
	    UG_LOG(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Re-parenting VOLUMES : " << " LEVEL " << i << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl);

	//	Loop over all octahedrons o on each level
	    for(VolumeIterator octIter = mg.begin<Octahedron>(i); octIter != mg.end<Octahedron>(i); ++octIter)
	    {
	        vTetChildrenOfOct.clear();
	        vTetChildrenOfOctParent.clear();

	        Volume* oct		 		= *octIter;
	        Volume* parentVol 		= dynamic_cast<Volume*>(mg.get_parent(oct));

		//	Get and collect child tetrahedrons tp of octahedral parent p
		//	Loop over all o's parent's children
			for(size_t j = 0; j < mg.num_children<Volume>(parentVol); ++j)
			{
			//	Store child tetrahedrons tp in vector vTetChildrenOfOctParent
				if(mg.get_child_volume(parentVol, j)->reference_object_id() == ROID_TETRAHEDRON)
				{
					Tetrahedron* tet = dynamic_cast<Tetrahedron*>(mg.get_child_volume(parentVol, j));
					vTetChildrenOfOctParent.push_back(tet);
				}
				else
					continue;
			}

		//	Get and collect child tetrahedrons tc of octahedron o
	        for(size_t j = 0; j < mg.num_children<Volume>(oct); ++j)
	        {
	            if(mg.get_child_volume(oct, j)->reference_object_id() == ROID_TETRAHEDRON)
	            {
	                Tetrahedron* tet = dynamic_cast<Tetrahedron*>(mg.get_child_volume(oct, j));
	                vTetChildrenOfOct.push_back(tet);
	            }
	            else
	                continue;
	        }

	        UG_LOG("vTetChildrenOfOct.size()		: " << vTetChildrenOfOct.size() 				<< std::endl);

	    //	Check containment of o's child tetrahedrons tc in o's parent's child tetrahedrons tp
		//	Iterate over all child tetrahedrons tp of octahedral parent p
			for(size_t j = 0; j < vTetChildrenOfOctParent.size(); ++j)
			{
				Tetrahedron* parentTet = vTetChildrenOfOctParent[j];

			//	Iterate over all remaining possible child tetrahedron candidates tc
				size_t size = vTetChildrenOfOct.size();
				for(size_t k = 0; k < size; ++k)
				{
					Tetrahedron* childTetCandidate = vTetChildrenOfOct[k];
					vector3 center = CalculateCenter(childTetCandidate, aaPos);

				//	If contained, associate family pair tp<->tc
					if(ContainsPoint(dynamic_cast<Volume*>(parentTet), center, aaPos))
					{
						mg.associate_parent(childTetCandidate, parentTet);
						//vTetChildrenOfOct.erase(vTetChildrenOfOct.begin()+k);
					}

				//	Correct vector size for erased element
					//if (size != (size_t)vTetChildrenOfOct.size())
					//{
					//	--k;
					//	size = vTetChildrenOfOct.size();
					//}
				}
			}

//	        for(size_t j = 0; j < vTetChildrenOfOctParent.size(); ++j)
//			{
//				Tetrahedron* parentTet = vTetChildrenOfOctParent[j];
//
//			//	Iterate over all remaining possible child tetrahedron candidates tc
//				for(size_t k = 0; k < vTetChildrenOfOct.size();)
//				{
//					Tetrahedron* childTetCandidate = vTetChildrenOfOct[k];
//					vector3 center = CalculateCenter(childTetCandidate, aaPos);
//
//				//	If contained, associate family pair tp<->tc
//					if(ContainsPoint(dynamic_cast<Volume*>(parentTet), center, aaPos))
//					{
//						mg.associate_parent(childTetCandidate, parentTet);
//						vTetChildrenOfOct.erase(vTetChildrenOfOct.begin()+k);
//					}
//					else
//						 ++k;
//				}
//			}

	        UG_LOG("vTetChildrenOfOctParent.size()		: " << vTetChildrenOfOctParent.size() 			<< std::endl);
	        UG_LOG(std::endl);
	    }
	}
	UG_LOG(std::endl);










//	>>>>>>>>>>>>>>>>>>>>>
//	Re-parenting EDGES, FACES
//	>>>>>>>>>>>>>>>>>>>>>

	std::vector<Edge*> vEdgeChildrenOfOctParent;
	std::vector<Face*> vFaceChildrenOfOctParent;

	std::vector<Vertex*> vVrtChildrenOfOct;
	std::vector<Edge*> vEdgeChildrenOfOct;
	std::vector<Face*> vFaceChildrenOfOct;
	std::vector<Face*> vFaceChildrenOfOctTmp;
	std::vector<Face*> vFaceChildrenOfOctWithEvenVertex;

//	Loop over all levels starting with toplevel-1 decrementing till level 1
//	and get current oct o, its children and its parent's children
	for(size_t i = mg.top_level()-1; i > 0; --i)
	{
		UG_LOG(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Re-parenting VRT, EDGES, FACES : " << " LEVEL " << i << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl);

	//	Loop over all octahedrons o on each level
		for(VolumeIterator octIter = mg.begin<Octahedron>(i); octIter != mg.end<Octahedron>(i); ++octIter)
		{
			vVrtChildrenOfOct.clear();
			vEdgeChildrenOfOct.clear();
			vFaceChildrenOfOct.clear();
			vFaceChildrenOfOctTmp.clear();
			vFaceChildrenOfOctWithEvenVertex.clear();

			vEdgeChildrenOfOctParent.clear();
			vFaceChildrenOfOctParent.clear();

			Volume* oct		 		= *octIter;
			Volume* parentVol 		= dynamic_cast<Volume*>(mg.get_parent(oct));

		//	Get and collect child edges of octahedral parent p
			for(size_t j = 0; j < mg.num_children<Edge>(parentVol); ++j)
			{
			//	Store child edges in vector vEdgeChildrenOfOctParent
				Edge* e = mg.get_child_edge(parentVol, j);
				vEdgeChildrenOfOctParent.push_back(e);
			}

		//	Get and collect child faces of octahedral parent p
			for(size_t j = 0; j < mg.num_children<Face>(parentVol); ++j)
			{
			//	Store child faces in vector vFaceChildrenOfOctParent
				Face* f = mg.get_child_face(parentVol, j);
				vFaceChildrenOfOctParent.push_back(f);
			}

		//	Get and collect child edges of octahedron o
			for(size_t j = 0; j < mg.num_children<Edge>(oct); ++j)
			{
				Edge* e = mg.get_child_edge(oct, j);
				vEdgeChildrenOfOct.push_back(e);
			}

		//	Get and collect child faces (FF) of octahedron o
			for(size_t j = 0; j < mg.num_children<Face>(oct); ++j)
			{
				Face* f = mg.get_child_face(oct, j);

			//	Variant for collecting a subset (12) of child faces (FF) of octahedron o - 4 remaining
				Vertex* v0 = f->vertex(0);
				Vertex* v1 = f->vertex(1);
				Vertex* v2 = f->vertex(2);

				if(	mg.get_parent(v0)->base_object_id() == VERTEX ||
					mg.get_parent(v1)->base_object_id() == VERTEX ||
					mg.get_parent(v2)->base_object_id() == VERTEX)
				{
					vFaceChildrenOfOct.push_back(f);
					vFaceChildrenOfOctWithEvenVertex.push_back(f);
				}
			}

		//	Gather all child faces of octahedron o BUT those 12 (FF) from above
			for(size_t j = 0; j < mg.num_children<Face>(oct); ++j)
			{
				bool bMatch = false;
				Face* f = mg.get_child_face(oct, j);

				for(size_t k = 0; k < vFaceChildrenOfOctWithEvenVertex.size(); ++k)
				{
					Face* FF = vFaceChildrenOfOctWithEvenVertex[k];

					if(f == FF)
					{
						bMatch = true;
					}
				}

				if(!bMatch)
					vFaceChildrenOfOctTmp.push_back(f);
			}

			UG_LOG("vFaceChildrenOfOctWithEvenVertex.size()	: " << vFaceChildrenOfOctWithEvenVertex.size() 				<< std::endl);
			UG_LOG("vFaceChildrenOfOctTmp.size()		: " << vFaceChildrenOfOctTmp.size() 				<< std::endl);
/*
		//	Get and collect the remaining child faces (FF) of octahedron o
		//	Iterate over all child faces but those previously collected 12 (FF) from above
			for(size_t k = 0; k < vFaceChildrenOfOctTmp.size(); ++k)
			{
				bool bFacesShareEdge = false;
				Face* f = vFaceChildrenOfOctTmp[k];

			//	Get edges of f
				EdgeDescriptor f_e0 = f->edge_desc(0);
				EdgeDescriptor f_e1 = f->edge_desc(1);
				EdgeDescriptor f_e2 = f->edge_desc(2);

			//	Iterate over all previously collected child faces (FF) with an even vertex (12)
			//	and additionally collect the remaining 4 child faces (FF) sharing an edge with the above
				for(size_t j = 0; j < vFaceChildrenOfOctWithEvenVertex.size(); ++j)
				{
					bool bFF_hasEvenVrt = false;
					bool bFF_hasInnerVrt = false;
					int sharedEdgeInd = -1;
					Face* FF = vFaceChildrenOfOctWithEvenVertex[j];

				//	Find the one edge, which doensn't contain an even vertex nor an inner vertex.
				//	If such an edge doesn't exist, continue to next face
					for(size_t l = 0; l < FF->num_edges(); ++l)
					{
						EdgeDescriptor e = FF->edge_desc(l);

						if(	mg.get_parent(e.vertex(0))->base_object_id() == VERTEX || mg.get_parent(e.vertex(1))->base_object_id() == VERTEX)
							bFF_hasEvenVrt = true;

						if(	mg.get_parent(e.vertex(0))->reference_object_id() == ROID_OCTAHEDRON || mg.get_parent(e.vertex(1))->reference_object_id() == ROID_OCTAHEDRON)
							bFF_hasInnerVrt = true;

						if(bFF_hasEvenVrt == false && bFF_hasInnerVrt == false)
							sharedEdgeInd = (int)l;
					}

				//	If such an edge was found, check if the possible child face candidate shares this edge
					if(sharedEdgeInd != -1)
					{
						if(	(	f_e0.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e0.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e1.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e1.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e2.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e2.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(1)	)
						)
						{
							bFacesShareEdge = true;
						}
					}
				}

			//	Push back child face candidate, if it shares the specific edge
				if(bFacesShareEdge)
				{
					UG_LOG("Match" << std::endl);
					vFaceChildrenOfOct.push_back(f);
				}
			}
*/


		//	Get and collect the remaining child faces (FF) of octahedron o
		//	Iterate over all previously collected child faces (FF) with an even vertex (12)
		//	and additionally collect the remaining 4 child faces (FF) sharing an edge with the above
			for(size_t j = 0; j < vFaceChildrenOfOctWithEvenVertex.size(); ++j)
			{
				int sharedEdgeInd = -1;
				Face* FF = vFaceChildrenOfOctWithEvenVertex[j];

			//	Find the one edge, which doensn't contain an even vertex nor an inner vertex.
			//	If such an edge doesn't exist, continue to next face
				for(size_t l = 0; l < FF->num_edges(); ++l)
				{
					bool bFF_hasEvenVrt = false;
					bool bFF_hasInnerVrt = false;
					EdgeDescriptor e = FF->edge_desc(l);

					if(	mg.get_parent(e.vertex(0))->base_object_id() == VERTEX || mg.get_parent(e.vertex(1))->base_object_id() == VERTEX)
						bFF_hasEvenVrt = true;

					if(	mg.get_parent(e.vertex(0))->reference_object_id() == ROID_OCTAHEDRON || mg.get_parent(e.vertex(1))->reference_object_id() == ROID_OCTAHEDRON)
						bFF_hasInnerVrt = true;

					if(bFF_hasEvenVrt == false && bFF_hasInnerVrt == false)
						sharedEdgeInd = (int)l;
				}

			//	If such an edge was found, look for a child face candidate which shares this edge
				if(sharedEdgeInd != -1)
				{
					UG_LOG("Found edge" << std::endl);

				//	Iterate over all child faces but those previously collected 12 (FF) from above
					for(size_t k = 0; k < vFaceChildrenOfOctTmp.size(); ++k)
					{
						Face* f = vFaceChildrenOfOctTmp[k];

					//	Get edges of f
						EdgeDescriptor f_e0 = f->edge_desc(0);
						EdgeDescriptor f_e1 = f->edge_desc(1);
						EdgeDescriptor f_e2 = f->edge_desc(2);

					//	Does f share the edge with FF ?
					//	Consider the different possible orientations
						if(	(	f_e0.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e0.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e0.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e0.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e1.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e1.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e1.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e1.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e2.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e2.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(1)	) ||
							(	f_e2.vertex(1) == FF->edge_desc(sharedEdgeInd).vertex(0) && f_e2.vertex(0) == FF->edge_desc(sharedEdgeInd).vertex(1)	)
						)
						{
							UG_LOG("Edge shared" << std::endl);

						// 	If yes, does it contain the inner octahedral vertex
							if(	mg.get_parent(f->vertex(0))->reference_object_id() == ROID_OCTAHEDRON ||
								mg.get_parent(f->vertex(1))->reference_object_id() == ROID_OCTAHEDRON ||
								mg.get_parent(f->vertex(2))->reference_object_id() == ROID_OCTAHEDRON )
							{
							//	Push back child face candidate
								UG_LOG("Match" << std::endl);
								vFaceChildrenOfOct.push_back(f);
							}
						}
					}
				}
				else
					continue;
			}














		//	Get and collect child vertices of octahedron o
			for(size_t j = 0; j < mg.num_children<Vertex>(oct); ++j)
			{
				Vertex* vrt = mg.get_child_vertex(oct);
				vVrtChildrenOfOct.push_back(vrt);
			}

			//UG_LOG("vVrtChildrenOfOct.size()		: " << vVrtChildrenOfOct.size() 				<< std::endl);
			//UG_LOG("vEdgeChildrenOfOct.size()		: " << vEdgeChildrenOfOct.size() 				<< std::endl);
			//UG_LOG("vFaceChildrenOfOct.size()		: " << vFaceChildrenOfOct.size() 				<< std::endl);

		//	Check containment of o's child vertices and edges in o's parent's child edges
		//	Iterate over all child edges of octahedral parent p
			/*
			for(size_t j = 0; j < vEdgeChildrenOfOctParent.size(); ++j)
			{
				Edge* parentEdge = vEdgeChildrenOfOctParent[j];

			//	Iterate over all remaining possible child vertex candidates
				size_t size_v = vVrtChildrenOfOct.size();
				for(size_t k = 0; k < size_v; ++k)
				{
					Vertex* childVrtCandidate = vVrtChildrenOfOct[k];

				//	If contained, associate family pair
					if(DistancePointToLine(aaPos[childVrtCandidate], aaPos[parentEdge->vertex(0)], aaPos[parentEdge->vertex(1)]) < tol)
					{
						mg.associate_parent(childVrtCandidate, parentEdge);
						vVrtChildrenOfOct.erase(vVrtChildrenOfOct.begin()+k);
					}

				//	Correct vector size for erased element
					if (size_v != (size_t)vVrtChildrenOfOct.size())
					{
						--k;
						size_v = vVrtChildrenOfOct.size();
					}
				}

			//	Iterate over all remaining possible child edge candidates
				size_t size_e = vEdgeChildrenOfOct.size();
				for(size_t k = 0; k < size_e; ++k)
				{
					Edge* childEdgeCandidate = vEdgeChildrenOfOct[k];
					vector3 center = CalculateCenter(childEdgeCandidate, aaPos);

				//	If contained, associate family pair
					if(DistancePointToLine(center, aaPos[parentEdge->vertex(0)], aaPos[parentEdge->vertex(1)]) < tol)
					{
						mg.associate_parent(childEdgeCandidate, parentEdge);
						vEdgeChildrenOfOct.erase(vEdgeChildrenOfOct.begin()+k);
					}

				//	Correct vector size for erased element
					if (size_e != (size_t)vEdgeChildrenOfOct.size())
					{
						--k;
						size_e = vEdgeChildrenOfOct.size();
					}
				}
			}
			*/

			UG_LOG(std::endl);
			//UG_LOG("vVrtChildrenOfOct.size()		: " << vVrtChildrenOfOct.size() 				<< std::endl);
			//UG_LOG("vEdgeChildrenOfOct.size()		: " << vEdgeChildrenOfOct.size() 				<< std::endl);
			UG_LOG("vFaceChildrenOfOct.size()		: " << vFaceChildrenOfOct.size() 				<< std::endl);
			UG_LOG("vFaceChildrenOfOctWithEvenVertex.size()	: " << vFaceChildrenOfOctWithEvenVertex.size() 				<< std::endl);
			//UG_LOG("vEdgeChildrenOfOctParent.size()		: " << vEdgeChildrenOfOctParent.size() 			<< std::endl);
			UG_LOG("vFaceChildrenOfOctParent.size()		: " << vFaceChildrenOfOctParent.size() 			<< std::endl);
			UG_LOG(std::endl);

		//	Check containment of o's child faces in o's parent's child faces (FF)
		//	Iterate over all remaining possible child face candidates
			for(size_t k = 0; k < vFaceChildrenOfOct.size(); ++k)
			{
				Face* childFace = vFaceChildrenOfOct[k];
				vector3 center = CalculateCenter(childFace, aaPos);

				number minDist = std::numeric_limits<number>::max();
				int closestFaceInd = -1;

			//	Iterate over all child faces of octahedral parent p
				for(size_t j = 0; j < vFaceChildrenOfOctParent.size(); ++j)
				{
					Face* parentFaceCandidate = vFaceChildrenOfOctParent[j];
					vector3 normal;
					CalculateNormal(normal, parentFaceCandidate, aaPos);
					vector3 closestPointOnTriangleToCenter;
					number bc1Out;
					number bc2Out;

					//UG_LOG("DistancePointToTriangle pre" << std::endl);
					number dist = DistancePointToTriangle(	closestPointOnTriangleToCenter, bc1Out, bc2Out, center,
								  	  	  	  	  	  	  	 aaPos[parentFaceCandidate->vertex(0)],
								  	  	  	  	  	  	  	 aaPos[parentFaceCandidate->vertex(1)],
								  	  	  	  	  	  	  	 aaPos[parentFaceCandidate->vertex(2)], normal);
					//UG_LOG("DistancePointToTriangle post" << std::endl);

				//	Find closest parent face candidate
					if(dist < minDist)
					{
						//UG_LOG("minDist" << std::endl);
						minDist = dist;
						closestFaceInd = (int)j;
					}
				}

//				if(closestFaceInd != -1)
//				{
//					UG_LOG("associate_parent; " << "closestFaceInd: " << closestFaceInd << "; minDist: " << minDist << std::endl);
//					mg.associate_parent(childFace, vFaceChildrenOfOctParent[closestFaceInd]);
//				}
//				else
//				{
//					UG_THROW("Couldn't find face parent to child face: " << ElementDebugInfo(mg, childFace));
//				}
			}

		//	Check containment of o's child edges in o's parent's child edges or faces
		//	Iterate over all remaining possible child face candidates
			/*
			for(size_t k = 0; k < vEdgeChildrenOfOct.size(); ++k)
			{
				Edge* childEdge = vEdgeChildrenOfOct[k];
				vector3 center = CalculateCenter(childEdge, aaPos);

				Vertex* v0 = childEdge->vertex(0);
				Vertex* v1 = childEdge->vertex(1);

				number minDist = std::numeric_limits<number>::max();
				int closestElemInd = -1;

				number dist;

			//	Case, where child edge can only have parent edge
				if(mg.get_parent(v0)->base_object_id() == VERTEX || mg.get_parent(v1)->base_object_id() == VERTEX)
				{
				//	Iterate over all child edges of octahedral parent p
					for(size_t j = 0; j < vEdgeChildrenOfOctParent.size(); ++j)
					{
						Edge* parentEdgeCandidate = vEdgeChildrenOfOctParent[j];
						dist = DistancePointToLine(center, aaPos[parentEdgeCandidate->vertex(0)], aaPos[parentEdgeCandidate->vertex(1)]);

					//	Find closest parent edge candidate
						if(dist < minDist)
						{
							minDist = dist;
							closestElemInd = (int)j;
						}
					}

					if(closestElemInd != -1)
					{
						mg.associate_parent(childEdge, vEdgeChildrenOfOctParent[closestElemInd]);
					}
					else
					{
						UG_THROW("Couldn't find face parent to child edge: " << ElementDebugInfo(mg, childEdge));
					}
				}

			//	Case, where child edge can only have parent face
				else
				{
				//	Iterate over all child faces of octahedral parent p
					for(size_t j = 0; j < vFaceChildrenOfOctParent.size(); ++j)
					{
						Face* parentFaceCandidate = vFaceChildrenOfOctParent[j];
						vector3 normal;
						CalculateNormal(normal, parentFaceCandidate, aaPos);
						vector3 closestPointOnTriangleToCenter;
						number bc1Out;
						number bc2Out;

						dist = DistancePointToTriangle(	closestPointOnTriangleToCenter, bc1Out, bc2Out, center,
																 aaPos[parentFaceCandidate->vertex(0)],
																 aaPos[parentFaceCandidate->vertex(1)],
																 aaPos[parentFaceCandidate->vertex(2)], normal);

					//	Find closest parent face candidate
						if(dist < minDist)
						{
							minDist = dist;
							closestElemInd = (int)j;
						}
					}

					if(closestElemInd != -1)
					{
						mg.associate_parent(childEdge, vFaceChildrenOfOctParent[closestElemInd]);
					}
					else
					{
						UG_THROW("Couldn't find face parent to child edge: " << ElementDebugInfo(mg, childEdge));
					}
				}
			}
			*/

			/*
		//	TEST GATHERING: Get and collect child faces of octahedron o
			vFaceChildrenOfOct.clear();
			for(size_t j = 0; j < mg.num_children<Face>(oct); ++j)
			{
				Face* f = mg.get_child_face(oct, j);
				vFaceChildrenOfOct.push_back(f);
			}

		//	Get and collect child edges of octahedron o
			vEdgeChildrenOfOct.clear();
			for(size_t j = 0; j < mg.num_children<Edge>(oct); ++j)
			{
				Edge* e = mg.get_child_edge(oct, j);
				vEdgeChildrenOfOct.push_back(e);
			}

			UG_LOG("vVrtChildrenOfOct.size()		: " << vVrtChildrenOfOct.size() 				<< std::endl);
			UG_LOG("vEdgeChildrenOfOct.size()		: " << vEdgeChildrenOfOct.size() 				<< std::endl);
			UG_LOG("vFaceChildrenOfOct.size()		: " << vFaceChildrenOfOct.size() 				<< std::endl);
			UG_LOG("vEdgeChildrenOfOctParent.size()		: " << vEdgeChildrenOfOctParent.size() 			<< std::endl);
			UG_LOG("vFaceChildrenOfOctParent.size()		: " << vFaceChildrenOfOctParent.size() 			<< std::endl);
			UG_LOG(std::endl);
			*/
		}
	}
	UG_LOG(std::endl);












/*
//	INNEREI TEST
	UG_LOG("**********************OCT" << std::endl);
	for(size_t i = 0; i < mg.num_levels(); ++i)
	{
		UG_LOG("************************* Level " << i << std::endl);
		for(VolumeIterator vIter = mg.begin<Octahedron>(i); vIter != mg.end<Octahedron>(i); ++vIter)
		{
			UG_LOG("**************************** Octahedron" << std::endl);
			Volume* oct = dynamic_cast<Volume*>(*vIter);
			UG_LOG("mg.num_child_vertices(oct): " << mg.num_child_vertices(oct) << std::endl);
			UG_LOG("mg.num_child_edges(oct): " << mg.num_child_edges(oct) << std::endl);
			UG_LOG("mg.num_child_faces(oct): " << mg.num_child_faces(oct) << std::endl);
		}
	}
	UG_LOG("**********************TET" << std::endl);


//	INNEREI TEST
	UG_LOG("**********************" << std::endl);
	for(size_t i = 0; i < mg.num_levels(); ++i)
	{
		UG_LOG("************************* Level " << i << std::endl);
		for(VolumeIterator vIter = mg.begin<Tetrahedron>(i); vIter != mg.end<Tetrahedron>(i); ++vIter)
		{
			UG_LOG("**************************** Tetrahedron" << std::endl);
			Volume* oct = dynamic_cast<Volume*>(*vIter);
			UG_LOG("mg.num_child_vertices(oct): " << mg.num_child_vertices(oct) << std::endl);
			UG_LOG("mg.num_child_edges(oct): " << mg.num_child_edges(oct) << std::endl);
			UG_LOG("mg.num_child_faces(oct): " << mg.num_child_faces(oct) << std::endl);
		}
	}
	UG_LOG("**********************" << std::endl);
*/


//	Erase all octahedrons o from hierarchy
	mg.erase(mg.begin<Octahedron>(), mg.end<Octahedron>());

/*
//	CHECK
	UG_LOG("----------------------- VOLUME CHECK -----------------------" << std::endl);

	int counter = 0;

	for(size_t i = 0; i < mg.num_levels(); ++i)
	{
		UG_LOG(">>>>>>>>>>>>>>>>>>>> Level " << i << ">>>>>>>>>>>>>>>>>>>>" << std::endl);
		counter = 0;

		for(VolumeIterator vIter = mg.begin<Volume>(i); vIter != mg.end<Volume>(i); ++vIter)
		{
			Volume* vol = *vIter;
			Volume* parent = dynamic_cast<Volume*>(mg.get_parent(vol));

			if(!mg.has_children(vol))
			{
				UG_LOG("Volume " << counter+1 << " " << vol->reference_object_id() << " in level " << i << " without child volume(s).");
			}
			else
			{
				UG_LOG("Volume " << counter+1 << " " << vol->reference_object_id() << " in level " << i << " with "<< mg.num_children<Volume>(vol) << " child volume(s).");
			}

			if(i > 0)
				UG_LOG(" Parent of type: " << parent->reference_object_id());

			UG_LOG(std::endl);

		counter++;
		}
	}


//	CHECK
	UG_LOG("----------------------- EDGE CHECK -----------------------" << std::endl);

	//int counter = 0;
	counter = 0;

	for(size_t i = 0; i < mg.num_levels(); ++i)
	{
		UG_LOG(">>>>>>>>>>>>>>>>>>>> Level " << i << ">>>>>>>>>>>>>>>>>>>>" << std::endl);
		counter = 0;

		for(EdgeIterator eIter = mg.begin<Edge>(i); eIter != mg.end<Edge>(i); ++eIter)
		{
			Edge* e = *eIter;
			GridObject* parent = mg.get_parent(e);
			//Edge* parent = dynamic_cast<Edge*>(mg.get_parent(e));

			if(!mg.has_children(e))
			{
				UG_LOG("Edge " << counter+1 << " " << e->reference_object_id() << " in level " << i << " without child edge(s).");
			}
			else
			{
				UG_LOG("Edge " << counter+1 << " " << e->reference_object_id() << " in level " << i << " with "<< mg.num_children<Edge>(e) << " child edge(s).");
			}

			if(i > 0)
				UG_LOG(" Parent of type: " << parent->reference_object_id());
				//UG_LOG(" Parent of type: " << parent->reference_object_id());

			UG_LOG(std::endl);

		counter++;
		}
	}
*/
}


////////////////////////////////////////////////////////////////////////////////////////////
//	ProjectToLimitSmoothTetGrid
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
//	MoveVertexToSubdivisionLoopPosition
void MoveVertexToSubdivisionLoopPosition(MultiGrid& mg, Vertex* vrt,
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


////////////////////////////////////////////////////////////////////////////////////////////
//	SubdivisionTetGridSmooth
/// (see Schaefer et al, "Smooth subdivision of tetrahedral meshes")
void SubdivisionTetGridSmooth(MultiGrid& mg, bool bPreserveBnd, bool bSubdivisionLoopBnd)
{
//	Ensure correct user input
	if(!bPreserveBnd)
		bSubdivisionLoopBnd = false;

	typedef APosition::ValueType pos_type;

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
		UG_THROW("SubdivisionTetGridSmooth: No volumes included in input grid for smooth TetGrid subdivision refinement.");

// 	Loop all vertices of top_level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		if(!bPreserveBnd)
		{
			MoveVertexToSmoothTetGridSubdivisionPosition(mg, vrt, aaPos, aaSmoothPos);
		}

	//	Surface preservation
		else
		{
		//	BOUNDARY VERTEX
			if(IsBoundaryVertex3D(mg, vrt))
			{
				if(bSubdivisionLoopBnd)
					MoveVertexToSubdivisionLoopPosition(mg, vrt, aaPos, aaSmoothPos);
				else
					aaSmoothPos[vrt] = aaPos[vrt];
			}

		//	INNER VERTEX
			else
			{
				MoveVertexToSmoothTetGridSubdivisionPosition(mg, vrt, aaPos, aaSmoothPos);
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
