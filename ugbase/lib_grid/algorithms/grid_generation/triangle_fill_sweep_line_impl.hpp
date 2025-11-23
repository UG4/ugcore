/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE_IMPL__
#define __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE_IMPL__

namespace ug
{
template <typename TIterator>
bool TriangleFill_SweepLine(Grid& grid, TIterator edgesBegin,
							TIterator edgesEnd, APosition& aPosVRT,
							AInt& aIntVRT,
							SubsetHandler* pSH,
							int newSubsetIndex)
{
	if(!grid.has_vertex_attachment(aPosVRT)){
		UG_LOG("position attachment missing in TriangleFill_SweepLine");
		return false;
	}

	if(!grid.has_vertex_attachment(aIntVRT)){
		UG_LOG("int attachment missing in TriangleFill_SweepLine");
		return false;
	}

	Grid::VertexAttachmentAccessor aaPos(grid, aPosVRT);
	Grid::VertexAttachmentAccessor aaInt(grid, aIntVRT);

//	set up the vertex and edge arrays and build some additional
//	information that is required when it comes to building the triangles.
	std::vector<Vertex*> vrtPtrs;
	std::vector<vector3> vrts;
	std::vector<int>	edges;

	grid.begin_marking();
	for(TIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
		Edge* e = *iter;

		for(size_t i = 0; i < 2; ++i){
			Vertex* vrt = e->vertex(i);
			if(!grid.is_marked(vrt)){
				aaInt[vrt] = (int)vrts.size();
				vrts.push_back(aaPos[vrt]);
				//vrts.push_back(vector2(aaPos[vrt].x(), aaPos[vrt].y()));
				vrtPtrs.push_back(vrt);
				grid.mark(vrt);
			}

			edges.push_back(aaInt[vrt]);
		}
	}
	grid.end_marking();

//	now call the original implementation
	size_t numInitialEdges = edges.size();
	std::vector<int> faces;
	if(TriangleFill_SweepLine(faces, vrts, edges)){

	//	now create the triangles
		for(size_t i = 0; i < faces.size(); i+=3){
			Triangle* tri = *grid.create<Triangle>(TriangleDescriptor(vrtPtrs[faces[i]],
																	  vrtPtrs[faces[i + 1]],
																	  vrtPtrs[faces[i + 2]]));
			if(pSH){
				pSH->assign_subset(tri, newSubsetIndex);
			}
		}

		return true;
	}

	else{
	//DEBUG ONLY
	//	create edges
		for(size_t i = numInitialEdges; i < edges.size(); i+=2){
			Edge* e = *grid.create<RegularEdge>(EdgeDescriptor(vrtPtrs[edges[i]], vrtPtrs[edges[i+1]]));
			if(pSH){
				pSH->assign_subset(e, newSubsetIndex);
			}
		}
	}

	return false;
}

}//	namespace ug

#endif
