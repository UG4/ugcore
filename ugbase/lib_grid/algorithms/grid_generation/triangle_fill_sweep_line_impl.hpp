//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m11 d19

#ifndef __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE_IMPL__
#define __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE_IMPL__

namespace ug
{
template <class TIterator>
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

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosVRT);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aIntVRT);

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
