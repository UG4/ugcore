// created by Sebastian Reiter
// y09 m11 d05
// s.b.reiter@googlemail.com

#include <cassert>
#include <vector>
#include "serialization.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	enumerations

/**
 * Don't change the constants, since they are used i.e. in external files too.
 * If you want to add constants, do so at the end of the enumeration.
 */
enum GeometricObjectID
{
	GOID_END_OF_GRID = -2,
	GOID_INVALID = -1,
	GOID_VERTEX = 0,
	GOID_HANGING_VERTEX = 1,
	GOID_EDGE = 2,
	GOID_CONSTRAINED_EDGE = 3,
	GOID_CONSTRAINING_EDGE = 4,
	GOID_TRIANGLE = 5,
	GOID_CONSTRAINED_TRIANGLE = 6,
	GOID_CONSTRAINING_TRIANGLE = 7,
	GOID_QUADRILATERAL = 8,
	GOID_CONSTRAINED_QUADRILATERAL = 9,
	GOID_CONSTRAINING_QUADRILATERAL = 10,
	GOID_TETRAHEDRON = 11,
	GOID_HEXAHEDRON = 12,
	GOID_PRISM = 13,
	GOID_PYRAMID = 14,
};

////////////////////////////////////////////////////////////////////////
//	SerializeGridElements
bool SerializeGridElements(Grid& grid, std::ostream& out)
{
//	call SerializeGridElements with the grids goc
	return SerializeGridElements(grid,
							grid.get_geometric_object_collection(),
							out);
}

////////////////////////////////////////////////////////////////////////
//	SerializeGridElements
bool SerializeGridElements(Grid& grid, GeometricObjectCollection goc,
						   std::ostream& out)
{
//	create the required int-attachment and call SerializeGridElements.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	bool retVal = SerializeGridElements(grid, goc, aInt, out);
	grid.detach_from_vertices(aInt);
	return retVal;
}

////////////////////////////////////////////////////////////////////////
//	SerializeGridElements
bool SerializeGridElements(Grid& grid, GeometricObjectCollection goc,
						   AInt& aIntVRT, std::ostream& out)
{
//TODO: add volume support
	assert(grid.has_vertex_attachment(aIntVRT) && "aIntVRT is not attached to the grid");
	if(!grid.has_vertex_attachment(aIntVRT))
		return false;

	Grid::VertexAttachmentAccessor<AInt> aaIntVRT(grid, aIntVRT);

	int tInt;
	number tNumber;

//	prepare vertices and set num-vertices and num-hanging-vertices.
	{
		int vrtInd = 0;
			
	//	init vertex-indices (only for Vertey type. Rest is done later on).
		for(VertexIterator iter = goc.begin<Vertex>();
			iter != goc.end<Vertex>(); ++iter)
		{
			aaIntVRT[*iter] = vrtInd++;
		}
		
	//	write vertices to the stream
		if(goc.num<Vertex>() > 0)
		{
			tInt = GOID_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Vertex>();
			out.write((char*)&tInt, sizeof(int));
		}
		
	//	write hanging vertices
		if(goc.num<HangingVertex>() > 0)
		{
			tInt = GOID_HANGING_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = goc.num<HangingVertex>();
			out.write((char*)&tInt, sizeof(int));
			
		//	write local-coords and assign indices
		//	we need a number stream for that
			for(HangingVertexIterator iter = goc.begin<HangingVertex>();
				iter != goc.end<HangingVertex>(); ++iter)
			{
				tNumber = (*iter)->get_local_coordinate_1();
				out.write((char*)&tNumber, sizeof(number));
				tNumber = (*iter)->get_local_coordinate_2();
				out.write((char*)&tNumber, sizeof(number));
				aaIntVRT[*iter] = vrtInd++;
			}
		}
	}

//	iterate through the edges and set up the edgeStream.
//int EDGE_GOID, int vrtInd1, int vrtInd2, [int numConstrainedVertices, {int constrainedVertexIndex}, int numConstrainedEdges, {int constrainedEdgeIndex}]
	{
		int edgeInd = 0;
		
	//	fill the stream
	//	normal edges first.
		if(goc.num<Edge>() > 0)
		{
			tInt = GOID_EDGE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Edge>();
			out.write((char*)&tInt, sizeof(int));

			for(EdgeIterator iter = goc.begin<Edge>();
				iter != goc.end<Edge>(); ++iter)
			{
				Edge* e = *iter;
				edgeInd++;
				out.write((char*)&aaIntVRT[e->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[e->vertex(1)], sizeof(int));
			}
		}
	//TODO: add support for hanging edges.
	}

//	faces
	{
	//TODO: add support for constrained faces etc...
		if(goc.num<Triangle>() > 0)
		{
			tInt = GOID_TRIANGLE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Triangle>();
			out.write((char*)&tInt, sizeof(int));
			
			for(TriangleIterator iter = goc.begin<Triangle>();
				iter != goc.end<Triangle>(); ++iter)
			{
				Triangle* t = *iter;
				out.write((char*)&aaIntVRT[t->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(2)], sizeof(int));
			}
		}
		
		if(goc.num<Quadrilateral>() > 0)
		{
			tInt = GOID_QUADRILATERAL;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Quadrilateral>();
			out.write((char*)&tInt, sizeof(int));

			for(QuadrilateralIterator iter = goc.begin<Quadrilateral>();
				iter != goc.end<Quadrilateral>(); ++iter)
			{
				Quadrilateral* q = *iter;
				out.write((char*)&aaIntVRT[q->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[q->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[q->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[q->vertex(3)], sizeof(int));
			}
		}
	}

//	mark the end of the grid-section
	tInt = GOID_END_OF_GRID;
	out.write((char*)&tInt, sizeof(int));

	return true;
}

////////////////////////////////////////////////////////////////////////
//	DeserializeGridElements
bool DeserializeGridElements(Grid& grid, std::istream& in)
{
//TODO: add volume support
	vector<VertexBase*>	vVrts;
	vector<EdgeBase*>	vEdges;
	vector<Face*>		vFaces;
	
//	create the vertices and store them in vVrts for later indexing.
	{
	//	iterate through the stream and create vertices
		while(!in.eof())
		{
		//	read the goid
			int goid = 0;
			in.read((char*)&goid, sizeof(int));

		//	check whether we reached the end of the grid-description.
			if(goid == GOID_END_OF_GRID)
				break;
	
		//	we have to read more elements. check how many.
			int numElems = 0;
			in.read((char*)&numElems, sizeof(int));

		//	depending on the goid we'll create new elements.
			switch(goid)
			{
				case GOID_VERTEX:
					{
						for(int i = 0; i < numElems; ++i)
							vVrts.push_back(*grid.create<Vertex>());
					}break;
					
				case GOID_HANGING_VERTEX:
					{
					//	create the hanging vertices and assign the local coordinates
						for(int i = 0; i < numElems; ++i)
						{
							HangingVertex* hv = *grid.create<HangingVertex>();
							number coord1, coord2;
							in.read((char*)&coord1, sizeof(number));
							in.read((char*)&coord2, sizeof(number));
							hv->set_local_coordinate_1(coord1);
							hv->set_local_coordinate_2(coord2);
							vVrts.push_back(hv);
						}
					}break;
				case GOID_EDGE:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));

							Edge* e = *grid.create<Edge>(EdgeDescriptor(vVrts[i1], vVrts[i2]));
							vEdges.push_back(e);
						}
					}break;
				case GOID_TRIANGLE:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));

							Triangle* t = *grid.create<Triangle>(TriangleDescriptor(
																	vVrts[i1],
																	vVrts[i2],
																	vVrts[i3]));
							vFaces.push_back(t);
						}
					}break;
				case GOID_QUADRILATERAL:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));

							Quadrilateral* q = *grid.create<Quadrilateral>(QuadrilateralDescriptor(
																	vVrts[i1],
																	vVrts[i2],
																	vVrts[i3],
																	vVrts[i4]));
							vFaces.push_back(q);
						}
					}break;
				
				default:
					LOG("Unknown geometric-object-id in grid-pack. Aborting reconstruction.\n");
					return false;
			}
		}
	}

	return true;
}
/*
bool SerializeMultiGridElements(MultiGrid& mg,
								GeometricObjectCollection goc,
								AInt& aIntVRT, AInt& aIntEDGE,
								AInt& aIntFACE, aInt& aIntVOL,
								std::ostream& out)
{
//TODO: add volume support
	assert(grid.has_vertex_attachment(aIntVRT) && "aIntVRT is not attached to the grid");
	if(!grid.has_vertex_attachment(aIntVRT))
		return false;
	assert(grid.has_vertex_attachment(aIntEDGE) && "aIntEDGE is not attached to the grid");
	if(!grid.has_vertex_attachment(aIntEDGE))
		return false;
	assert(grid.has_vertex_attachment(aIntFACE) && "aIntFACE is not attached to the grid");
	if(!grid.has_vertex_attachment(aIntFACE))
		return false;
	assert(grid.has_vertex_attachment(aIntVOL) && "aIntVOL is not attached to the grid");
	if(!grid.has_vertex_attachment(aIntVOL))
		return false;
		
	Grid::VertexAttachmentAccessor<AInt> aaIntVRT(grid, aIntVRT);
	Grid::VertexAttachmentAccessor<AInt> aaIntEDGE(grid, aIntEDGE);
	Grid::VertexAttachmentAccessor<AInt> aaIntFACE(grid, aIntFACE);
	Grid::VertexAttachmentAccessor<AInt> aaIntVOL(grid, aIntVOL);

	int tInt;
	number tNumber;

//	write the vertices. Assign numbers
	{
		int vrtInd = 0;
			
	//	init vertex-indices (only for Vertey type. Rest is done later on).
		for(VertexIterator iter = goc.begin<Vertex>();
			iter != goc.end<Vertex>(); ++iter)
		{
			aaIntVRT[*iter] = vrtInd++;
		}
		
	//	write vertices to the stream
		if(goc.num<Vertex>() > 0)
		{
			tInt = GOID_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Vertex>();
			out.write((char*)&tInt, sizeof(int));
		}
		
	//	write hanging vertices
		if(goc.num<HangingVertex>() > 0)
		{
			tInt = GOID_HANGING_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = goc.num<HangingVertex>();
			out.write((char*)&tInt, sizeof(int));
			
		//	write local-coords and assign indices
		//	we need a number stream for that
			for(HangingVertexIterator iter = goc.begin<HangingVertex>();
				iter != goc.end<HangingVertex>(); ++iter)
			{
				tNumber = (*iter)->get_local_coordinate_1();
				out.write((char*)&tNumber, sizeof(number));
				tNumber = (*iter)->get_local_coordinate_2();
				out.write((char*)&tNumber, sizeof(number));
				aaIntVRT[*iter] = vrtInd++;
			}
		}
	}

//	iterate through the edges and set up the edgeStream.
//int EDGE_GOID, int vrtInd1, int vrtInd2, [int numConstrainedVertices, {int constrainedVertexIndex}, int numConstrainedEdges, {int constrainedEdgeIndex}]
	{
		int edgeInd = 0;
		
	//	fill the stream
	//	normal edges first.
		if(goc.num<Edge>() > 0)
		{
			tInt = GOID_EDGE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Edge>();
			out.write((char*)&tInt, sizeof(int));

			for(EdgeIterator iter = goc.begin<Edge>();
				iter != goc.end<Edge>(); ++iter)
			{
				Edge* e = *iter;
				edgeInd++;
				out.write((char*)&aaIntVRT[e->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[e->vertex(1)], sizeof(int));
			}
		}
	//TODO: add support for hanging edges.
	}

//	faces
	{
	//TODO: add support for constrained faces etc...
		if(goc.num<Triangle>() > 0)
		{
			tInt = GOID_TRIANGLE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Triangle>();
			out.write((char*)&tInt, sizeof(int));
			
			for(TriangleIterator iter = goc.begin<Triangle>();
				iter != goc.end<Triangle>(); ++iter)
			{
				Triangle* t = *iter;
				out.write((char*)&aaIntVRT[t->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(2)], sizeof(int));
			}
		}
		
		if(goc.num<Quadrilateral>() > 0)
		{
			tInt = GOID_QUADRILATERAL;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Quadrilateral>();
			out.write((char*)&tInt, sizeof(int));

			for(QuadrilateralIterator iter = goc.begin<Quadrilateral>();
				iter != goc.end<Quadrilateral>(); ++iter)
			{
				Quadrilateral* q = *iter;
				out.write((char*)&aaIntVRT[q->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[q->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[q->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[q->vertex(3)], sizeof(int));
			}
		}
	}

//	mark the end of the grid-section
	tInt = GOID_END_OF_GRID;
	out.write((char*)&tInt, sizeof(int));

	return true;
}
*/

////////////////////////////////////////////////////////////////////////
//	WriteSubsetIndicesToStream
//	helper method for SerializeSubsetHandler
template <class TElemIter>
static
void WriteSubsetIndicesToStream(TElemIter iterBegin, TElemIter iterEnd,
								SubsetHandler& sh, std::ostream& out)
{
	for(;iterBegin != iterEnd; ++iterBegin)
	{
		int si = sh.get_subset_index(*iterBegin);
		out.write((char*)&si, sizeof(int));
	}
}

////////////////////////////////////////////////////////////////////////
//	SerializeSubsetHandler
bool SerializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							GeometricObjectCollection goc,
							std::ostream& out)
{
//	write a magic number at the beginning and at the end.
	int magicNumber = 654664;
	out.write((char*)&magicNumber, sizeof(int));

//	serialize subset-infos
	int tInt;
	tInt = (int)sh.num_subsets();
	out.write((char*)&tInt, sizeof(int));

	for(uint i = 0; i < sh.num_subsets(); ++i)
	{
		SubsetInfo& si = sh.subset_info(i);
	//	write the name (first the size, then the rest)
		int nameSize = si.name.size() + 1;
		out.write((char*)&nameSize, sizeof(int));
		out.write(si.name.c_str(), nameSize);

	//	write the material index
		out.write((char*)&si.materialIndex, sizeof(int));
	//	write the color
		out.write((char*)&si.color, sizeof(vector4));
	//	write the subset-state
		out.write((char*)&si.subsetState, sizeof(uint));
	}

//	serialize vertex-subsets
	WriteSubsetIndicesToStream(goc.begin<VertexBase>(),
								goc.end<VertexBase>(),
								sh, out);

//	serialize edge-subsets
	WriteSubsetIndicesToStream(goc.begin<EdgeBase>(),
								goc.end<EdgeBase>(),
								sh, out);

//	serialize face-subsets
	WriteSubsetIndicesToStream(goc.begin<Face>(),
								goc.end<Face>(),
								sh, out);

//	serialize volume-subsets
	WriteSubsetIndicesToStream(goc.begin<Volume>(),
								goc.end<Volume>(),
								sh, out);

	out.write((char*)&magicNumber, sizeof(int));

	return true;
}

////////////////////////////////////////////////////////////////////////
//	SerializeSubsetHandler
bool SerializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							std::ostream& out)
{
	return SerializeSubsetHandler(grid, sh,
							grid.get_geometric_object_collection(),
							out);
}

////////////////////////////////////////////////////////////////////////
//	ReadSubsetIndicesFromStream
//	helper method for DeserializeSubsetHandler
template <class TElemIter>
static
void ReadSubsetIndicesFromStream(TElemIter iterBegin, TElemIter iterEnd,
								SubsetHandler& sh, std::istream& in)
{
	for(;iterBegin != iterEnd; ++iterBegin)
	{
		int si;
		in.read((char*)&si, sizeof(int));
		sh.assign_subset(*iterBegin, si);
	}
}

////////////////////////////////////////////////////////////////////////
//	DeserializeSubsetHandler
bool DeserializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							GeometricObjectCollection goc,
							std::istream& in)
{
//	read a magic number at the beginning and at the end.
	int magicNumber = 654664;
	int tInd;

//	make sure that the magic number matches
	in.read((char*)&tInd, sizeof(int));
	if(tInd != magicNumber)
		return false;

//	deserialize subset-infos
	int numSubsets;
	in.read((char*)&numSubsets, sizeof(int));

//	a buffer to read the name
	vector<char> vBuff(256);
	for(int i = 0; i < numSubsets; ++i)
	{
		SubsetInfo& si = sh.subset_info(i);
	//	read the name (first the size, then the rest)
		int nameSize;
		in.read((char*)&nameSize, sizeof(int));
	//	check whether the buffer has to be resized
		if(nameSize > vBuff.size())
			vBuff.resize(nameSize);
	//	read the name
		in.read(&vBuff.front(), nameSize);
		si.name = &vBuff.front();

	//	read the material index
		in.read((char*)&si.materialIndex, sizeof(int));
	//	read the color
		in.read((char*)&si.color, sizeof(vector4));
	//	read the subset-state
		in.read((char*)&si.subsetState, sizeof(uint));
	}

//	serialize vertex-subsets
	ReadSubsetIndicesFromStream(goc.begin<VertexBase>(),
								goc.end<VertexBase>(),
								sh, in);

//	serialize edge-subsets
	ReadSubsetIndicesFromStream(goc.begin<EdgeBase>(),
								goc.end<EdgeBase>(),
								sh, in);

//	serialize face-subsets
	ReadSubsetIndicesFromStream(goc.begin<Face>(),
								goc.end<Face>(),
								sh, in);

//	serialize volume-subsets
	ReadSubsetIndicesFromStream(goc.begin<Volume>(),
								goc.end<Volume>(),
								sh, in);

	//	make sure that the magic number matches
	in.read((char*)&tInd, sizeof(int));
	if(tInd != magicNumber)
		return false;

	return true;
}

////////////////////////////////////////////////////////////////////////
//	DeserializeSubsetHandler
bool DeserializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							std::istream& in)
{
	return DeserializeSubsetHandler(grid, sh,
							grid.get_geometric_object_collection(),
							in);
}

}//	end of namespace
