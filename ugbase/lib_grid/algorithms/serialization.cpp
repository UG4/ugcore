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
	GOID_GEOMETRIC_OBJECT = 0,
	GOID_VERTEX_BASE = 10,
	GOID_VERTEX = 11,
	GOID_HANGING_VERTEX = 12,
	GOID_EDGE_BASE = 20,
	GOID_EDGE = 21,
	GOID_CONSTRAINED_EDGE = 22,
	GOID_CONSTRAINING_EDGE = 23,
	GOID_FACE = 30,
	GOID_TRIANGLE = 31,
	GOID_CONSTRAINED_TRIANGLE = 32,
	GOID_CONSTRAINING_TRIANGLE = 33,
	GOID_QUADRILATERAL = 40,
	GOID_CONSTRAINED_QUADRILATERAL = 41,
	GOID_CONSTRAINING_QUADRILATERAL = 42,
	GOID_VOLUME = 60,
	GOID_TETRAHEDRON = 61,
	GOID_HEXAHEDRON = 70,
	GOID_PRISM = 80,
	GOID_PYRAMID = 90,
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

//	volumes
	{
		if(goc.num<Tetrahedron>() > 0)
		{
			tInt = GOID_TETRAHEDRON;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Tetrahedron>();
			out.write((char*)&tInt, sizeof(int));
			
			for(TetrahedronIterator iter = goc.begin<Tetrahedron>();
				iter != goc.end<Tetrahedron>(); ++iter)
			{
				Tetrahedron* t = *iter;
				out.write((char*)&aaIntVRT[t->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(3)], sizeof(int));
			}
		}
		
		if(goc.num<Hexahedron>() > 0)
		{
			tInt = GOID_HEXAHEDRON;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Hexahedron>();
			out.write((char*)&tInt, sizeof(int));
			
			for(HexahedronIterator iter = goc.begin<Hexahedron>();
				iter != goc.end<Hexahedron>(); ++iter)
			{
				Hexahedron* h = *iter;
				out.write((char*)&aaIntVRT[h->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(3)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(4)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(5)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(6)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(7)], sizeof(int));
			}
		}
		
		if(goc.num<Prism>() > 0)
		{
			tInt = GOID_PRISM;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Prism>();
			out.write((char*)&tInt, sizeof(int));
			
			for(PrismIterator iter = goc.begin<Prism>();
				iter != goc.end<Prism>(); ++iter)
			{
				Prism* p = *iter;
				out.write((char*)&aaIntVRT[p->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(3)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(4)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(5)], sizeof(int));
			}
		}
		
		if(goc.num<Pyramid>() > 0)
		{
			tInt = GOID_PYRAMID;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<Pyramid>();
			out.write((char*)&tInt, sizeof(int));
			
			for(PyramidIterator iter = goc.begin<Pyramid>();
				iter != goc.end<Pyramid>(); ++iter)
			{
				Pyramid* p = *iter;
				out.write((char*)&aaIntVRT[p->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(3)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(4)], sizeof(int));
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
				case GOID_TETRAHEDRON:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));

							grid.create<Tetrahedron>(TetrahedronDescriptor(
													vVrts[i1], vVrts[i2],
													vVrts[i3], vVrts[i4]));
						}
					}break;
				case GOID_HEXAHEDRON:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5, i6, i7, i8;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							in.read((char*)&i6, sizeof(int));
							in.read((char*)&i7, sizeof(int));
							in.read((char*)&i8, sizeof(int));

							grid.create<Hexahedron>(HexahedronDescriptor(
													vVrts[i1], vVrts[i2],
													vVrts[i3], vVrts[i4],
													vVrts[i5], vVrts[i6],
													vVrts[i7], vVrts[i8]));
						}
					}break;
				case GOID_PRISM:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5, i6;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							in.read((char*)&i6, sizeof(int));

							grid.create<Prism>(PrismDescriptor(
													vVrts[i1], vVrts[i2],
													vVrts[i3], vVrts[i4],
													vVrts[i5], vVrts[i6]));
						}
					}break;
				case GOID_PYRAMID:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));

							grid.create<Pyramid>(PyramidDescriptor(
													vVrts[i1], vVrts[i2],
													vVrts[i3], vVrts[i4],
													vVrts[i5]));
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



////////////////////////////////////////////////////////////////////////
//	writes the parent of the given element - with type and index
template<class TElem>
static void WriteParent(MultiGrid& mg, TElem* pElem,
						Grid::VertexAttachmentAccessor<AInt>& aaIntVRT,
						Grid::EdgeAttachmentAccessor<AInt>& aaIntEDGE,
						Grid::FaceAttachmentAccessor<AInt>& aaIntFACE,
						Grid::VolumeAttachmentAccessor<AInt>& aaIntVOL,
						std::ostream& out)
{
	char type;
	int index;
	GeometricObject* pParent = mg.get_parent(pElem);
	
	if(pParent)
	{
		int parentType = pParent->base_object_type_id();
		
		switch(parentType)
		{
			case VERTEX:
				type = GOID_VERTEX_BASE;
				index = aaIntVRT[(VertexBase*)pParent];
				out.write((char*)&type, sizeof(char));
				out.write((char*)&index, sizeof(int));
				return;
			case EDGE:
				type = GOID_EDGE_BASE;
				index = aaIntEDGE[(EdgeBase*)pParent];
				out.write((char*)&type, sizeof(char));
				out.write((char*)&index, sizeof(int));
				return;
			case FACE:
				type = GOID_FACE;
				index = aaIntFACE[(Face*)pParent];
				out.write((char*)&type, sizeof(char));
				out.write((char*)&index, sizeof(int));
				return;
			case VOLUME:
				type = GOID_VOLUME;
				index = aaIntVOL[(Volume*)pParent];
				out.write((char*)&type, sizeof(char));
				out.write((char*)&index, sizeof(int));
				return;
		}
	}

//	if we reach this point the parent is invalid.
	type = GOID_INVALID;
	index = -1;
	out.write((char*)&type, sizeof(char));
	out.write((char*)&index, sizeof(int));

}

////////////////////////////////////////////////////////////////////////
bool SerializeMultiGridElements(MultiGrid& mg,
								MultiLevelGeometricObjectCollection mgoc,
								AInt& aIntVRT, AInt& aIntEDGE,
								AInt& aIntFACE, AInt& aIntVOL,
								std::ostream& out)
{
//TODO: add volume support
	assert(mg.has_vertex_attachment(aIntVRT) && "aIntVRT is not attached to the grid");
	if(!mg.has_vertex_attachment(aIntVRT))
		return false;
	assert(mg.has_vertex_attachment(aIntEDGE) && "aIntEDGE is not attached to the grid");
	if(!mg.has_vertex_attachment(aIntEDGE))
		return false;
	assert(mg.has_vertex_attachment(aIntFACE) && "aIntFACE is not attached to the grid");
	if(!mg.has_vertex_attachment(aIntFACE))
		return false;
	assert(mg.has_vertex_attachment(aIntVOL) && "aIntVOL is not attached to the grid");
	if(!mg.has_vertex_attachment(aIntVOL))
		return false;

	Grid::VertexAttachmentAccessor<AInt> aaIntVRT(mg, aIntVRT);
	Grid::EdgeAttachmentAccessor<AInt> aaIntEDGE(mg, aIntEDGE);
	Grid::FaceAttachmentAccessor<AInt> aaIntFACE(mg, aIntFACE);
	Grid::VolumeAttachmentAccessor<AInt> aaIntVOL(mg, aIntVOL);

	int tInt;
	number tNumber;

//	iterate through the different levels
	uint numLevels = mgoc.num_levels();
	int vrtInd = 0;
	int edgeInd = 0;
	int faceInd = 0;
	int volInd = 0;
	
	for(uint iLevel = 0; iLevel < numLevels; ++iLevel)
	{
	//	prepare vertices and set num-vertices and num-hanging-vertices.
		{		
		//	write vertices
			if(mgoc.num<Vertex>(iLevel) > 0)
			{
				tInt = GOID_VERTEX;
				out.write((char*)&tInt, sizeof(int));
				tInt = (int)mgoc.num<Vertex>(iLevel);
				out.write((char*)&tInt, sizeof(int));
				
				for(VertexIterator iter = mgoc.begin<Vertex>(iLevel);
					iter != mgoc.end<Vertex>(iLevel); ++iter)
				{
					aaIntVRT[*iter] = vrtInd++;
				//	write the parent
					WriteParent(mg, *iter, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
				}
				
			}
			
		//	write hanging vertices
			if(mgoc.num<HangingVertex>(iLevel) > 0)
			{
				tInt = GOID_HANGING_VERTEX;
				out.write((char*)&tInt, sizeof(int));
				tInt = mgoc.num<HangingVertex>(iLevel);
				out.write((char*)&tInt, sizeof(int));
				
			//	write local-coords and assign indices
			//	we need a number stream for that
				for(HangingVertexIterator iter = mgoc.begin<HangingVertex>(iLevel);
					iter != mgoc.end<HangingVertex>(iLevel); ++iter)
				{
					tNumber = (*iter)->get_local_coordinate_1();
					out.write((char*)&tNumber, sizeof(number));
					tNumber = (*iter)->get_local_coordinate_2();
					out.write((char*)&tNumber, sizeof(number));
					aaIntVRT[*iter] = vrtInd++;
				//	write the parent
					WriteParent(mg, *iter, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
				}
			}
		}

	//	iterate through the edges and set up the edgeStream.
	//int EDGE_GOID, int vrtInd1, int vrtInd2, [int numConstrainedVertices, {int constrainedVertexIndex}, int numConstrainedEdges, {int constrainedEdgeIndex}]
		{
		//	fill the stream
		//	normal edges first.
			if(mgoc.num<Edge>(iLevel) > 0)
			{
				tInt = GOID_EDGE;
				out.write((char*)&tInt, sizeof(int));
				tInt = (int)mgoc.num<Edge>(iLevel);
				out.write((char*)&tInt, sizeof(int));

				for(EdgeIterator iter = mgoc.begin<Edge>(iLevel);
					iter != mgoc.end<Edge>(iLevel); ++iter)
				{
					Edge* e = *iter;
					out.write((char*)&aaIntVRT[e->vertex(0)], sizeof(int));
					out.write((char*)&aaIntVRT[e->vertex(1)], sizeof(int));
					aaIntEDGE[*iter] = edgeInd++;
				//	write the parent
					WriteParent(mg, e, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
				}
			}
		//TODO: add support for hanging edges.
		}

	//	faces
		{
		//TODO: add support for constrained faces etc...
			if(mgoc.num<Triangle>(iLevel) > 0)
			{
				tInt = GOID_TRIANGLE;
				out.write((char*)&tInt, sizeof(int));
				tInt = (int)mgoc.num<Triangle>(iLevel);
				out.write((char*)&tInt, sizeof(int));
				
				for(TriangleIterator iter = mgoc.begin<Triangle>(iLevel);
					iter != mgoc.end<Triangle>(iLevel); ++iter)
				{
					Triangle* t = *iter;
					out.write((char*)&aaIntVRT[t->vertex(0)], sizeof(int));
					out.write((char*)&aaIntVRT[t->vertex(1)], sizeof(int));
					out.write((char*)&aaIntVRT[t->vertex(2)], sizeof(int));
					aaIntFACE[*iter] = faceInd++;
				//	write the parent
					WriteParent(mg, t, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
				}
			}
			
			if(mgoc.num<Quadrilateral>(iLevel) > 0)
			{
				tInt = GOID_QUADRILATERAL;
				out.write((char*)&tInt, sizeof(int));
				tInt = (int)mgoc.num<Quadrilateral>(iLevel);
				out.write((char*)&tInt, sizeof(int));

				for(QuadrilateralIterator iter = mgoc.begin<Quadrilateral>(iLevel);
					iter != mgoc.end<Quadrilateral>(iLevel); ++iter)
				{
					Quadrilateral* q = *iter;
					out.write((char*)&aaIntVRT[q->vertex(0)], sizeof(int));
					out.write((char*)&aaIntVRT[q->vertex(1)], sizeof(int));
					out.write((char*)&aaIntVRT[q->vertex(2)], sizeof(int));
					out.write((char*)&aaIntVRT[q->vertex(3)], sizeof(int));
					aaIntFACE[*iter] = faceInd++;
				//	write the parent
					WriteParent(mg, q, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
				}
			}
		}
	
	//	volumes
		if(mgoc.num<Tetrahedron>(iLevel) > 0)
		{
			tInt = GOID_TETRAHEDRON;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<Tetrahedron>(iLevel);
			out.write((char*)&tInt, sizeof(int));
			
			for(TetrahedronIterator iter = mgoc.begin<Tetrahedron>(iLevel);
				iter != mgoc.end<Tetrahedron>(iLevel); ++iter)
			{
				Tetrahedron* t = *iter;
				out.write((char*)&aaIntVRT[t->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[t->vertex(3)], sizeof(int));
				aaIntVOL[*iter] = volInd++;
			//	write the parent
				WriteParent(mg, t, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
			}
		}
		
		if(mgoc.num<Hexahedron>(iLevel) > 0)
		{
			tInt = GOID_HEXAHEDRON;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<Hexahedron>(iLevel);
			out.write((char*)&tInt, sizeof(int));
			
			for(HexahedronIterator iter = mgoc.begin<Hexahedron>(iLevel);
				iter != mgoc.end<Hexahedron>(iLevel); ++iter)
			{
				Hexahedron* h = *iter;
				out.write((char*)&aaIntVRT[h->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(3)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(4)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(5)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(6)], sizeof(int));
				out.write((char*)&aaIntVRT[h->vertex(7)], sizeof(int));
				aaIntVOL[*iter] = volInd++;
			//	write the parent
				WriteParent(mg, h, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
			}
		}
		
		if(mgoc.num<Prism>(iLevel) > 0)
		{
			tInt = GOID_PRISM;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<Prism>(iLevel);
			out.write((char*)&tInt, sizeof(int));
			
			for(PrismIterator iter = mgoc.begin<Prism>(iLevel);
				iter != mgoc.end<Prism>(iLevel); ++iter)
			{
				Prism* p = *iter;
				out.write((char*)&aaIntVRT[p->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(3)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(4)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(5)], sizeof(int));
				aaIntVOL[*iter] = volInd++;
			//	write the parent
				WriteParent(mg, p, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
			}
		}
		
		if(mgoc.num<Pyramid>(iLevel) > 0)
		{
			tInt = GOID_PYRAMID;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<Pyramid>(iLevel);
			out.write((char*)&tInt, sizeof(int));
			
			for(PyramidIterator iter = mgoc.begin<Pyramid>(iLevel);
				iter != mgoc.end<Pyramid>(iLevel); ++iter)
			{
				Pyramid* p = *iter;
				out.write((char*)&aaIntVRT[p->vertex(0)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(1)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(2)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(3)], sizeof(int));
				out.write((char*)&aaIntVRT[p->vertex(4)], sizeof(int));
				aaIntVOL[*iter] = volInd++;
			//	write the parent
				WriteParent(mg, p, aaIntVRT, aaIntEDGE, aaIntFACE, aaIntVOL, out);
			}
		}
	}
	
//	mark the end of the grid-section
	tInt = GOID_END_OF_GRID;
	out.write((char*)&tInt, sizeof(int));

	return true;
}

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
bool SerializeMultiGridElements(MultiGrid& mg,
								MultiLevelGeometricObjectCollection mgoc,
								std::ostream& out)
{
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);
	
	bool retVal = SerializeMultiGridElements(mg,
						mg.get_multi_level_geometric_object_collection(),
						aInt, aInt, aInt, aInt,
						out);
						
	mg.detach_from_vertices(aInt);
	mg.detach_from_vertices(aInt);
	mg.detach_from_vertices(aInt);
	mg.detach_from_vertices(aInt);
	
	return retVal;
}

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
bool SerializeMultiGridElements(MultiGrid& mg,
								std::ostream& out)
{
	return SerializeMultiGridElements(mg,
						mg.get_multi_level_geometric_object_collection(),
						out);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
static GeometricObject*
GetParent(std::istream& in, vector<VertexBase*>& vVrts,
		vector<EdgeBase*>& vEdges, vector<Face*>& vFaces,
		vector<Volume*> vVols)
{
	char type;
	int index;
	in.read((char*)&type, sizeof(char));
	in.read((char*)&index, sizeof(int));
	
	switch(type)
	{
		case GOID_VERTEX_BASE:
			assert(index < vVrts.size() && "bad index!");
			return vVrts[index];
		case GOID_EDGE_BASE:
			assert(index < vEdges.size() && "bad index!");
			return vEdges[index];
		case GOID_FACE:
			assert(index < vFaces.size() && "bad index!");
			return vFaces[index];
		case GOID_VOLUME:
			assert(index < vVols.size() && "bad index!");
			return vVols[index];
	}
	
	return NULL;
}

////////////////////////////////////////////////////////////////////////
//	DeserializeMultiGridElements
bool DeserializeMultiGridElements(MultiGrid& mg, std::istream& in,
									std::vector<VertexBase*>* pvVrts,
									std::vector<EdgeBase*>* pvEdges,
									std::vector<Face*>* pvFaces,
									std::vector<Volume*>* pvVols)
{
//TODO: add volume support

//	if the user specified element-vectors, we will use them.
//	if not we'll use our own.
	vector<VertexBase*>	vVrtsTMP;
	vector<EdgeBase*>	vEdgesTMP;
	vector<Face*>		vFacesTMP;
	vector<Volume*>		vVolsTMP;
	
	if(!pvVrts)
		pvVrts = &vVrtsTMP;
	if(!pvEdges)
		pvEdges = &vEdgesTMP;
	if(!pvFaces)
		pvFaces = &vFacesTMP;
	if(!pvVols)
		pvVols = &vVolsTMP;
		
	vector<VertexBase*>& vVrts = *pvVrts;
	vector<EdgeBase*>& vEdges = *pvEdges;
	vector<Face*>& vFaces = *pvFaces;
	vector<Volume*>& vVols = *pvVols;
	
	vVrts.clear();
	vEdges.clear();
	vFaces.clear();
	vVols.clear();
	
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
						{
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							vVrts.push_back(*mg.create<Vertex>(parent));
						}
					}break;
					
				case GOID_HANGING_VERTEX:
					{
					//	create the hanging vertices and assign the local coordinates
						for(int i = 0; i < numElems; ++i)
						{
							number coord1, coord2;
							in.read((char*)&coord1, sizeof(number));
							in.read((char*)&coord2, sizeof(number));
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							HangingVertex* hv = *mg.create<HangingVertex>(parent);
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
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Edge* e = *mg.create<Edge>(EdgeDescriptor(vVrts[i1], vVrts[i2]), parent);
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

							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Triangle* t = *mg.create<Triangle>(TriangleDescriptor(
																	vVrts[i1],
																	vVrts[i2],
																	vVrts[i3]),
																	parent);
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
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Quadrilateral* q = *mg.create<Quadrilateral>(QuadrilateralDescriptor(
																	vVrts[i1],
																	vVrts[i2],
																	vVrts[i3],
																	vVrts[i4]),
																	parent);
							vFaces.push_back(q);
						}
					}break;
				case GOID_TETRAHEDRON:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Tetrahedron* t = *mg.create<Tetrahedron>(TetrahedronDescriptor(
																	vVrts[i1], vVrts[i2],
																	vVrts[i3], vVrts[i4]));
							vVols.push_back(t);
						}
					}break;
				case GOID_HEXAHEDRON:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5, i6, i7, i8;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							in.read((char*)&i6, sizeof(int));
							in.read((char*)&i7, sizeof(int));
							in.read((char*)&i8, sizeof(int));
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Hexahedron* h = *mg.create<Hexahedron>(HexahedronDescriptor(
																	vVrts[i1], vVrts[i2],
																	vVrts[i3], vVrts[i4],
																	vVrts[i5], vVrts[i6],
																	vVrts[i7], vVrts[i8]));
							vVols.push_back(h);
						}
					}break;
				case GOID_PRISM:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5, i6;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							in.read((char*)&i6, sizeof(int));
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Prism* p = *mg.create<Prism>(PrismDescriptor(
															vVrts[i1], vVrts[i2],
															vVrts[i3], vVrts[i4],
															vVrts[i5], vVrts[i6]));
							vVols.push_back(p);
						}
					}break;
				case GOID_PYRAMID:
					{
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							GeometricObject* parent = GetParent(in, vVrts, vEdges, vFaces, vVols);
							Pyramid* p = *mg.create<Pyramid>(PyramidDescriptor(
																vVrts[i1], vVrts[i2],
																vVrts[i3], vVrts[i4],
																vVrts[i5]));
							vVols.push_back(p);
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
