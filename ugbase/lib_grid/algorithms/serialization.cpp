// created by Sebastian Reiter
// y09 m11 d05
// s.b.reiter@googlemail.com

#include <cassert>
#include <vector>
#include <algorithm>
#include "serialization.h"
#include "common/serialization.h"
#include "debug_util.h"
#include "common/util/hash.h"

using namespace std;

#define PROFILE_GRID_SERIALIZATION
#ifdef PROFILE_GRID_SERIALIZATION
	#define SRLZ_PROFILE_FUNC()	PROFILE_FUNC_GROUP("serialization")
	#define SRLZ_PROFILE(name)		PROFILE_BEGIN_GROUP(name, "serialization")
	#define SRLZ_PROFILE_END()		PROFILE_END()
#else
	#define SRLZ_PROFILE_FUNC()
	#define SRLZ_PROFILE(name)
	#define SRLZ_PROFILE_END()
#endif


namespace ug
{

////////////////////////////////////////////////////////////////////////
//	Implementation
void GridDataSerializationHandler::add(SPVertexDataSerializer cb)
{
	m_vrtSerializers.push_back(cb);
}

void GridDataSerializationHandler::add(SPEdgeDataSerializer cb)
{
	m_edgeSerializers.push_back(cb);
}

void GridDataSerializationHandler::add(SPFaceDataSerializer cb)
{
	m_faceSerializers.push_back(cb);
}

void GridDataSerializationHandler::add(SPVolumeDataSerializer cb)
{
	m_volSerializers.push_back(cb);
}

void GridDataSerializationHandler::add(SPGridDataSerializer cb)
{
	m_gridSerializers.push_back(cb);
}

template<class TSerializers>
void GridDataSerializationHandler::
write_info(BinaryBuffer& out, TSerializers& serializers) const
{
	for(size_t i = 0; i < serializers.size(); ++i)
		serializers[i]->write_info(out);
}

template<class TSerializers>
void GridDataSerializationHandler::
read_info(BinaryBuffer& in, TSerializers& serializers)
{
	for(size_t i = 0; i < serializers.size(); ++i)
		serializers[i]->read_info(in);
}

void GridDataSerializationHandler::write_infos(BinaryBuffer& out) const
{
	write_info(out, m_vrtSerializers);
	write_info(out, m_edgeSerializers);
	write_info(out, m_faceSerializers);
	write_info(out, m_volSerializers);
	write_info(out, m_gridSerializers);
}

void GridDataSerializationHandler::read_infos(BinaryBuffer& in)
{
	read_info(in, m_vrtSerializers);
	read_info(in, m_edgeSerializers);
	read_info(in, m_faceSerializers);
	read_info(in, m_volSerializers);
	read_info(in, m_gridSerializers);
}

void GridDataSerializationHandler::
serialize(BinaryBuffer& out, GridObjectCollection goc) const
{
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		serialize(out, goc.begin<VertexBase>(lvl), goc.end<VertexBase>(lvl));
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		serialize(out, goc.begin<EdgeBase>(lvl), goc.end<EdgeBase>(lvl));
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		serialize(out, goc.begin<Face>(lvl), goc.end<Face>(lvl));
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		serialize(out, goc.begin<Volume>(lvl), goc.end<Volume>(lvl));
}

void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, GridObjectCollection goc)
{
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		deserialize(in, goc.begin<VertexBase>(lvl), goc.end<VertexBase>(lvl));
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		deserialize(in, goc.begin<EdgeBase>(lvl), goc.end<EdgeBase>(lvl));
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		deserialize(in, goc.begin<Face>(lvl), goc.end<Face>(lvl));
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl)
		deserialize(in, goc.begin<Volume>(lvl), goc.end<Volume>(lvl));
}

template<class TSerializers>
void GridDataSerializationHandler::
deserialization_starts(TSerializers& serializers)
{
	for(size_t i = 0; i < serializers.size(); ++i)
		serializers[i]->deserialization_starts();
}

void GridDataSerializationHandler::
deserialization_starts()
{
	deserialization_starts(m_vrtSerializers);
	deserialization_starts(m_edgeSerializers);
	deserialization_starts(m_faceSerializers);
	deserialization_starts(m_volSerializers);
	deserialization_starts(m_gridSerializers);
}

template<class TSerializers>
void GridDataSerializationHandler::
deserialization_done(TSerializers& serializers)
{
	for(size_t i = 0; i < serializers.size(); ++i)
		serializers[i]->deserialization_done();
}

void GridDataSerializationHandler::
deserialization_done()
{
	deserialization_done(m_vrtSerializers);
	deserialization_done(m_edgeSerializers);
	deserialization_done(m_faceSerializers);
	deserialization_done(m_volSerializers);
	deserialization_done(m_gridSerializers);
}

////////////////////////////////////////////////////////////////////////
SubsetHandlerSerializer::
SubsetHandlerSerializer(ISubsetHandler& sh) :
	m_sh(sh)
{
}

void SubsetHandlerSerializer::
write_info(BinaryBuffer& out) const
{
//	serialize the subset infos
	Serialize(out, m_sh.num_subsets());
	for(int i = 0; i < m_sh.num_subsets(); ++i){
		SubsetInfo& si = m_sh.subset_info(i);
		Serialize(out, si.name);
		Serialize(out, si.color);
		Serialize(out, si.m_propertyMap);
	}
}

void SubsetHandlerSerializer::
read_info(BinaryBuffer& in)
{
//	deserialize the subset infos
	int num;
	Deserialize(in, num);

	for(int i = 0; i < num; ++i){
		SubsetInfo& si = m_sh.subset_info(i);
		Deserialize(in, si.name);
		Deserialize(in, si.color);
		Deserialize(in, si.m_propertyMap);
	}
}

void SubsetHandlerSerializer::
write_data(BinaryBuffer& out, VertexBase* o) const
{
	Serialize(out, m_sh.get_subset_index(o));
}

void SubsetHandlerSerializer::
write_data(BinaryBuffer& out, EdgeBase* o) const
{
	Serialize(out, m_sh.get_subset_index(o));
}

void SubsetHandlerSerializer::
write_data(BinaryBuffer& out, Face* o) const
{
	Serialize(out, m_sh.get_subset_index(o));
}

void SubsetHandlerSerializer::
write_data(BinaryBuffer& out, Volume* o) const
{
	Serialize(out, m_sh.get_subset_index(o));
}

void SubsetHandlerSerializer::
read_data(BinaryBuffer& in, VertexBase* o)
{
	int si;
	Deserialize(in, si);
	m_sh.assign_subset(o, si);
}

void SubsetHandlerSerializer::
read_data(BinaryBuffer& in, EdgeBase* o)
{
	int si;
	Deserialize(in, si);
	m_sh.assign_subset(o, si);
}

void SubsetHandlerSerializer::
read_data(BinaryBuffer& in, Face* o)
{
	int si;
	Deserialize(in, si);
	m_sh.assign_subset(o, si);
}

void SubsetHandlerSerializer::
read_data(BinaryBuffer& in, Volume* o)
{
	int si;
	Deserialize(in, si);
	m_sh.assign_subset(o, si);
}



////////////////////////////////////////////////////////////////////////
//	enumerations

/**
 * Don't change the constants, since they are used i.e. in external files too.
 * If you want to add constants, do so at the end of the enumeration.
 */
enum GridObjectID
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

	GOID_NEW_LEVEL = 1000
};

////////////////////////////////////////////////////////////////////////
//	GRID HEADER
enum GridHeaderConstants{
	GHC_HEADER_BEGIN = 1,
	GHC_HEADER_END = 2,
	GHC_READ_OPTIONS = 3,
};

enum GridHeaderReadOptions{
	GHRO_READ_DEFAULT = 0,
	GHRO_READ_LEVELS =	1 << 0,
	GHRO_READ_PARENTS =	1 << 1
};

struct GridHeader{
	GridHeader() :
		m_readOptions(GHRO_READ_DEFAULT) {}
	GridHeader(uint readOptions) :
		m_readOptions(readOptions)	{}

	bool contains_option(uint option){
		return (m_readOptions & option) == option;
	}

	uint m_readOptions;
};

static void WriteGridHeader(const GridHeader& gridHeader, BinaryBuffer& out)
{
//	we use a temporary integer
//	the header begins
	int t = GHC_HEADER_BEGIN;
	out.write((char*)&t, sizeof(int));

//	we now write the read options
	t = GHC_READ_OPTIONS;
	out.write((char*)&t, sizeof(int));
	out.write((char*)&gridHeader.m_readOptions, sizeof(uint));

//	the header ends
	t = GHC_HEADER_END;
	out.write((char*)&t, sizeof(int));
}

static bool ReadGridHeader(GridHeader& gridHeader, BinaryBuffer& in)
{
//	initialize the header to its defaults
	gridHeader = GridHeader();

//	make sure that the header begins properly
	int t;
	in.read((char*)&t, sizeof(int));

	if(t != GHC_HEADER_BEGIN)
		return false;

	bool bHeaderOpen = true;
	while(!in.eof() && bHeaderOpen){
	//	read the next symbol
		in.read((char*)&t, sizeof(int));

		switch(t){
			case GHC_READ_OPTIONS:{
				int opt;
				in.read((char*)&opt, sizeof(uint));
				gridHeader.m_readOptions = opt;
			}break;

			case GHC_HEADER_END:
				bHeaderOpen = false;
				break;
		}
	}

	if(bHeaderOpen){
	//	the header was not closed properly
		return false;
	}

	return true;
}


////////////////////////////////////////////////////////////////////////
//	PARENT INFO
///	Stores a tuple (type, index), identifying a parent.
typedef std::pair<byte, int> ParentInfo;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SerializeGridElements
bool SerializeGridElements(Grid& grid, BinaryBuffer& out)
{
//	call SerializeGridElements with the grids goc
	return SerializeGridElements(grid,
							grid.get_grid_objects(),
							out);
}

////////////////////////////////////////////////////////////////////////
//	SerializeGridElements
bool SerializeGridElements(Grid& grid, GridObjectCollection goc,
						   BinaryBuffer& out)
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
bool SerializeGridElements(Grid& grid, GridObjectCollection goc,
						   AInt& aIntVRT, BinaryBuffer& out)
{	
//TODO: add volume support
	assert(grid.has_vertex_attachment(aIntVRT) && "aIntVRT is not attached to the grid");
	if(!grid.has_vertex_attachment(aIntVRT))
		return false;

	Grid::VertexAttachmentAccessor<AInt> aaIntVRT(grid, aIntVRT);

	int tInt;
	number tNumber;

//	first we'll write the grid header.
//	since we're writing a normal grid, we use the standard header.
	WriteGridHeader(GridHeader(), out);

//	prepare vertices and set num-vertices and num-hanging-vertices.
	{
		int vrtInd = 0;
			
	//	init vertex-indices (only for Vertey type. Rest is done later on).
		for(RegularVertexIterator iter = goc.begin<RegularVertex>();
			iter != goc.end<RegularVertex>(); ++iter)
		{
			aaIntVRT[*iter] = vrtInd++;
		}
		
	//	write vertices to the stream
		if(goc.num<RegularVertex>() > 0)
		{
			tInt = GOID_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)goc.num<RegularVertex>();
			out.write((char*)&tInt, sizeof(int));
		}
		
	//	write hanging vertices
		if(goc.num<ConstrainedVertex>() > 0)
		{
			tInt = GOID_HANGING_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = goc.num<ConstrainedVertex>();
			out.write((char*)&tInt, sizeof(int));
			
		//	write local-coords and assign indices
		//	we need a number stream for that
			for(ConstrainedVertexIterator iter = goc.begin<ConstrainedVertex>();
				iter != goc.end<ConstrainedVertex>(); ++iter)
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
bool DeserializeGridElements(Grid& grid, BinaryBuffer& in,
							bool readGridHeader)
{
//TODO: add volume support
	vector<VertexBase*>	vVrts;
	vector<EdgeBase*>	vEdges;
	vector<Face*>		vFaces;
	
	GridHeader gridHeader;
	if(readGridHeader){
		if(!ReadGridHeader(gridHeader, in)){
			UG_LOG("Invalid GridHeader.");
			return false;
		}
	}

	if(gridHeader.contains_option(GHRO_READ_LEVELS)){
		UG_LOG("ERROR in DeserializeGridElements: READ_LEVELS not supported for flat grids.");
		return false;
	}
	if(gridHeader.contains_option(GHRO_READ_PARENTS)){
		UG_LOG("ERROR in DeserializeGridElements: READ_PARENTS not supported for flat grids.");
		return false;
	}

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
							vVrts.push_back(*grid.create<RegularVertex>());
					}break;
					
				case GOID_HANGING_VERTEX:
					{
					//	create the hanging vertices and assign the local coordinates
						for(int i = 0; i < numElems; ++i)
						{
							ConstrainedVertex* hv = *grid.create<ConstrainedVertex>();
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
//	This method relies on the fact, that mg is in marking mode and
//	that all and only parents which have already been written to
//	the stream are marked.
template<class TElem>
static void WriteParent(MultiGrid& mg, TElem* pElem,
						MultiElementAttachmentAccessor<AInt>&	aaInt,
						BinaryBuffer& out)
{
	GridObject* pParent = mg.get_parent(pElem);
	char type = mg.parent_type(pElem);
	int index = -1;

	if(pParent && mg.is_marked(pParent)){
		UG_ASSERT(pParent->base_object_id() == type,
				  "parent->base_object_id() and MultiGrid::parent_type mismatch!");
		index = aaInt[pParent];
	}

	out.write((char*)&type, sizeof(char));
	out.write((char*)&index, sizeof(int));

}


////////////////////////////////////////////////////////////////////////
bool SerializeMultiGridElements(MultiGrid& mg,
								GridObjectCollection mgoc,
								MultiElementAttachmentAccessor<AInt>&	aaInt,
								BinaryBuffer& out,
								MultiElementAttachmentAccessor<AGeomObjID>* paaID)
{
	SRLZ_PROFILE_FUNC();

	int tInt;
	number tNumber;

//	first we'll write the header. we have to enable level- and parent-reads
	WriteGridHeader(GridHeader(GHRO_READ_LEVELS | GHRO_READ_PARENTS), out);

//	iterate through the different levels
	uint numLevels = mgoc.num_levels();
	int vrtInd = 0;
	int edgeInd = 0;
	int faceInd = 0;
	int volInd = 0;
	
//	we have to mark all elements which were already written
	mg.begin_marking();

////////////////////////////////
//	vertices
	for(uint iLevel = 0; iLevel < numLevels; ++iLevel)
	{
	//	write the level
		tInt = GOID_NEW_LEVEL;
		out.write((char*)&tInt, sizeof(int));
		out.write((char*)&iLevel, sizeof(uint));

	//	prepare vertices and set num-vertices and num-hanging-vertices.
	//	write vertices
		if(mgoc.num<RegularVertex>(iLevel) > 0)
		{
			tInt = GOID_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<RegularVertex>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(RegularVertexIterator iter = mgoc.begin<RegularVertex>(iLevel);
				iter != mgoc.end<RegularVertex>(iLevel); ++iter)
			{
				aaInt[*iter] = vrtInd++;
				mg.mark(*iter);
				WriteParent(mg, *iter, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}

	//	write hanging vertices
		if(mgoc.num<ConstrainedVertex>(iLevel) > 0)
		{
			tInt = GOID_HANGING_VERTEX;
			out.write((char*)&tInt, sizeof(int));
			tInt = mgoc.num<ConstrainedVertex>(iLevel);
			out.write((char*)&tInt, sizeof(int));
			
		//	write local-coords and assign indices
		//	we need a number stream for that
			for(ConstrainedVertexIterator iter = mgoc.begin<ConstrainedVertex>(iLevel);
				iter != mgoc.end<ConstrainedVertex>(iLevel); ++iter)
			{
				ConstrainedVertex* v = *iter;
				mg.mark(v);
				tNumber = v->get_local_coordinate_1();
				out.write((char*)&tNumber, sizeof(number));
				tNumber = v->get_local_coordinate_2();
				out.write((char*)&tNumber, sizeof(number));
				aaInt[v] = vrtInd++;

			//	write constraining object
				int type = -1;
				int ind = -1;
				if(GridObject* cobj = v->get_constraining_object()){
					type = cobj->base_object_id();
					if(mg.is_marked(cobj)){
						switch(type){
							case EDGE:
								ind = aaInt[static_cast<EdgeBase*>(cobj)];
								break;
							case FACE:
								ind = aaInt[static_cast<Face*>(cobj)];
								break;
						}
					}
				}

				out.write((char*)&type, sizeof(int));
				out.write((char*)&ind, sizeof(int));
				tInt = v->get_parent_base_object_id();
				out.write((char*)&tInt, sizeof(int));

				UG_ASSERT(v->get_parent_base_object_id() != -1,
						  "Bad constraining element id in constrained vertex encountered:"
						   << ElementDebugInfo(mg, v));

				WriteParent(mg, v, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}

////////////////////////////////
	//	iterate through the edges and set up the edgeStream.
	//int EDGE_GOID, int vrtInd1, int vrtInd2, [int numConstrainedVertices, {int constrainedVertexIndex}, int numConstrainedEdges, {int constrainedEdgeIndex}]

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
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				aaInt[*iter] = edgeInd++;
				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}

	//	now constrained edges
		if(mgoc.num<ConstrainedEdge>(iLevel) > 0)
		{
			tInt = GOID_CONSTRAINED_EDGE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<ConstrainedEdge>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(ConstrainedEdgeIterator iter = mgoc.begin<ConstrainedEdge>(iLevel);
				iter != mgoc.end<ConstrainedEdge>(iLevel); ++iter)
			{
				ConstrainedEdge* e = *iter;
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				aaInt[*iter] = edgeInd++;

			//	write constraining object
				int type = -1;
				int ind = -1;
				if(GridObject* cobj = e->get_constraining_object()){
					if(mg.is_marked(cobj)){
						type = cobj->base_object_id();
						switch(type){
							case EDGE:
								ind = aaInt[static_cast<EdgeBase*>(cobj)];
								break;
							case FACE:
								ind = aaInt[static_cast<Face*>(cobj)];
								break;
						}
					}
				}

				out.write((char*)&type, sizeof(int));
				out.write((char*)&ind, sizeof(int));
				tInt = e->get_parent_base_object_id();
				out.write((char*)&tInt, sizeof(int));

				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}

	//	now constraining edges
		if(mgoc.num<ConstrainingEdge>(iLevel) > 0)
		{
			tInt = GOID_CONSTRAINING_EDGE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<ConstrainingEdge>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(ConstrainingEdgeIterator iter = mgoc.begin<ConstrainingEdge>(iLevel);
				iter != mgoc.end<ConstrainingEdge>(iLevel); ++iter)
			{
				ConstrainingEdge* e = *iter;
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				aaInt[*iter] = edgeInd++;
				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}

////////////////////////////////
	//	faces
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
				mg.mark(t);
				out.write((char*)&aaInt[t->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[t->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[t->vertex(2)], sizeof(int));
				aaInt[*iter] = faceInd++;
				WriteParent(mg, t, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
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
				mg.mark(q);
				out.write((char*)&aaInt[q->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[q->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[q->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[q->vertex(3)], sizeof(int));
				aaInt[*iter] = faceInd++;
				WriteParent(mg, q, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}
	
		if(mgoc.num<ConstrainedTriangle>(iLevel) > 0)
		{
			tInt = GOID_CONSTRAINED_TRIANGLE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<ConstrainedTriangle>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(ConstrainedTriangleIterator iter = mgoc.begin<ConstrainedTriangle>(iLevel);
				iter != mgoc.end<ConstrainedTriangle>(iLevel); ++iter)
			{
				ConstrainedTriangle* e = *iter;
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(2)], sizeof(int));
				aaInt[e] = faceInd++;

			//	write constraining object
				int ind = -1;
				if(GridObject* cobj = e->get_constraining_object()){
					if(mg.is_marked(cobj)){
						UG_ASSERT(cobj->base_object_id() == FACE,
								  "A constrained face can only be constrained by "
								  "a constraining face!");
						if(cobj->base_object_id() == FACE)
							ind = aaInt[static_cast<Face*>(cobj)];
					}
				}

				out.write((char*)&ind, sizeof(int));
				tInt = e->get_parent_base_object_id();
				out.write((char*)&tInt, sizeof(int));

				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[e]);
			}
		}

		if(mgoc.num<ConstrainingTriangle>(iLevel) > 0)
		{
			tInt = GOID_CONSTRAINING_TRIANGLE;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<ConstrainingTriangle>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(ConstrainingTriangleIterator iter = mgoc.begin<ConstrainingTriangle>(iLevel);
				iter != mgoc.end<ConstrainingTriangle>(iLevel); ++iter)
			{
				ConstrainingTriangle* e = *iter;
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(2)], sizeof(int));
				aaInt[e] = faceInd++;
				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[e]);
			}
		}

		if(mgoc.num<ConstrainedQuadrilateral>(iLevel) > 0)
		{
			tInt = GOID_CONSTRAINED_QUADRILATERAL;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<ConstrainedQuadrilateral>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(ConstrainedQuadrilateralIterator iter = mgoc.begin<ConstrainedQuadrilateral>(iLevel);
				iter != mgoc.end<ConstrainedQuadrilateral>(iLevel); ++iter)
			{
				ConstrainedQuadrilateral* e = *iter;
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(3)], sizeof(int));
				aaInt[e] = faceInd++;

			//	write constraining object
				int ind = -1;
				if(GridObject* cobj = e->get_constraining_object()){
					if(mg.is_marked(cobj)){
						UG_ASSERT(cobj->base_object_id() == FACE,
								  "A constrained face can only be constrained by "
								  "a constraining face!");
						if(cobj->base_object_id() == FACE)
							ind = aaInt[static_cast<Face*>(cobj)];
					}
				}

				out.write((char*)&ind, sizeof(int));
				tInt = e->get_parent_base_object_id();
				out.write((char*)&tInt, sizeof(int));

				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[e]);
			}
		}

		if(mgoc.num<ConstrainingQuadrilateral>(iLevel) > 0)
		{
			tInt = GOID_CONSTRAINING_QUADRILATERAL;
			out.write((char*)&tInt, sizeof(int));
			tInt = (int)mgoc.num<ConstrainingQuadrilateral>(iLevel);
			out.write((char*)&tInt, sizeof(int));

			for(ConstrainingQuadrilateralIterator iter = mgoc.begin<ConstrainingQuadrilateral>(iLevel);
				iter != mgoc.end<ConstrainingQuadrilateral>(iLevel); ++iter)
			{
				ConstrainingQuadrilateral* e = *iter;
				mg.mark(e);
				out.write((char*)&aaInt[e->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[e->vertex(3)], sizeof(int));
				aaInt[e] = faceInd++;
				WriteParent(mg, e, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[e]);
			}
		}

////////////////////////////////
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
				mg.mark(t);
				out.write((char*)&aaInt[t->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[t->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[t->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[t->vertex(3)], sizeof(int));
				aaInt[*iter] = volInd++;
				WriteParent(mg, t, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
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
				mg.mark(h);
				out.write((char*)&aaInt[h->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(3)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(4)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(5)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(6)], sizeof(int));
				out.write((char*)&aaInt[h->vertex(7)], sizeof(int));
				aaInt[*iter] = volInd++;
				WriteParent(mg, h, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
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
				mg.mark(p);
				out.write((char*)&aaInt[p->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(3)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(4)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(5)], sizeof(int));
				aaInt[*iter] = volInd++;
				WriteParent(mg, p, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
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
				mg.mark(p);
				out.write((char*)&aaInt[p->vertex(0)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(1)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(2)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(3)], sizeof(int));
				out.write((char*)&aaInt[p->vertex(4)], sizeof(int));
				aaInt[*iter] = volInd++;
				WriteParent(mg, p, aaInt, out);
				if(paaID)	Serialize(out, (*paaID)[*iter]);
			}
		}
	}
	
	mg.end_marking();

//	mark the end of the grid-section
	tInt = GOID_END_OF_GRID;
	out.write((char*)&tInt, sizeof(int));

	return true;
}

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
bool SerializeMultiGridElements(MultiGrid& mg,
								GridObjectCollection goc,
								BinaryBuffer& out)
{
	AInt aInt;
	mg.attach_to_all(aInt);
	MultiElementAttachmentAccessor<AInt> aaInt(mg, aInt);
	
	bool retVal = SerializeMultiGridElements(mg, goc, aaInt, out);
						
	mg.detach_from_all(aInt);
	return retVal;
}

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
bool SerializeMultiGridElements(MultiGrid& mg,
								BinaryBuffer& out)
{
	return SerializeMultiGridElements(mg,
						mg.get_grid_objects(),
						out);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
static pair<GridObject*, char>
GetParent(BinaryBuffer& in, const vector<VertexBase*>& vVrts,
		const vector<EdgeBase*>& vEdges, const vector<Face*>& vFaces,
		const vector<Volume*>& vVols)
{
	char type;
	int index;
	in.read((char*)&type, sizeof(char));
	in.read((char*)&index, sizeof(int));

	if(index >= 0){
		switch(type){
			case VERTEX:
				assert(index < (int)vVrts.size() && "bad index!");
				return  make_pair(vVrts[index], type);
			case EDGE:
				assert(index < (int)vEdges.size() && "bad index!");
				return make_pair(vEdges[index], type);
			case FACE:
				assert(index < (int)vFaces.size() && "bad index!");
				return make_pair(vFaces[index], type);
			case VOLUME:
				assert(index < (int)vVols.size() && "bad index!");
				return make_pair(vVols[index], type);
		}
	}
	
	return make_pair<GridObject*, char>(NULL, type);
}

////////////////////////////////////////////////////////////////////////
//	DeserializeMultiGridElements
bool DeserializeMultiGridElements(MultiGrid& mg, BinaryBuffer& in,
									std::vector<VertexBase*>* pvVrts,
									std::vector<EdgeBase*>* pvEdges,
									std::vector<Face*>* pvFaces,
									std::vector<Volume*>* pvVols,
									MultiElementAttachmentAccessor<AGeomObjID>* paaID)
{
	SRLZ_PROFILE_FUNC();

//todo	A parents global id should be serialized and used to identify a parent
//		if it was not sent along with an element but was already contained on
//		the target process.

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

//	Read the header first
	GridHeader gridHeader;
	if(!ReadGridHeader(gridHeader, in)){
		UG_LOG("Invalid GridHeader.");
		return false;
	}

	if(!gridHeader.contains_option(GHRO_READ_LEVELS)){
		UG_LOG("ERROR in DeserializeMultiGridElements: READ_LEVELS required for MultiGrids.");
		return false;
	}
	if(!gridHeader.contains_option(GHRO_READ_PARENTS)){
		UG_LOG("ERROR in DeserializeMultiGridElements: READ_PARENTS required for MultiGrids.");
		return false;
	}


	GeomObjID id;

	SRLZ_PROFILE(srlz_settingUpHashes);
//	create hashes for existing geometric objects
	Hash<GeomObjID, VertexBase*>	vrtHash((int)(1.1f * (float)mg.num<VertexBase>()));
	Hash<GeomObjID, EdgeBase*>	edgeHash((int)(1.1f * (float)mg.num<EdgeBase>()));
	Hash<GeomObjID, Face*>		faceHash((int)(1.1f * (float)mg.num<Face>()));
	Hash<GeomObjID, Volume*>		volHash((int)(1.1f * (float)mg.num<Volume>()));

	vrtHash.reserve(mg.num<VertexBase>());
	edgeHash.reserve(mg.num<EdgeBase>());
	faceHash.reserve(mg.num<Face>());
	volHash.reserve(mg.num<Volume>());


	if(paaID){
	//	add existing elements to the hashes
		for(VertexBaseIterator iter = mg.begin<VertexBase>();
			iter != mg.end<VertexBase>(); ++iter)
		{vrtHash.insert((*paaID)[*iter], *iter);}

		for(EdgeBaseIterator iter = mg.begin<EdgeBase>();
			iter != mg.end<EdgeBase>(); ++iter)
		{edgeHash.insert((*paaID)[*iter], *iter);}

		for(FaceIterator iter = mg.begin<Face>();
			iter != mg.end<Face>(); ++iter)
		{faceHash.insert((*paaID)[*iter], *iter);}

		for(VolumeIterator iter = mg.begin<Volume>();
			iter != mg.end<Volume>(); ++iter)
		{volHash.insert((*paaID)[*iter], *iter);}
	}
	SRLZ_PROFILE_END();

//	create the vertices and store them in vVrts for later indexing.
	{
		SRLZ_PROFILE(srlz_readingData);
		uint currentLevel = 0;
	//	iterate through the stream and create vertices
		while(!in.eof())
		{
		//	read the goid
			int goid = 0;
			in.read((char*)&goid, sizeof(int));

		//	check whether we reached the end of the grid-description.
			if(goid == GOID_END_OF_GRID)
				break;

			if(goid == GOID_NEW_LEVEL){
			//	read the current level and start at the beginning of the loop
				in.read((char*)&currentLevel, sizeof(uint));
				continue;
			}

		//	we have to read more elements. check how many.
			int numElems = 0;
			in.read((char*)&numElems, sizeof(int));

		//	depending on the goid we'll create new elements.
			switch(goid)
			{
				case GOID_VERTEX:
					{
						SRLZ_PROFILE(srlz_vertices);
						for(int i = 0; i < numElems; ++i)
						{
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								VertexBase* oldVrt;
								if(vrtHash.get_entry(oldVrt, id)){
									assert(dynamic_cast<RegularVertex*>(oldVrt));
									vVrts.push_back(oldVrt);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldVrt)))
										mg.associate_parent(oldVrt, parent);
									continue;
								}
								UG_ASSERT(!(parent && mg.num_children<VertexBase>(parent)),
									  "Parent has a child vertex already. "
									  << "ID of parent: " << (*paaID)[parent]
									  << ", ID of existing child: "
									  	  << (*paaID)[mg.get_child<VertexBase>(parent, 0)]
									  << ", ID of new element: " << id);
							}

							if(parent)
								vVrts.push_back(*mg.create<RegularVertex>(parent));
							else{
								vVrts.push_back(*mg.create<RegularVertex>(currentLevel));
								mg.set_parent_type(vVrts.back(), pInfo.second);
							}

							if(paaID)
								(*paaID)[vVrts.back()] = id;
						}
					}break;
					
				case GOID_HANGING_VERTEX:
					{
						SRLZ_PROFILE(srlz_hangingVertices);
					//	create the hanging vertices and assign the local coordinates
						for(int i = 0; i < numElems; ++i)
						{
							number coord1, coord2;
							in.read((char*)&coord1, sizeof(number));
							in.read((char*)&coord2, sizeof(number));
							int cgType;
							int cgInd;
							in.read((char*)&cgType, sizeof(int));
							in.read((char*)&cgInd, sizeof(int));
							int parentBaseObjId;
							in.read((char*)&parentBaseObjId, sizeof(int));

							UG_ASSERT(parentBaseObjId != -1,
									  "Bad constraining element id in constrained vertex encountered");

							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								VertexBase* oldVrt;
								if(vrtHash.get_entry(oldVrt, id)){
									assert(dynamic_cast<ConstrainedVertex*>(oldVrt));
									vVrts.push_back(oldVrt);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldVrt))){
										mg.associate_parent(oldVrt, parent);
									//	make sure that constrained/constraining relations are fine
										switch(parent->base_object_id()){
											case EDGE:{
												ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(parent);
												UG_ASSERT(cge, "Constraining edge has to be of type ConstrainingEdge");
												cge->add_constrained_object(oldVrt);
												static_cast<ConstrainedVertex*>(oldVrt)->set_constraining_object(cge);
											}break;

											case FACE:{
												ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(parent);
												UG_ASSERT(cgf, "Constraining face has to be of type ConstrainingFace");
												cgf->add_constrained_object(oldVrt);
												static_cast<ConstrainedVertex*>(oldVrt)->set_constraining_object(cgf);
											}break;

											default:
												UG_THROW("Constraining object has to be an edge or a face");
												break;
										}
									}
									continue;
								}
							}

							ConstrainedVertex* hv;
							if(parent)
								hv = *mg.create<ConstrainedVertex>(parent);
							else{
								hv = *mg.create<ConstrainedVertex>(currentLevel);
								mg.set_parent_type(hv, pInfo.second);
							}
							hv->set_local_coordinate_1(coord1);
							hv->set_local_coordinate_2(coord2);
							hv->set_parent_base_object_id(parentBaseObjId);
							vVrts.push_back(hv);
							if(paaID)
								(*paaID)[hv] = id;

							if(cgInd != -1){
								switch(cgType){
								case EDGE:
									assert(cgInd < (int)vEdges.size());
									assert(dynamic_cast<ConstrainingEdge*>(vEdges[cgInd]));
									hv->set_constraining_object(vEdges[cgInd]);
									static_cast<ConstrainingEdge*>(vEdges[cgInd])
													->add_constrained_object(hv);
									break;
								case FACE:
									assert(cgInd < (int)vFaces.size());
									assert(dynamic_cast<ConstrainingFace*>(vFaces[cgInd]));
									hv->set_constraining_object(vFaces[cgInd]);
									static_cast<ConstrainingFace*>(vFaces[cgInd])
													->add_constrained_object(hv);
									break;
								}
							}
						}
					}break;
				case GOID_EDGE:
					{
						SRLZ_PROFILE(srlz_edges);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								EdgeBase* oldEdge;
								if(edgeHash.get_entry(oldEdge, id)){
									assert(dynamic_cast<Edge*>(oldEdge));
									vEdges.push_back(oldEdge);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldEdge)))
										mg.associate_parent(oldEdge, parent);
									continue;
								}
							}
							Edge* e;
							if(parent)
								e = *mg.create<Edge>(
										EdgeDescriptor(vVrts[i1], vVrts[i2]), parent);
							else{
								e = *mg.create<Edge>(
										EdgeDescriptor(vVrts[i1], vVrts[i2]), currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}
							vEdges.push_back(e);
							if(paaID)
								(*paaID)[e] = id;
						}
					}break;
				case GOID_CONSTRAINING_EDGE:
					{
						SRLZ_PROFILE(srlz_constrainingEdges);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								EdgeBase* oldEdge;
								if(edgeHash.get_entry(oldEdge, id)){
									assert(dynamic_cast<ConstrainingEdge*>(oldEdge));
									vEdges.push_back(oldEdge);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldEdge)))
										mg.associate_parent(oldEdge, parent);
									continue;
								}
							}

							ConstrainingEdge* e;
							if(parent)
								e = *mg.create<ConstrainingEdge>(
										EdgeDescriptor(vVrts[i1], vVrts[i2]), parent);
							else{
								e = *mg.create<ConstrainingEdge>(
										EdgeDescriptor(vVrts[i1], vVrts[i2]), currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}
							vEdges.push_back(e);
							if(paaID)
								(*paaID)[e] = id;
						}
					}break;
				case GOID_CONSTRAINED_EDGE:
					{
						SRLZ_PROFILE(srlz_constrainedEdges);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							int cgType;
							int cgInd;
							in.read((char*)&cgType, sizeof(int));
							in.read((char*)&cgInd, sizeof(int));
							int parentBaseObjId;
							in.read((char*)&parentBaseObjId, sizeof(int));

							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								EdgeBase* oldEdge;
								if(edgeHash.get_entry(oldEdge, id)){
									assert(dynamic_cast<ConstrainedEdge*>(oldEdge));
									vEdges.push_back(oldEdge);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldEdge))){
										mg.associate_parent(oldEdge, parent);
									//	make sure that constrained/constraining relations are fine
										switch(parent->base_object_id()){
											case EDGE:{
												ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(parent);
												UG_ASSERT(cge, "Constraining edge has to be of type ConstrainingEdge");
												cge->add_constrained_object(oldEdge);
												static_cast<ConstrainedEdge*>(oldEdge)->set_constraining_object(cge);
											}break;

											case FACE:{
												ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(parent);
												UG_ASSERT(cgf, "Constraining face has to be of type ConstrainingFace");
												cgf->add_constrained_object(oldEdge);
												static_cast<ConstrainedEdge*>(oldEdge)->set_constraining_object(cgf);
											}break;

											default:
												UG_THROW("Constraining object has to be an edge or a face");
												break;
										}

									}
									continue;
								}
							}

							ConstrainedEdge* e;
							if(parent)
								e = *mg.create<ConstrainedEdge>(
										EdgeDescriptor(vVrts[i1], vVrts[i2]), parent);
							else{
								e = *mg.create<ConstrainedEdge>(
										EdgeDescriptor(vVrts[i1], vVrts[i2]), currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}

							e->set_parent_base_object_id(parentBaseObjId);

							vEdges.push_back(e);
							if(paaID)
								(*paaID)[e] = id;

							if(cgInd != -1){
								switch(cgType){
								case EDGE:
									assert(cgInd < (int)vEdges.size());
									assert(dynamic_cast<ConstrainingEdge*>(vEdges[cgInd]));
									e->set_constraining_object(vEdges[cgInd]);
									static_cast<ConstrainingEdge*>(vEdges[cgInd])
													->add_constrained_object(e);
									break;
								case FACE:
									assert(cgInd < (int)vFaces.size());
									assert(dynamic_cast<ConstrainingFace*>(vFaces[cgInd]));
									e->set_constraining_object(vFaces[cgInd]);
									static_cast<ConstrainingFace*>(vFaces[cgInd])
													->add_constrained_object(e);
									break;
								}
							}
						}
					}break;
				case GOID_TRIANGLE:
					{
						SRLZ_PROFILE(srlz_triangles);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Face* oldFace;
								if(faceHash.get_entry(oldFace, id)){
									assert(dynamic_cast<Triangle*>(oldFace));
									vFaces.push_back(oldFace);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldFace)))
										mg.associate_parent(oldFace, parent);
									continue;
								}
							}

							Triangle* t;
							if(parent)
								t = *mg.create<Triangle>(TriangleDescriptor(
															vVrts[i1], vVrts[i2],
															vVrts[i3]), parent);
							else{
								t = *mg.create<Triangle>(TriangleDescriptor(
															vVrts[i1], vVrts[i2],
															vVrts[i3]), currentLevel);
								mg.set_parent_type(t, pInfo.second);
							}
							vFaces.push_back(t);
							if(paaID)
								(*paaID)[t] = id;
						}
					}break;
				case GOID_QUADRILATERAL:
					{
						SRLZ_PROFILE(srlz_quadrilaterals);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Face* oldFace;
								if(faceHash.get_entry(oldFace, id)){
									assert(dynamic_cast<Quadrilateral*>(oldFace));
									vFaces.push_back(oldFace);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldFace)))
										mg.associate_parent(oldFace, parent);
									continue;
								}
							}

							Quadrilateral* q;
							if(parent)
								q = *mg.create<Quadrilateral>(QuadrilateralDescriptor(
															vVrts[i1], vVrts[i2], vVrts[i3],
															vVrts[i4]), parent);
							else{
								q = *mg.create<Quadrilateral>(QuadrilateralDescriptor(
															vVrts[i1], vVrts[i2], vVrts[i3],
															vVrts[i4]), currentLevel);
								mg.set_parent_type(q, pInfo.second);
							}
							vFaces.push_back(q);
							if(paaID)
								(*paaID)[q] = id;
						}
					}break;

				case GOID_CONSTRAINING_TRIANGLE:
					{
						SRLZ_PROFILE(srlz_constrainingTriangles);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Face* oldFace;
								if(faceHash.get_entry(oldFace, id)){
									assert(dynamic_cast<ConstrainingFace*>(oldFace));
									vFaces.push_back(oldFace);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldFace)))
										mg.associate_parent(oldFace, parent);
									continue;
								}
							}

							ConstrainingFace* e;
							if(parent)
								e = *mg.create<ConstrainingTriangle>(
										TriangleDescriptor(vVrts[i1], vVrts[i2], vVrts[i3]),
										parent);
							else{
								e = *mg.create<ConstrainingTriangle>(
										TriangleDescriptor(vVrts[i1], vVrts[i2], vVrts[i3]),
										currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}
							vFaces.push_back(e);
							if(paaID)
								(*paaID)[e] = id;
						}
					}break;

				case GOID_CONSTRAINED_TRIANGLE:
					{
						SRLZ_PROFILE(srlz_constrainedTriangles);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							int cgInd;
							in.read((char*)&cgInd, sizeof(int));
							int parentBaseObjId;
							in.read((char*)&parentBaseObjId, sizeof(int));

							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Face* oldFace;
								if(faceHash.get_entry(oldFace, id)){
									assert(dynamic_cast<ConstrainedFace*>(oldFace));
									vFaces.push_back(oldFace);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldFace))){
										mg.associate_parent(oldFace, parent);
									//	make sure that constrained/constraining relations are fine
										UG_ASSERT(parent->base_object_id() == FACE,
												  "Only faces may constrain faces");
										ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(parent);
										UG_ASSERT(cgf, "Constraining face has to be of type ConstrainingFace");
										cgf->add_constrained_object(oldFace);
										static_cast<ConstrainedFace*>(oldFace)->set_constraining_object(cgf);
									}
									continue;
								}
							}

							ConstrainedFace* e;
							if(parent)
								e = *mg.create<ConstrainedTriangle>(
										TriangleDescriptor(vVrts[i1], vVrts[i2], vVrts[i3]),
										parent);
							else{
								e = *mg.create<ConstrainedTriangle>(
										TriangleDescriptor(vVrts[i1], vVrts[i2], vVrts[i3]),
										currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}

							e->set_parent_base_object_id(parentBaseObjId);

							vFaces.push_back(e);
							if(paaID)
								(*paaID)[e] = id;

							if(cgInd != -1){
								assert(cgInd < (int)vFaces.size());
								assert(dynamic_cast<ConstrainingFace*>(vFaces[cgInd]));
								e->set_constraining_object(vFaces[cgInd]);
								static_cast<ConstrainingFace*>(vFaces[cgInd])
												->add_constrained_object(e);
							}
						}
					}break;

				case GOID_CONSTRAINING_QUADRILATERAL:
					{
						SRLZ_PROFILE(srlz_constrainingQuadrilaterals);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Face* oldFace;
								if(faceHash.get_entry(oldFace, id)){
									assert(dynamic_cast<ConstrainingFace*>(oldFace));
									vFaces.push_back(oldFace);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldFace)))
										mg.associate_parent(oldFace, parent);
									continue;
								}
							}

							ConstrainingFace* e;
							if(parent)
								e = *mg.create<ConstrainingQuadrilateral>(
										QuadrilateralDescriptor(vVrts[i1], vVrts[i2],
																vVrts[i3], vVrts[i4]),
										parent);
							else{
								e = *mg.create<ConstrainingQuadrilateral>(
										QuadrilateralDescriptor(vVrts[i1], vVrts[i2],
																vVrts[i3], vVrts[i4]),
										currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}

							vFaces.push_back(e);
							if(paaID)
								(*paaID)[e] = id;
						}
					}break;

				case GOID_CONSTRAINED_QUADRILATERAL:
					{
						SRLZ_PROFILE(srlz_constrainedQuadrilaterals);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							int cgInd;
							in.read((char*)&cgInd, sizeof(int));
							int parentBaseObjId;
							in.read((char*)&parentBaseObjId, sizeof(int));

							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Face* oldFace;
								if(faceHash.get_entry(oldFace, id)){
									UG_ASSERT(dynamic_cast<ConstrainedFace*>(oldFace),
											"Face should be constrained! gid: " << id);
									vFaces.push_back(oldFace);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldFace))){
										mg.associate_parent(oldFace, parent);
									//	make sure that constrained/constraining relations are fine
										UG_ASSERT(parent->base_object_id() == FACE,
												  "Only faces may constrain faces");
										ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(parent);
										UG_ASSERT(cgf, "Constraining face has to be of type ConstrainingFace");
										cgf->add_constrained_object(oldFace);
										static_cast<ConstrainedFace*>(oldFace)->set_constraining_object(cgf);
									}
									continue;
								}
							}

							ConstrainedFace* e;
							if(parent)
								e = *mg.create<ConstrainedQuadrilateral>(
										QuadrilateralDescriptor(vVrts[i1], vVrts[i2],
																vVrts[i3], vVrts[i4]),
										parent);
							else{
								e = *mg.create<ConstrainedQuadrilateral>(
										QuadrilateralDescriptor(vVrts[i1], vVrts[i2],
																vVrts[i3], vVrts[i4]),
										currentLevel);
								mg.set_parent_type(e, pInfo.second);
							}

							e->set_parent_base_object_id(parentBaseObjId);

							vFaces.push_back(e);
							if(paaID)
								(*paaID)[e] = id;

							if(cgInd != -1){
								assert(cgInd < (int)vFaces.size());
								assert(dynamic_cast<ConstrainingFace*>(vFaces[cgInd]));
								e->set_constraining_object(vFaces[cgInd]);
								static_cast<ConstrainingFace*>(vFaces[cgInd])
												->add_constrained_object(e);
							}
						}
					}break;

				case GOID_TETRAHEDRON:
					{
						SRLZ_PROFILE(srlz_tetrahedrons);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Volume* oldVol;
								if(volHash.get_entry(oldVol, id)){
									assert(dynamic_cast<Tetrahedron*>(oldVol));
									vVols.push_back(oldVol);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldVol)))
										mg.associate_parent(oldVol, parent);
									continue;
								}
							}

							Tetrahedron* t;
							if(parent)
								t = *mg.create<Tetrahedron>(TetrahedronDescriptor(
															vVrts[i1], vVrts[i2],
															vVrts[i3], vVrts[i4]),
															parent);
							else{
								t = *mg.create<Tetrahedron>(TetrahedronDescriptor(
															vVrts[i1], vVrts[i2],
															vVrts[i3], vVrts[i4]),
															currentLevel);
								mg.set_parent_type(t, pInfo.second);
							}
							vVols.push_back(t);
							if(paaID)
								(*paaID)[t] = id;
						}
					}break;
				case GOID_HEXAHEDRON:
					{
						SRLZ_PROFILE(srlz_hexahedrons);
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
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Volume* oldVol;
								if(volHash.get_entry(oldVol, id)){
									assert(dynamic_cast<Hexahedron*>(oldVol));
									vVols.push_back(oldVol);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldVol)))
										mg.associate_parent(oldVol, parent);
									continue;
								}
							}

							Hexahedron* h;
							if(parent)
								h = *mg.create<Hexahedron>(HexahedronDescriptor(
													vVrts[i1], vVrts[i2], vVrts[i3],
													vVrts[i4], vVrts[i5], vVrts[i6],
													vVrts[i7], vVrts[i8]), parent);
							else{
								h = *mg.create<Hexahedron>(HexahedronDescriptor(
													vVrts[i1], vVrts[i2], vVrts[i3],
													vVrts[i4], vVrts[i5], vVrts[i6],
													vVrts[i7], vVrts[i8]), currentLevel);
								mg.set_parent_type(h, pInfo.second);
							}
							vVols.push_back(h);
							if(paaID)
								(*paaID)[h] = id;
						}
					}break;
				case GOID_PRISM:
					{
						SRLZ_PROFILE(srlz_prisms);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5, i6;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							in.read((char*)&i6, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Volume* oldVol;
								if(volHash.get_entry(oldVol, id)){
									assert(dynamic_cast<Prism*>(oldVol));
									vVols.push_back(oldVol);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldVol)))
										mg.associate_parent(oldVol, parent);
									continue;
								}
							}

							Prism* p;
							if(parent)
								p = *mg.create<Prism>(PrismDescriptor(
												vVrts[i1], vVrts[i2], vVrts[i3],
												vVrts[i4], vVrts[i5], vVrts[i6]),
												parent);
							else{
								p = *mg.create<Prism>(PrismDescriptor(
												vVrts[i1], vVrts[i2], vVrts[i3],
												vVrts[i4], vVrts[i5], vVrts[i6]),
												currentLevel);
								mg.set_parent_type(p, pInfo.second);
							}
							vVols.push_back(p);
							if(paaID)
								(*paaID)[p] = id;
						}
					}break;
				case GOID_PYRAMID:
					{
						SRLZ_PROFILE(srlz_pyramids);
						for(int i = 0; i < numElems; ++i)
						{
							int i1, i2, i3, i4, i5;
							in.read((char*)&i1, sizeof(int));
							in.read((char*)&i2, sizeof(int));
							in.read((char*)&i3, sizeof(int));
							in.read((char*)&i4, sizeof(int));
							in.read((char*)&i5, sizeof(int));
							pair<GridObject*, char> pInfo = GetParent(in, vVrts, vEdges,
																			vFaces, vVols);
							GridObject* parent = pInfo.first;

							if(paaID){
								Deserialize(in, id);
								Volume* oldVol;
								if(volHash.get_entry(oldVol, id)){
									assert(dynamic_cast<Pyramid*>(oldVol));
									vVols.push_back(oldVol);
								//	make sure that its parent is registered
									if(parent && (!mg.get_parent(oldVol)))
										mg.associate_parent(oldVol, parent);
									continue;
								}
							}

							Pyramid* p;
							if(parent)
								p = *mg.create<Pyramid>(PyramidDescriptor(
													vVrts[i1], vVrts[i2], vVrts[i3],
													vVrts[i4], vVrts[i5]), parent);
							else{
								p = *mg.create<Pyramid>(PyramidDescriptor(
													vVrts[i1], vVrts[i2], vVrts[i3],
													vVrts[i4], vVrts[i5]), currentLevel);
								mg.set_parent_type(p, pInfo.second);
							}
							vVols.push_back(p);
							if(paaID)
								(*paaID)[p] = id;
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
								ISubsetHandler& sh, BinaryBuffer& out)
{
	for(;iterBegin != iterEnd; ++iterBegin)
	{
		int si = sh.get_subset_index(*iterBegin);
		out.write((char*)&si, sizeof(int));
	}
}

////////////////////////////////////////////////////////////////////////
bool SerializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							GridObjectCollection goc,
							BinaryBuffer& out)
{
//	write a magic number at the beginning and at the end.
	int magicNumber = 654664;
	out.write((char*)&magicNumber, sizeof(int));

//	serialize subset-infos
	int numSubsets = (int)sh.num_subsets();
	out.write((char*)&numSubsets, sizeof(int));

	for(int i = 0; i < numSubsets; ++i)
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
	//	write the property map
		Serialize(out, si.m_propertyMap);
	}

	for(size_t i = 0; i < goc.num_levels(); ++i)
	{
	//	serialize vertex-subsets
		WriteSubsetIndicesToStream(goc.begin<VertexBase>(i),
									goc.end<VertexBase>(i),
									sh, out);

	//	serialize edge-subsets
		WriteSubsetIndicesToStream(goc.begin<EdgeBase>(i),
									goc.end<EdgeBase>(i),
									sh, out);

	//	serialize face-subsets
		WriteSubsetIndicesToStream(goc.begin<Face>(i),
									goc.end<Face>(i),
									sh, out);

	//	serialize volume-subsets
		WriteSubsetIndicesToStream(goc.begin<Volume>(i),
									goc.end<Volume>(i),
									sh, out);
	}
	
	out.write((char*)&magicNumber, sizeof(int));

	return true;

}

////////////////////////////////////////////////////////////////////////
//	SerializeSubsetHandler
bool SerializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							BinaryBuffer& out)
{
	return SerializeSubsetHandler(grid, sh,
							grid.get_grid_objects(),
							out);
}

////////////////////////////////////////////////////////////////////////
//	ReadSubsetIndicesFromStream
//	helper method for DeserializeSubsetHandler
template <class TElemIter>
static
void ReadSubsetIndicesFromStream(TElemIter iterBegin, TElemIter iterEnd,
								ISubsetHandler& sh, BinaryBuffer& in)
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
bool DeserializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							GridObjectCollection goc,
							BinaryBuffer& in, bool readPropertyMap)
{
//	read a magic number at the beginning and at the end.
	int magicNumber = 654664;
	int tInd;

//	make sure that the magic number matches
	in.read((char*)&tInd, sizeof(int));
	if(tInd != magicNumber){
		UG_LOG(" WARNING: magic-number mismatch before read in DeserializeSubsetHandler. Data-salad possible!\n");
		return false;
	}

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
		if(nameSize > (int)vBuff.size())
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
	//	read the property map
		if(readPropertyMap)
			Deserialize(in, si.m_propertyMap);
	}

	for(size_t i = 0; i < goc.num_levels(); ++i)
	{
	//	serialize vertex-subsets
		ReadSubsetIndicesFromStream(goc.begin<VertexBase>(i),
									goc.end<VertexBase>(i),
									sh, in);

	//	serialize edge-subsets
		ReadSubsetIndicesFromStream(goc.begin<EdgeBase>(i),
									goc.end<EdgeBase>(i),
									sh, in);

	//	serialize face-subsets
		ReadSubsetIndicesFromStream(goc.begin<Face>(i),
									goc.end<Face>(i),
									sh, in);

	//	serialize volume-subsets
		ReadSubsetIndicesFromStream(goc.begin<Volume>(i),
									goc.end<Volume>(i),
									sh, in);
	}

	//	make sure that the magic number matches
	in.read((char*)&tInd, sizeof(int));
	if(tInd != magicNumber){
				UG_LOG(" WARNING: magic-number mismatch after read in DeserializeSubsetHandler. Data-salad possible!\n");
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
//	DeserializeSubsetHandler
bool DeserializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							BinaryBuffer& in, bool readPropertyMap)
{
	return DeserializeSubsetHandler(grid, sh,
							grid.get_grid_objects(),
							in, readPropertyMap);
}

}//	end of namespace
