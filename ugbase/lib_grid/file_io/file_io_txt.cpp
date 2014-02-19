//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d14

#include <fstream>
#include <vector>
#include "file_io_txt.h"
#include "../lg_base.h"

using namespace std;

namespace ug
{

bool LoadGridFromTXT(Grid& grid, const char* filename, AVector3& aPos)
{
	ifstream in(filename);

	if(!in)
		return false;

	//grid.clear();

	int numVrts, numTris;

	in >> numVrts;
	in >> numTris;

//	create points
//	store pointers to the vertices on the fly in a vector.
	vector<Vertex*>	vVrts;
	vVrts.reserve(numVrts);

	for(int i = 0; i < numVrts; ++i)
		vVrts[i] = *grid.create<RegularVertex>();

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	read the points
	{
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter)
		{
			int Index;
			in >> Index;
			in >> aaPos[*iter].x();
			in >> aaPos[*iter].y();
			in >> aaPos[*iter].z();
		}
	}

//	read the triangles
	{
		for(int i = 0; i < numTris; ++i)
		{
			int Index, i1, i2, i3;
			in >> Index;
			in >> i1;
			in >> i2;
			in >> i3;
			grid.create<Triangle>(TriangleDescriptor(vVrts[i1], vVrts[i2], vVrts[i3]));
		}
	}

	in.close();
	return true;
}

bool SaveGridToTXT(Grid& grid, const char* filename, AVector3& aPos)
{
	if(!grid.has_vertex_attachment(aPos))
		return false;

	ofstream out(filename);

	if(!out)
		return false;

//	write the header
	out << grid.num_vertices() << " " << grid.num<Triangle>() << endl;

//	write the vertices
//	store in each vertex at which position it has been written to the file.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);

	{
		int counter = 0;

		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++)
		{
			out << counter << " " << 	aaPos[*iter].x() << " " <<
										aaPos[*iter].y() << " " <<
										aaPos[*iter].z() << endl;

			aaInt[*iter] = counter++;
		}
	}

//	write the faces.
	{
		int counter = 0;

		for(TriangleIterator iter = grid.begin<Triangle>(); iter != grid.end<Triangle>(); ++iter, ++counter)
		{
			Triangle* t = *iter;
			out << counter << " " 	<< aaInt[t->vertex(0)] << " "
									<< aaInt[t->vertex(1)] << " "
									<< aaInt[t->vertex(2)] << endl;
		}
	}

	grid.detach_from_vertices(aInt);
	return true;
}

}//	end of namespace
