//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d01

#include <fstream>
#include <vector>
#include <iostream>
#include <cstring>
#include "file_io_art.h"
#include "../lib_grid.h"
#include "lib_grid/subset_handler.h"

using namespace std;

namespace ug
{

bool LoadGridFromART(Grid& grid, const char* filename,
					 ISubsetHandler* pSH,
					 AVector3& aPos)
{
//	open the stream
	ifstream in(filename);

	if(!in)
		return false;

//	The subset handler:
	SubsetHandler shTmp;
	if(!pSH)
	{
		shTmp.assign_grid(grid);
		pSH = &shTmp;
	}
	ISubsetHandler& sh = *pSH;

//	make sure that the position attachment is attached properly
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);

//	this buffer will be used to store all the lines.
	char* buf = new char[512];
	char* delim = " :";

//	read the vertices
	vector<Vertex*>	vVrts;

	while((!in.eof()))
	{
		in.getline(buf, 512);
		char* tok = strtok(buf, delim);
		if(!tok)
			continue;

	//	ignore all lines that start with a %
		if(*tok == '%')
			continue;
	//	if we found a $ this element-type is done
		if(*tok == '$')
			break;

	//	create a new vertex
		Vertex* v = *grid.create<Vertex>();
	//	read the coordinates
//TODO:	make sure that everything is ok.
		aaPos[v].x = atof(tok);
		tok = strtok(NULL, delim);
		aaPos[v].y = atof(tok);
		tok = strtok(NULL, delim);
		aaPos[v].z = atof(tok);

	//	store the vertex in an array
		vVrts.push_back(v);
	}

//	read the edges
	vector<Edge*>	vEdges;

	while((!in.eof()))
	{
		in.getline(buf, 512);
		char* tok = strtok(buf, delim);
		if(!tok)
			continue;

	//	ignore all lines that start with a %
		if(*tok == '%')
			continue;
	//	if we found a $ this element-type is done
		if(*tok == '$')
			break;

	//	read the indices
//TODO:	make sure that everything is ok.
		int si, i1, i2;
		si = atoi(tok);
		tok = strtok(NULL, delim);
		i1 = atoi(tok);
		tok = strtok(NULL, delim);
		i2 = atoi(tok);

	//	create a new edge
		Edge* e = *grid.create<Edge>(EdgeDescriptor(vVrts[i1], vVrts[i2]));
//	edges won't be assigned to subsets in the moment, since things are a little chaotic in the files...
/*
		if(si < 1)
			si = -1;
		else if(si == 10000)
			si = 0;
		
		if(si != -1)
			sh.assign_subset(e, si);
*/
	//	store the edge in an array
		vEdges.push_back(e);
	}

//	read the faces
	vector<Triangle*>	vTris;

	while((!in.eof()))
	{
		in.getline(buf, 512);
		char* tok = strtok(buf, delim);
		if(!tok)
			continue;

	//	ignore all lines that start with a %
		if(*tok == '%')
			continue;
	//	if we found a $ this element-type is done
		if(*tok == '$')
			break;

	//	read the indices
//TODO:	make sure that everything is ok.
		int si, i1, i2, i3;
		si = atoi(tok);
		tok = strtok(NULL, delim);
		i1 = atoi(tok);
		tok = strtok(NULL, delim);
		i2 = atoi(tok);
		tok = strtok(NULL, delim);
		i3 = atoi(tok);

	//	we have to find the vertices that make up the face.
	//	We can get them from two of the edges.
		Edge* e1 = vEdges[i1];
		Edge* e2 = vEdges[i2];

		VertexBase* vrt[3];
		vrt[0] = e1->vertex(0);
		vrt[1] = e1->vertex(1);
		vrt[2] = e2->vertex(0);//just a guess.
		if((vrt[2] == vrt[0]) || (vrt[2] == vrt[1]))
			vrt[2] = e2->vertex(1);

	//	create a new face
		Triangle* t = *grid.create<Triangle>(TriangleDescriptor(vrt[0], vrt[1], vrt[2]));

		if(si < 1)
			si = -1;
		else if(si == 10000)
			si = 0;
		
		if(si != -1)
			sh.assign_subset(t, si);

	//	store the face in an array
		vTris.push_back(t);
	}

//	read the volumes
	while((!in.eof()))
	{
		in.getline(buf, 512);
		char* tok = strtok(buf, delim);
		if(!tok)
			continue;

	//	ignore all lines that start with a %
		if(*tok == '%')
			continue;
	//	if we found a $ this element-type is done
		if(*tok == '$')
			break;

	//	read the indices
//TODO:	make sure that everything is ok.
		int si, i1, i2, i3, i4;
		si = atoi(tok);
		tok = strtok(NULL, delim);
		i1 = atoi(tok);
		tok = strtok(NULL, delim);
		i2 = atoi(tok);
		tok = strtok(NULL, delim);
		i3 = atoi(tok);
		tok = strtok(NULL, delim);
		i4 = atoi(tok);

	//	we have to find the vertices that make up the volume.
	//	We can get them from two of the faces.
		Triangle* t1 = vTris[i1];
		Triangle* t2 = vTris[i2];

		VertexBase* vrt[4];
		vrt[0] = t1->vertex(0);
		vrt[1] = t1->vertex(1);
		vrt[2] = t1->vertex(2);
		vrt[3] = t2->vertex(0);//just a guess.
		if((vrt[3] == vrt[0]) || (vrt[3] == vrt[1]) || (vrt[3] == vrt[2]))
			vrt[3] = t2->vertex(1);
		if((vrt[3] == vrt[0]) || (vrt[3] == vrt[1]) || (vrt[3] == vrt[2]))
			vrt[3] = t2->vertex(2);

	//	create a new face
		Volume* v = *grid.create<Tetrahedron>(TetrahedronDescriptor(vrt[0], vrt[1], vrt[2], vrt[3]));

		si -= 1;
		sh.assign_subset(v, si);
	}

//	delete the buffer
	delete[] buf;

	return true;
}

/*
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
	vector<VertexBase*>	vVrts;
	vVrts.reserve(numVrts);

	for(int i = 0; i < numVrts; ++i)
		vVrts[i] = *grid.create<Vertex>();

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	read the points
	{
		for(VertexBaseIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter)
		{
			int Index;
			in >> Index;
			in >> aaPos[*iter].x;
			in >> aaPos[*iter].y;
			in >> aaPos[*iter].z;
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

		for(VertexBaseIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++)
		{
			out << counter << " " << 	aaPos[*iter].x << " " <<
										aaPos[*iter].y << " " <<
										aaPos[*iter].z << endl;

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
*/
}//	end of namespace
