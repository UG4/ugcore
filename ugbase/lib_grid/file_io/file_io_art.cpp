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

/**
 * If you havn't already used strtok on the given buffer,
 * pass true to newTokBuff (default == true).
 *
 * Make sure that the buffer only contains index values!
 */
static void ReadIndices(vector<int>& vIndsOut, char* buffer,
						const char* delim, bool newTokBuff = true)
{
	vIndsOut.clear();
	
	char* tok;
	if(newTokBuff)
		tok = strtok(buffer, delim);
	else
		tok = strtok(NULL, delim);
		
	while(tok)
	{
		vIndsOut.push_back(atoi(tok));
		tok = strtok(NULL, delim);
	}
}

template <class TType>

static void CollectUniqueObjects(vector<TType>& vecOut,
								const vector<TType>& vecIn)
{
//	copy indices to avoid problems when vIndsOut and vIndsInd operate on
//	the same storage.
	bool identicContainers = (&vecOut == &vecIn);

	vector<TType> vTmp;
	vector<TType>* pDest = &vecOut;
	if(identicContainers)
		pDest = &vTmp;
	
//	we're operating on destInds
	vector<TType>& vecDest = *pDest;
	vecDest.clear();

//	iterate over all elements of vInds
	for(size_t i = 0; i < vecIn.size(); ++i){
		const TType& val = vecIn[i];
	//	check whether val is already contained in vecDest
		bool gotOne = false;
		for(size_t j = 0; j < vecDest.size(); ++j){
			if(vecDest[j] == val){
				gotOne = true;
				break;
			}
		}
		
	//	if we havn't found one, we'll insert val into vecDest
		if(!gotOne)
			vecDest.push_back(val);
	}
	
//	swap storage if required
	if(identicContainers)
		vecOut.swap(vTmp);
}

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
	const char* delim = " :";

//	those vectors will be reused at different sections
	vector<int> vInds;
	vector<VertexBase*> vVrtDump;
	vector<VertexBase*> vLocalVrts;
	
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
	vector<Face*>	vFaces;

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

		int si = atoi(tok);
		
	//	read the indices
		ReadIndices(vInds, buf, delim, false);

	//	there have to be at least three edges
		if(vInds.size() < 3)
			continue;
			
	//	copy all vertices of the edges to vVrtDump
//TODO:	make sure that all indices are in the correct ranges
		vVrtDump.clear();
		for(size_t i = 0; i < vInds.size(); ++i){
			EdgeBase* e = vEdges[vInds[i]];
			vVrtDump.push_back(e->vertex(0));
			vVrtDump.push_back(e->vertex(1));
		}

	//	now collect the unique vertices in vVrtDump
		CollectUniqueObjects(vLocalVrts, vVrtDump);
		
//TODO: check orientation
	//	create the faces
		Face* f = NULL;
		size_t numVrts = vLocalVrts.size();
		switch(numVrts)
		{
			case 3:
				f = *grid.create<Triangle>(TriangleDescriptor(vLocalVrts[0], vLocalVrts[1],
																vLocalVrts[2]));
				break;
			case 4:
				f = *grid.create<Quadrilateral>(QuadrilateralDescriptor(vLocalVrts[0], vLocalVrts[1],
																	vLocalVrts[2], vLocalVrts[3]));
				break;
			default:
				LOG("  LoadGridFromART: bad number of vertices in face: " << numVrts << " (3 or 4 supported).\n");
				continue;
		};

		if(si < 1)
			si = -1;
		else if(si == 10000)
			si = 0;
		
		if(si != -1)
			sh.assign_subset(f, si);

	//	store the face in an array
		vFaces.push_back(f);
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

		int si = atoi(tok);
		
	//	read the indices
		ReadIndices(vInds, buf, delim, false);

	//	there have to be at least three edges
		if(vInds.size() < 3)
			continue;
			
	//	copy all vertices of the edges to vVrtDump
//TODO:	make sure that all indices are in the correct ranges
		vVrtDump.clear();
		for(size_t i = 0; i < vInds.size(); ++i){
			Face* f = vFaces[vInds[i]];
			for(size_t j = 0; j < f->num_vertices(); ++j)
				vVrtDump.push_back(f->vertex(j));
		}

	//	now collect the unique vertices in vVrtDump
		CollectUniqueObjects(vLocalVrts, vVrtDump);

//TODO: check orientation
	//	create the volume
		Volume* v = NULL;
		size_t numVrts = vLocalVrts.size();
		switch(numVrts)
		{
			case 4:
				v = *grid.create<Tetrahedron>(TetrahedronDescriptor(vLocalVrts[0], vLocalVrts[1],
																vLocalVrts[2], vLocalVrts[3]));
				break;
			case 5:
				v = *grid.create<Pyramid>(PyramidDescriptor(vLocalVrts[0], vLocalVrts[1],
															vLocalVrts[2], vLocalVrts[3],
															vLocalVrts[4]));
				break;
			case 6:
				v = *grid.create<Prism>(PrismDescriptor(vLocalVrts[0], vLocalVrts[1],
														vLocalVrts[2], vLocalVrts[3],
														vLocalVrts[4], vLocalVrts[5]));
				break;
			case 8:
				v = *grid.create<Hexahedron>(HexahedronDescriptor(vLocalVrts[0], vLocalVrts[1],
																vLocalVrts[2], vLocalVrts[3],
																vLocalVrts[4], vLocalVrts[5],
																vLocalVrts[6], vLocalVrts[7]));
				break;
			default:
				LOG("  LoadGridFromART: bad number of vertices in volume: " << numVrts << " (4, 5, 6 or 8 supported).\n");
				continue;
		};

		si -= 1;
		sh.assign_subset(v, si);
	}

//	delete the buffer
	delete[] buf;

	return true;
}

////////////////////////////////////////////////////////////////////////
bool SaveGridToART(Grid& srcGrid, const char* filename,
				   SubsetHandler* pSH, AVector3& aPos)
{
	if(!srcGrid.has_vertex_attachment(aPos)){
		LOG("  Aborting SaveGridToART: position attachment is missing.\n");
		return false;
	}
	
//	open the file
	ofstream out(filename);
	out.precision(20);
	out.setf(ios_base::scientific);
	
	if(!out)
		return false;
		
//	Options FACEOPT_AUTOGENERATE_EDGES and VOLOPT_AUTOGENERATE_FACES
//	have to be enabled. If they are not, we'll create a copy of the
//	grid and enable those options there.
//	note that we need a second subset handler too
	Grid tGrid;
	SubsetHandler tSH(tGrid);
	Grid* pGrid;	// Used to initialise the reference Grid& grid later on.
	if(!srcGrid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES
								| VOLOPT_AUTOGENERATE_FACES))
	{
	//	copy grid and subset-handler
		tGrid = srcGrid;
		if(pSH){
			tSH = *pSH;
			pSH = &tSH;
		}
		
		tGrid.enable_options(VOLOPT_AUTOGENERATE_FACES);
		tGrid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
		pGrid = &tGrid;
		
	//	make sure that the position attachment has been copied
		if(!tGrid.has_vertex_attachment(aPos)){
		//	copy it manually
			if(!CopyAttachments<VertexBase>(srcGrid, aPos, tGrid, aPos))
				return false;
		}
	}
	else {
		pGrid = &srcGrid;
	}
	
//	We'll work on this reference.
	Grid& grid = *pGrid;
	
//	write the header
	out << "%% Version 3.0" << endl;
	out << "%% VertexNumber: " << grid.num<VertexBase>() << endl;
	out << "%% EdgeNumber: " << grid.num<EdgeBase>() << endl;
	out << "%% FaceNumber: " << grid.num<Face>() << endl;
	out << "%% ElementNumber: " << grid.num<Volume>() << endl;
	out << "%% Translation: (0, 0, 0)" << endl;
	out << "%% Cornermark: 10001" << endl;
	out << "%% DO NOT CHANGE LINES ABOVE !!!!" << endl;
	out << "% NET: Vertices <-> Edges <-> Faces <-> Elements" << endl;
	
	
//	write vertices
	{
		out << "% Vertices: x y z" << endl;
		Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
		for(VertexBaseIterator iter = grid.begin<VertexBase>();
			iter != grid.end<VertexBase>(); ++iter)
		{
			out << aaPos[*iter].x << " ";
			out << aaPos[*iter].y << " ";
			out << aaPos[*iter].z << endl;
		}
	//	done - write end sign
		out << "$" << endl;
	}

//	write edges
	{
	//	vertices have to be associated with indices
		AInt aInt;
		grid.attach_to_vertices(aInt);
		Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
		AssignIndices<VertexBase>(grid.begin<VertexBase>(), grid.end<VertexBase>(), aaInt);

	//	the subset index
		int si = 0;
		
	//	begin writing
		out << "% Edges (Indices to List of Points):" << endl;
	//	iterate over edges
		for(EdgeBaseIterator iter = grid.begin<EdgeBase>();
			iter != grid.end<EdgeBase>(); ++iter)
		{
			EdgeBase* e = *iter;

			if(pSH){
				si = pSH->get_subset_index(e);
				if(si == 0)
					si = 10000;
				else if(si == -1)
					si = 0;
			}

			out << si << ": " << aaInt[e->vertex(0)] << " " << aaInt[e->vertex(1)] << endl;
		}

	//	done - write end sign
		out << "$" << endl;

	//	detach indices
		grid.detach_from_vertices(aInt);
	}

//	write faces
	{
		
	//	edges have to be associated with indices
		AInt aInt;
		grid.attach_to_edges(aInt);
		Grid::EdgeAttachmentAccessor<AInt> aaInt(grid, aInt);
		AssignIndices<EdgeBase>(grid.begin<EdgeBase>(), grid.end<EdgeBase>(), aaInt);
		
	//	the subset index
		int si = 0;
		
	//	we'll collect edges in here
		vector<EdgeBase*> vEdges;
		
	//	begin writing
		out << "% Faces (Indices to List of Edges):" << endl;
	//	iterate over faces
		for(FaceIterator iter = grid.begin<Face>();
			iter != grid.end<Face>(); ++iter)
		{
			Face* f = *iter;
			if(pSH){
				si = pSH->get_subset_index(f);
				if(si == 0)
					si = 10000;
				else if(si == -1)
					si = 0;
			}
			
		//	collect edges
			CollectEdges(vEdges, grid, f);
			
		//	write edges
			out << si << ":";
			for(size_t i = 0; i < vEdges.size(); ++i)
				out << " " << aaInt[vEdges[i]];
			out << endl;
		}
		
	//	done - write end sign
		out << "$" << endl;
		
	//	detach indices
		grid.detach_from_edges(aInt);
	}

//	write volumes
	{
		
	//	faces have to be associated with indices
		AInt aInt;
		grid.attach_to_faces(aInt);
		Grid::FaceAttachmentAccessor<AInt> aaInt(grid, aInt);
		AssignIndices<Face>(grid.begin<Face>(), grid.end<Face>(), aaInt);
		
	//	the subset index
		int si = 1;
		
	//	we'll collect faces in here
		vector<Face*> vFaces;
		
	//	begin writing
		out << "% Elements (Indices to List of Faces):" << endl;
	//	iterate over volumes
		for(VolumeIterator iter = grid.begin<Volume>();
			iter != grid.end<Volume>(); ++iter)
		{
			Volume* v = *iter;
			if(pSH){
				si = pSH->get_subset_index(v) + 1;
			}
			
		//	collect faces
			CollectFaces(vFaces, grid, v);
			
		//	write faces
			out << si << ":";
			for(size_t i = 0; i < vFaces.size(); ++i)
				out << " " << aaInt[vFaces[i]];
			out << endl;
		}
		
	//	done - write end sign
		out << "$" << endl;
		
	//	detach indices
		grid.detach_from_faces(aInt);
	}

	out.close();
	return true;
}
				   
}//	end of namespace
