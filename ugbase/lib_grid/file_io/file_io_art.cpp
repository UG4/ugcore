#include <fstream>
#include <vector>
#include <iostream>
#include <cstring>
#include "file_io_art.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/attachment_util.h"
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

//	uses grid::mark
static Tetrahedron*
CreateTetrahedron(Grid& grid, vector<Triangle*>& vTris,
				Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	const char* errorMsg = "ERROR in file_io_art.cpp CreateTetrahedron. ";
	
	if(vTris.size() != 4){
		UG_LOG(errorMsg << "Bad number of triangles: " << vTris.size() << endl);
		return NULL;
	}
	
	Vertex* vrts[4];
	int vrtCount = 0;
	
//	get the 4 points of the tetrahedron
	
	grid.begin_marking();

	for(size_t i = 0; i < 4; ++i){
		for(size_t j = 0; j < 3; ++j){
			Vertex* vrt = vTris[i]->vertex(j);
			if(!grid.is_marked(vrt)){
			//	make sure that we won't collect too many vertices.
				if(vrtCount == 4){
					UG_LOG(errorMsg << "Triangles do not build a tetrahedron.\n");
					return NULL;
				}
				
				grid.mark(vrt);
				vrts[vrtCount++] = vrt;
			}
		}
	}
	
	grid.end_marking();

	if(vrtCount < 4){
		UG_LOG(errorMsg << "Triangles have less than 4 distinct vertices.\n");
		return NULL;
	}
	
//	we have to check the orientation of the tetrahedron
//	compare the normal of the first triangle to the top-vertex
	if(PointFaceTest(aaPos[vrts[3]], vTris[0], aaPos) < 0){
	//	we have to change the order of the first three vrts
		swap(vrts[0], vrts[1]);
	}
	
//	create the tetrahedron
	return *grid.create<Tetrahedron>(
				TetrahedronDescriptor(vrts[0], vrts[1],
									  vrts[2], vrts[3]));
	
}

//	uses grid::mark
static Pyramid*
CreatePyramid(Grid& grid, vector<Triangle*>& vTris,
				vector<Quadrilateral*>& vQuads,
				Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	const char* errorMsg = "ERROR in file_io_art.cpp CreatePyramid. ";
	
	if(vTris.size() != 4){
		UG_LOG(errorMsg << "Bad number of triangles: " << vTris.size() << endl);
		return NULL;
	}

	if(vQuads.size() != 1){
		UG_LOG(errorMsg << "Bad number of quadrilaterals: " << vQuads.size() << endl);
		return NULL;
	}
		
	Vertex* vrts[5];
	int vrtCount = 0;
	
	grid.begin_marking();
	
//	get the 5 points of the pyramid
//	take the vertices of the base-quad first
	for(size_t j = 0; j < 4; ++j){
		Vertex* vrt = vQuads[0]->vertex(j);
		if(!grid.is_marked(vrt)){
			grid.mark(vrt);
			vrts[vrtCount++] = vrt;
		}
	}
	
//	now find the top vertex by checking the first triangle
	for(size_t j = 0; j < 3; ++j){
		Vertex* vrt = vTris[0]->vertex(j);
		if(!grid.is_marked(vrt)){
			grid.mark(vrt);
			vrts[vrtCount++] = vrt;
			break;
		}
	}
	
	grid.end_marking();
	
//	make sure that all vertices have been found
	if(vrtCount != 5){
		UG_LOG(errorMsg << "Faces do not build a pyramid.\n");
		return NULL;
	}
	
//	we have to check the orientation of the prism
//	compare the normal of the base-quad to the top
	if(PointFaceTest(aaPos[vrts[4]], vQuads[0], aaPos) < 0){
	//	we have to change the order of the first four vrts
		swap(vrts[1], vrts[2]);
		swap(vrts[0], vrts[3]);
	}

//	create the pyramid
	return *grid.create<Pyramid>(
				PyramidDescriptor(vrts[0], vrts[1],
							vrts[2], vrts[3], vrts[4]));
}


//	this method assumes that the given faces already exist in the grid
//	together with their associated edges
static Prism*
CreatePrism(Grid& grid, vector<Triangle*>& vTris,
				vector<Quadrilateral*>& vQuads,
				Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	const char* errorMsg = "ERROR in file_io_art.cpp CreatePrism. ";
	
	if(vTris.size() != 2){
		UG_LOG(errorMsg << "Bad number of triangles: " << vTris.size() << endl);
		return NULL;
	}

	if(vQuads.size() != 3){
		UG_LOG(errorMsg << "Bad number of quadrilaterals: " << vQuads.size() << endl);
		return NULL;
	}
		
	Vertex* vrts[6];
	int vrtCount = 0;
	
//	get the 6 points of the prism
//	calculate the center on the fly
	vector3 vCenter(0, 0, 0);
	
	for(size_t i = 0; i < 2; ++i){
		for(size_t j = 0; j < 3; ++j){
			vrts[vrtCount++] = vTris[i]->vertex(j);
			VecAdd(vCenter, vCenter, aaPos[vTris[i]->vertex(j)]);
		}
	}
	
	VecScale(vCenter, vCenter, 1. / 6.);
	
//	we have to check the orientation of the prism
//	compare the normal of the first triangle to the center
	if(PointFaceTest(vCenter, vTris[0], aaPos) < 0){
	//	we have to change the order of the first three vrts
		swap(vrts[0], vrts[1]);
	}

//	compare the normal of the second triangle to the center
	if(PointFaceTest(vCenter, vTris[1], aaPos) > 0){
	//	we have to change the order of the last three vrts
		swap(vrts[3], vrts[4]);
	}
	
//	we have to find the vertex of the second triangle that
//	is connected to the first vertex of the first triangle
	int index = 0;
	for(; index < 3; ++index){
		if(grid.get_edge(vrts[0], vrts[3 + index]))
			break;
	}
	
	if(index == 3)
		return NULL;

//	create the tetrahedron
	return *grid.create<Prism>(
				PrismDescriptor(vrts[0], vrts[1], vrts[2],
								vrts[3 + index],
								vrts[3 + (index + 1) % 3],
								vrts[3 + (index + 2) % 3]));	
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
	vector<Vertex*> vVrtDump;
	vector<Vertex*> vLocalVrts;
	
//	read the vertices
	vector<RegularVertex*>	vVrts;

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
		RegularVertex* v = *grid.create<RegularVertex>();
	//	read the coordinates
//TODO:	make sure that everything is ok.
		aaPos[v].x() = atof(tok);
		tok = strtok(NULL, delim);
		aaPos[v].y() = atof(tok);
		tok = strtok(NULL, delim);
		aaPos[v].z() = atof(tok);

	//	store the vertex in an array
		vVrts.push_back(v);
	}

//	read the edges
	vector<RegularEdge*>	vEdges;

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
		int i1, i2;
		//int si = atoi(tok);
		tok = strtok(NULL, delim);
		i1 = atoi(tok);
		tok = strtok(NULL, delim);
		i2 = atoi(tok);

	//	create a new edge
		RegularEdge* e = *grid.create<RegularEdge>(EdgeDescriptor(vVrts[i1], vVrts[i2]));
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
			Edge* e = vEdges[vInds[i]];
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
	vector<Triangle*> vTris;
	vector<Quadrilateral*> vQuads;
	
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

	//	there have to be at least three faces
		if(vInds.size() < 3)
			continue;
			
	//	collect faces in vTris and vQuads
		vTris.clear();
		vQuads.clear();
		for(size_t i = 0; i < vInds.size(); ++i){
			Face* f = vFaces[vInds[i]];
			if(f->num_vertices() == 3){	
				Triangle* t = dynamic_cast<Triangle*>(f);
				if(t)
					vTris.push_back(t);
				else{
					UG_LOG("  bad triangle in LoadGridFromART!");
				}
			}
			else if(f->num_vertices() == 4){
				Quadrilateral* q = dynamic_cast<Quadrilateral*>(f);
				if(q)
					vQuads.push_back(q);
				else{
					UG_LOG("  bad quadrilateral in LoadGridFromART!");
				}
			}
			else {
				UG_LOG("  Bad face in LoadGridFromART! Aborting...");
				return false;
			}
		}

	//	get the type of volume
		int volType = ROID_UNKNOWN;
		if(vTris.size() == 4 && vQuads.size() == 0)
			volType = ROID_TETRAHEDRON;
		else if(vTris.size() == 4 && vQuads.size() == 1)
			volType = ROID_PYRAMID;
		else if(vTris.size() == 2 && vQuads.size() == 3)
			volType = ROID_PRISM;
		else if(vTris.size() == 0 && vQuads.size() == 6)
			volType = ROID_HEXAHEDRON;


	//	create the volume
		Volume* v = NULL;

		switch(volType)
		{
			case ROID_TETRAHEDRON:
				v = CreateTetrahedron(grid, vTris, aaPos);
				break;
			case ROID_PYRAMID:
				v = CreatePyramid(grid, vTris, vQuads, aaPos);
				break;
			case ROID_PRISM:
				v = CreatePrism(grid, vTris, vQuads, aaPos);
				break;
			case ROID_HEXAHEDRON:
//TODO: create hexahedron
				break;
			default:
				LOG("  LoadGridFromART: bad volume type. volume has "
					<< vTris.size() << " triangles and "
					<< vQuads.size() << " quadrilaterals.\n");
				continue;
		};

		si -= 1;
		if(v)
			sh.assign_subset(v, si);
		else {
			LOG("  LoadGridFromART: could not create volume element.\n");
		}
	}

//	delete the buffer
	delete[] buf;

	return true;
}

////////////////////////////////////////////////////////////////////////
bool SaveGridToART(Grid& srcGrid, const char* filename,
				   ISubsetHandler* pSH, AVector3& aPos)
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
			if(!CopyAttachments<Vertex>(srcGrid, aPos, tGrid, aPos))
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
	out << "%% VertexNumber: " << grid.num<Vertex>() << endl;
	out << "%% EdgeNumber: " << grid.num<Edge>() << endl;
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
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter)
		{
			out << aaPos[*iter].x() << " ";
			out << aaPos[*iter].y() << " ";
			out << aaPos[*iter].z() << endl;
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
		AssignIndices<Vertex>(grid.begin<Vertex>(), grid.end<Vertex>(), aaInt);

	//	the subset index
		int si = 0;
		
	//	begin writing
		out << "% Edges (Indices to List of Points):" << endl;
	//	iterate over edges
		for(EdgeIterator iter = grid.begin<Edge>();
			iter != grid.end<Edge>(); ++iter)
		{
			Edge* e = *iter;

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
		AssignIndices<Edge>(grid.begin<Edge>(), grid.end<Edge>(), aaInt);
		
	//	the subset index
		int si = 0;
		
	//	we'll collect edges in here
		vector<Edge*> vEdges;
		
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
