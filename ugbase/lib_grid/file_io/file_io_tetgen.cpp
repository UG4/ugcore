//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d18

#include <fstream>
#include "file_io_tetgen.h"
#include "common/util/string_util.h"
#include "../lg_base.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
bool LoadGridFromELE(Grid& grid, const char* filename, ISubsetHandler* pSH,
					APosition& aPos)
{
//	build the correct filenames
	string sElements = filename;
	string sNodes, sFaces;
	if(sElements.find(".ele", sElements.size() - 4) != string::npos)
	{
	//	the filename ends as we would expect
		sNodes = sFaces = sElements.substr(0, sElements.size() - 4);
		sNodes.append(".node");
		sFaces.append(".face");
	}
	else
	{
		LOG("Problem in LoadGridFromELE with file " << filename << ". filename has to end with .ele\n");
		return false;
	}
		
	return ImportGridFromTETGEN(grid, sNodes.c_str(), sFaces.c_str(), sElements.c_str(),
								aPos, pSH);
}

////////////////////////////////////////////////////////////////////////					
bool SaveGridToELE(Grid& grid, const char* filename, ISubsetHandler* pSH,
					APosition& aPos)
{
//	build the correct filenames
	string sElements = filename;
	string sNodes, sFaces;
	if(sElements.find(".ele", sElements.size() - 4) != string::npos)
	{
	//	the filename ends as we would expect
		sNodes = sFaces = sElements.substr(0, sElements.size() - 4);
		sNodes.append(".node");
		sFaces.append(".face");
	}
	else
	{
		sNodes = sFaces = sElements;
		sElements.append(".ele");
		sNodes.append(".node");
		sFaces.append(".face");
	}
	
//	give a warning if the grid doesn't consist of tetrahedrons only.
	if(grid.num<Tetrahedron>() < grid.num<Face>()){
		UG_LOG("  INFO in SaveGridToELE: Non-tetrahedral elements will be skipped.\n");
	}
	
//	if a subset-handler was supplied we have to copy its subset-indices
//	to a simple int-attachment.
	if(pSH)
	{
		AInt aInt;
		grid.attach_to_volumes(aInt);
		Grid::VolumeAttachmentAccessor<AInt> aaInt(grid, aInt);
		
		for(VolumeIterator iter = grid.begin<Volume>(); iter != grid.end<Volume>(); ++iter)
			aaInt[*iter] = pSH->get_subset_index(*iter);	
		
	//	face subset-indices
		grid.attach_to_faces(aInt);
		Grid::FaceAttachmentAccessor<AInt> aaIntFace(grid, aInt);
		for(FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>(); ++iter)
			aaIntFace[*iter] = pSH->get_subset_index(*iter);	

		bool retVal = ExportGridToTETGEN(grid, sNodes.c_str(), sFaces.c_str(),
								sElements.c_str(), aPos, NULL, NULL, &aInt, &aInt);
								
		grid.detach_from_faces(aInt);
		grid.detach_from_volumes(aInt);
		
		return retVal;
	}
	
	return ExportGridToTETGEN(grid, sNodes.c_str(), sFaces.c_str(), sElements.c_str(), aPos);
	
}
					
////////////////////////////////////////////////////////////////////////
//	ImportGridFromTETGEN
bool ImportGridFromTETGEN(Grid& grid,
						const char* nodesFilename, const char* facesFilename,
						const char* elemsFilename, AVector3& aPos,
						std::vector<AFloat>* pvNodeAttributes,
						AInt* paNodeBoundaryMarker,
						AInt* paFaceBoundaryMarker,
						AInt* paElementAttribute)
{
//	read nodes and store them in an array for index access

	vector<RegularVertex*>	vVertices;

	{
		ifstream in(nodesFilename);
		if(!in)
		{
			LOG("WARNING in ImportGridFromTETGEN: nodes file not found: " << nodesFilename << endl);
			return false;
		}

		uint numNodes, dim, numAttribs, numBoundaryMarkers;
		in >> numNodes;
		in >> dim;
		in >> numAttribs;
		in >> numBoundaryMarkers;

	//	set up attachment accessors
		if(!grid.has_vertex_attachment(aPos))
			grid.attach_to_vertices(aPos);
		Grid::VertexAttachmentAccessor<AVector3> aaPosVRT(grid, aPos);

		Grid::VertexAttachmentAccessor<AInt> aaBMVRT;
		if(paNodeBoundaryMarker != NULL)
			aaBMVRT.access(grid, *paNodeBoundaryMarker);

		vector<Grid::VertexAttachmentAccessor<AFloat> > vaaAttributesVRT;
		if(pvNodeAttributes != NULL)
		{
			vaaAttributesVRT.resize(pvNodeAttributes->size());
			for(uint i = 0; i < pvNodeAttributes->size(); ++i)
				vaaAttributesVRT[i].access(grid, (*pvNodeAttributes)[i]);
		}

	//	read the vertices
		int index;
		for(uint i = 0; i < numNodes; ++i)
		{
			RegularVertex* v = *grid.create<RegularVertex>();
			vVertices.push_back(v);

		//	read index and coords
			in >> index;
			in >> aaPosVRT[v].x();
			in >> aaPosVRT[v].y();
			in >> aaPosVRT[v].z();

		//	read attributes
			if(numAttribs > 0)
			{
				for(uint j = 0; j < numAttribs; ++j)
				{
					float tmp;
					in >> tmp;
					if(j < vaaAttributesVRT.size())
						(vaaAttributesVRT[j])[v] = tmp;
				}
			}

		//	read boundary marker
			if(numBoundaryMarkers > 0)
			{
				int bm;
				in >> bm;
				if(paNodeBoundaryMarker != NULL)
					aaBMVRT[v] = bm;
			}
		}

		in.close();
	}

//	read faces
	if(facesFilename != NULL)
	{
		ifstream in(facesFilename);
		if(in)
		{
			int numFaces, numBoundaryMarkers;
			in >> numFaces;
			in >> numBoundaryMarkers;

			Grid::FaceAttachmentAccessor<AInt> aaBMFACE;
			if(paFaceBoundaryMarker != NULL)
				aaBMFACE.access(grid, *paFaceBoundaryMarker);

			for(int i = 0; i < numFaces; ++i)
			{
				int index, i1, i2, i3;
				in >> index;
				in >> i1;
				in >> i2;
				in >> i3;

				Triangle* t = *grid.create<Triangle>(TriangleDescriptor(vVertices[i1], vVertices[i2], vVertices[i3]));



				if(numBoundaryMarkers > 0)
				{
					int bm;
					in >> bm;
					if(aaBMFACE.valid())
						aaBMFACE[t] = bm;
				}
				else
				{
					if(aaBMFACE.valid())
						aaBMFACE[t] = -1;
				}
			}
		}
		else
			LOG("WARNING in ImportGridFromTETGEN: faces file not found: " << facesFilename << endl);

		in.close();
	}

//	read volumes
	if(elemsFilename != NULL)
	{
		ifstream in(elemsFilename);
		if(in)
		{
			int numTets, numNodesPerTet, numAttribs;
			in >> numTets;
			in >> numNodesPerTet;
			in >> numAttribs;

		//	attachment accessors:
			Grid::VolumeAttachmentAccessor<AInt> aaAttributeVOL;
			if(paElementAttribute)
				aaAttributeVOL.access(grid, *paElementAttribute);

			vector<int> vTetNodes(numNodesPerTet);
			for(int i = 0; i < numTets; ++i)
			{
				int index;
				in >> index;
				for(int j = 0; j < numNodesPerTet; ++j)
					in >> vTetNodes[j];

				Tetrahedron* t = *grid.create<Tetrahedron>(TetrahedronDescriptor(
											vVertices[vTetNodes[0]], vVertices[vTetNodes[1]],
											vVertices[vTetNodes[2]], vVertices[vTetNodes[3]]));

				if(numAttribs > 0)
				{
					int a;
					in >> a;
					if(paElementAttribute != NULL)
						aaAttributeVOL[t] = a;
				}
			}
		}
		else
			LOG("WARNING in ImportGridFromTETGEN: elems file not found: " << elemsFilename << endl);

		in.close();
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
//	ImportGridFromTETGEN
bool ImportGridFromTETGEN(Grid& grid,
						const char* nodesFilename, const char* facesFilename,
						const char* elemsFilename, AVector3& aPos,
						ISubsetHandler* psh,
						std::vector<AFloat>* pvNodeAttributes)
{
//	read nodes and store them in an array for index access

	vector<RegularVertex*>	vVertices;

	{
		ifstream in(nodesFilename);
		if(!in)
		{
			LOG("WARNING in ImportGridFromTETGEN: nodes file not found: " << nodesFilename << endl);
			return false;
		}

		uint numNodes, dim, numAttribs, numBoundaryMarkers;
		in >> numNodes;
		in >> dim;
		in >> numAttribs;
		in >> numBoundaryMarkers;

	//	set up attachment accessors
		if(!grid.has_vertex_attachment(aPos))
			grid.attach_to_vertices(aPos);
		Grid::VertexAttachmentAccessor<AVector3> aaPosVRT(grid, aPos);

		vector<Grid::VertexAttachmentAccessor<AFloat> > vaaAttributesVRT;
		if(pvNodeAttributes != NULL)
		{
			vaaAttributesVRT.resize(pvNodeAttributes->size());
			for(uint i = 0; i < pvNodeAttributes->size(); ++i)
				vaaAttributesVRT[i].access(grid, (*pvNodeAttributes)[i]);
		}

	//	read the vertices
		int index;
		for(uint i = 0; i < numNodes; ++i)
		{
			RegularVertex* v = *grid.create<RegularVertex>();
			vVertices.push_back(v);

		//	read index and coords
			in >> index;
			in >> aaPosVRT[v].x();
			in >> aaPosVRT[v].y();
			in >> aaPosVRT[v].z();

		//	read attributes
			if(numAttribs > 0)
			{
				for(uint j = 0; j < numAttribs; ++j)
				{
					float tmp;
					in >> tmp;
					if(j < vaaAttributesVRT.size())
						(vaaAttributesVRT[j])[v] = tmp;
				}
			}

		//	read boundary marker
			if(numBoundaryMarkers > 0)
			{
				int bm;
				in >> bm;
				if(psh != NULL)
					psh->assign_subset(v, bm);
			}
		}

		in.close();
	}

//	read faces
	if(facesFilename != NULL)
	{
		ifstream in(facesFilename);
		if(in)
		{
			int numFaces, numBoundaryMarkers;
			in >> numFaces;
			in >> numBoundaryMarkers;

			for(int i = 0; i < numFaces; ++i)
			{
				int index, i1, i2, i3;
				in >> index;
				in >> i1;
				in >> i2;
				in >> i3;

				Triangle* t = *grid.create<Triangle>(TriangleDescriptor(vVertices[i1], vVertices[i2], vVertices[i3]));

				if(numBoundaryMarkers > 0){
					int bm;
					in >> bm;
					if(psh != NULL){
						psh->assign_subset(t, bm);
					}
				}
			}
		}
		else
			LOG("WARNING in ImportGridFromTETGEN: faces file not found: " << facesFilename << endl);

		in.close();
	}

//	read volumes
	if(elemsFilename != NULL)
	{
		ifstream in(elemsFilename);
		if(in)
		{
			int numTets, numNodesPerTet, numAttribs;
			in >> numTets;
			in >> numNodesPerTet;
			in >> numAttribs;

			vector<int> vTetNodes(numNodesPerTet);
			for(int i = 0; i < numTets; ++i)
			{
				int index;
				in >> index;
				for(int j = 0; j < numNodesPerTet; ++j)
					in >> vTetNodes[j];

				Tetrahedron* t = *grid.create<Tetrahedron>(TetrahedronDescriptor(
											vVertices[vTetNodes[0]], vVertices[vTetNodes[1]],
											vVertices[vTetNodes[2]], vVertices[vTetNodes[3]]));

				int a = 0;
				if(numAttribs > 0)
					in >> a;
				if(psh != NULL)
					psh->assign_subset(t, a);
			}
		}
		else
			LOG("WARNING in ImportGridFromTETGEN: elems file not found: " << elemsFilename << endl);

		in.close();
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
//	ExportGridToSMESH
bool ExportGridToSMESH(Grid& grid, const char* filename, AVector3& aPos,
						std::vector<AFloat>* pvNodeAttributes,
						AInt* paNodeBoundaryMarker,
						AInt* paFaceBoundaryMarker,
						std::vector<vector3>* pvHoles,
						std::vector<vector3>* pvRegionPositions,
						std::vector<int>* pvRegionAttributes,
						std::vector<float>* pvRegionVolumeConstraints)

{
	if(!grid.has_vertex_attachment(aPos))
			return false;

	ofstream out(filename);

	if(!out)
		return false;

//	check if regions are specified in the correct way.
	if(pvRegionPositions != NULL)
	{
		if(pvRegionAttributes != NULL)
			if(pvRegionAttributes->size() != pvRegionPositions->size())
				return false;

		if(pvRegionVolumeConstraints != NULL)
			if(pvRegionVolumeConstraints->size() != pvRegionPositions->size())
				return false;
	}

	vector<int> vTmpRegionAttributes;
	if((pvRegionAttributes == NULL) && (pvRegionPositions != NULL) && (pvRegionVolumeConstraints != NULL))
	{
	//	attributes have to be supplied if constraints are given.
	//	since no attributes were passed, generate your own.
		pvRegionAttributes = &vTmpRegionAttributes;
		for(uint i = 0; i < pvRegionPositions->size(); ++i)
			pvRegionAttributes->push_back(0);
	}

//	write points. store indices with each point for later array-like access.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaIntVRT(grid, aInt);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);

	{
		int numAttribs = 0;
		int numBoundaryMarkers = 0;

	//	attachment-accessors for the nodes attributes
		vector<Grid::VertexAttachmentAccessor<AFloat> > vaaFloatVRT;
		if(pvNodeAttributes != NULL)
		{
			numAttribs = pvNodeAttributes->size();
			vaaFloatVRT.resize(pvNodeAttributes->size());
			for(uint i = 0; i < pvNodeAttributes->size(); ++i)
				vaaFloatVRT[i].access(grid, (*pvNodeAttributes)[i]);
		}

	//	attachment-accessor for the nodes boundary-marker
		Grid::VertexAttachmentAccessor<AInt> aaBMVRT;
		if(paNodeBoundaryMarker != NULL)
		{
			numBoundaryMarkers = 1;
			aaBMVRT.access(grid, *paNodeBoundaryMarker);
		}

	//	write number of nodes, dimension, number of attributes, boundary markers (0 or 1)
		out << grid.num_vertices() << " 3 " << numAttribs << " " << numBoundaryMarkers << endl;

		int counter = 0;
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++, counter++)
		{
			aaIntVRT[*iter] = counter;
			out << counter << " " <<	aaPos[*iter].x() << " " <<
										aaPos[*iter].y() << " " <<
										aaPos[*iter].z();

		//	write attributes:
			for(uint i = 0; i < vaaFloatVRT.size(); ++i)
			{
				out << " " << (vaaFloatVRT[i])[*iter];
			}

		//	write boundary markers:
			if(paNodeBoundaryMarker != NULL)
				out << " " << aaBMVRT[*iter];

			out << endl;

		}
	}

	out << endl;

//	write facets
	{
		Grid::FaceAttachmentAccessor<AInt> aaBMFACE;
		if(paFaceBoundaryMarker != NULL)
		{
			out << grid.num_faces() << " 1" << endl;
			aaBMFACE.access(grid, *paFaceBoundaryMarker);
		}
		else
			out << grid.num_faces() << " 0" << endl;

		for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter)
		{
			Face* f = *iter;
			out << f->num_vertices();

			for(uint i = 0; i < f->num_vertices(); ++i)
				out << " " << aaIntVRT[f->vertex(i)];

			if(paFaceBoundaryMarker != NULL)
				out << " " << aaBMFACE[f];

			out << endl;
		}
	}

	out << endl;

//	write holes
	if(pvHoles != NULL)
	{
		vector<vector3>& vHoles = *pvHoles;
		out << vHoles.size() << endl;
		for(uint i = 0; i < vHoles.size(); ++i)
		{
			out << i << " " <<	vHoles[i].x() << " " <<
								vHoles[i].y() << " " <<
								vHoles[i].z() << endl;
		}
	}
	else
	{
		out << "0" << endl;
	}

	out << endl;

//	write regions (optional)
	if(pvRegionPositions != NULL)
	{
		vector<vector3>& vRegionPositions = *pvRegionPositions;
		out << vRegionPositions.size() << endl;
		for(uint i = 0; i < vRegionPositions.size(); ++i)
		{
			out << i << " " << 	vRegionPositions[i].x() << " " <<
								vRegionPositions[i].y() << " " <<
								vRegionPositions[i].z();

			if(pvRegionAttributes != NULL)
				out << " " << (*pvRegionAttributes)[i];

			if(pvRegionVolumeConstraints != NULL)
				out << " " << (*pvRegionVolumeConstraints)[i];

			out << endl;
		}
	}

	out.close();

	grid.detach_from_vertices(aInt);

	return true;
}


////////////////////////////////////////////////////////////////////////
//	ExportGridToSMESH
bool ExportGridToSMESH(Grid& grid, const char* filename, AVector3& aPos,
						ISubsetHandler* psh,
						std::vector<AFloat>* pvNodeAttributes,
						std::vector<vector3>* pvHoles,
						std::vector<vector3>* pvRegionPositions,
						std::vector<int>* pvRegionAttributes,
						std::vector<float>* pvRegionVolumeConstraints)

{
	if(!grid.has_vertex_attachment(aPos))
			return false;

	ofstream out(filename);

	if(!out)
		return false;

//	check if regions are specified in the correct way.
	if(pvRegionPositions != NULL)
	{
		if(pvRegionAttributes != NULL)
			if(pvRegionAttributes->size() != pvRegionPositions->size())
				return false;

		if(pvRegionVolumeConstraints != NULL)
			if(pvRegionVolumeConstraints->size() != pvRegionPositions->size())
				return false;
	}

	vector<int> vTmpRegionAttributes;
	if((pvRegionAttributes == NULL) && (pvRegionPositions != NULL) && (pvRegionVolumeConstraints != NULL))
	{
	//	attributes have to be supplied if constraints are given.
	//	since no attributes were passed, generate your own.
		pvRegionAttributes = &vTmpRegionAttributes;
		for(uint i = 0; i < pvRegionPositions->size(); ++i)
			pvRegionAttributes->push_back(0);
	}

//	write points. store indices with each point for later array-like access.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaIntVRT(grid, aInt);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);

	{
		int numAttribs = 0;
		int numBoundaryMarkers = 0;

	//	attachment-accessors for the nodes attributes
		vector<Grid::VertexAttachmentAccessor<AFloat> > vaaFloatVRT;
		if(pvNodeAttributes != NULL)
		{
			numAttribs = pvNodeAttributes->size();
			vaaFloatVRT.resize(pvNodeAttributes->size());
			for(uint i = 0; i < pvNodeAttributes->size(); ++i)
				vaaFloatVRT[i].access(grid, (*pvNodeAttributes)[i]);
		}

	//	attachment-accessor for the nodes boundary-marker
		if(psh != NULL)
			numBoundaryMarkers = 1;

	//	write number of nodes, dimension, number of attributes, boundary markers (0 or 1)
		out << grid.num_vertices() << " 3 " << numAttribs << " " << numBoundaryMarkers << endl;

		int counter = 0;
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++, counter++)
		{
			aaIntVRT[*iter] = counter;
			out << counter << " " <<	aaPos[*iter].x() << " " <<
										aaPos[*iter].y() << " " <<
										aaPos[*iter].z();

		//	write attributes:
			for(uint i = 0; i < vaaFloatVRT.size(); ++i)
				out << " " << (vaaFloatVRT[i])[*iter];

		//	write boundary markers:
			if(psh != NULL)
				out << " " << psh->get_subset_index(*iter);

			out << endl;

		}
	}

	out << endl;

//	write facets
	{
		if(psh != NULL)
			out << grid.num_faces() << " 1" << endl;
		else
			out << grid.num_faces() << " 0" << endl;

		for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter)
		{
			Face* f = *iter;
			out << f->num_vertices();

			for(uint i = 0; i < f->num_vertices(); ++i)
				out << " " << aaIntVRT[f->vertex(i)];

			if(psh != NULL)
				out << " " << psh->get_subset_index(f);

			out << endl;
		}
	}

	out << endl;

//	write holes
	if(pvHoles != NULL)
	{
		vector<vector3>& vHoles = *pvHoles;
		out << vHoles.size() << endl;
		for(uint i = 0; i < vHoles.size(); ++i)
		{
			out << i << " " <<	vHoles[i].x() << " " <<
								vHoles[i].y() << " " <<
								vHoles[i].z() << endl;
		}
	}
	else
	{
		out << "0" << endl;
	}

	out << endl;

//	write regions (optional)
	if(pvRegionPositions != NULL)
	{
		vector<vector3>& vRegionPositions = *pvRegionPositions;
		out << vRegionPositions.size() << endl;
		for(uint i = 0; i < vRegionPositions.size(); ++i)
		{
			out << i << " " << 	vRegionPositions[i].x() << " " <<
								vRegionPositions[i].y() << " " <<
								vRegionPositions[i].z();

			if(pvRegionAttributes != NULL)
				out << " " << (*pvRegionAttributes)[i];

			if(pvRegionVolumeConstraints != NULL)
				out << " " << (*pvRegionVolumeConstraints)[i];

			out << endl;
		}
	}

	out.close();

	grid.detach_from_vertices(aInt);

	return true;
}

////////////////////////////////////////////////////////////////////////
bool LoadGridFromSMESH(Grid& grid, const char* filename, AVector3& aPos,
						ISubsetHandler* psh)
{
	ifstream file(filename);
	UG_COND_THROW(!file, "Couldn't open " << filename << " for reading.");

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);

	enum ReadState{
		READ_VRT_HEADER,
		READ_VERTICES,
		READ_FACE_HEADER,
		READ_FACES
	};

	ReadState readState = READ_VRT_HEADER;

	size_t numVrts = 0;
	size_t dim = 0;
	size_t numAttribs = 0;
	size_t bndMarker = 0;
	size_t numFaces = 0;
	size_t numFacesRead = 0;
	vector<Vertex*> vrts;

	size_t curLine = 0;

	while(!file.eof()){
		++curLine;

		string line;
		getline(file, line);
		string trimmedLine = TrimString(line);
		
		if(trimmedLine.empty())
			continue;

		if(trimmedLine[0] == '#')
			continue;
		
		stringstream in(trimmedLine);

		try{
			switch(readState){
				case READ_VRT_HEADER:{
					in >> numVrts >> dim >> numAttribs >> bndMarker;
					UG_COND_THROW(!in, "Couldn't read vertex header");
					if(numVrts == 0)
						return true;
					readState = READ_VERTICES;
				}break;


				case READ_VERTICES:{
					size_t vrtInd;
					int tmp;
					Vertex* vrt = *grid.create<RegularVertex>();
					in >> vrtInd >> aaPos[vrt].x() >> aaPos[vrt].y() >> aaPos[vrt].z();
					for(size_t i = 0; i < numAttribs; ++i)
						in >> tmp;
					
					if(psh && bndMarker){
						in >> tmp;
						psh->assign_subset(vrt, tmp);
					}


					UG_COND_THROW(!in, "Couldn't read vertex");

					vrts.push_back(vrt);
					if(vrts.size() >= numVrts)
						readState = READ_FACE_HEADER;
				}break;


				case READ_FACE_HEADER:{
					in >> numFaces >> bndMarker;
					UG_COND_THROW(!in, "Couldn't read face header");
					if(numFaces == 0)
						return true;
					readState = READ_FACES;
				}break;


				case READ_FACES:{
					size_t numCorners;
					size_t inds[4];
					
					in >> numCorners;
					for(size_t i = 0; i < numCorners; ++i)
						in >> inds[i];

					Face* f = NULL;
					switch(numCorners){
						case 3:{
							f = *grid.create<Triangle>(
										TriangleDescriptor(
											vrts[inds[0]], vrts[inds[1]],
											vrts[inds[2]]));
						}break;

						case 4:{
							f = *grid.create<Quadrilateral>(
										QuadrilateralDescriptor(
											vrts[inds[0]], vrts[inds[1]],
											vrts[inds[2]], vrts[inds[3]]));
						}break;

						default:
							UG_THROW("Faces with " << numCorners << " corners are currently not supported");
					}

					if(psh && bndMarker){
						int si;
						in >> si;
						psh->assign_subset(f, si);
					}

					++numFacesRead;
					if(numFacesRead >= numFaces)
						return true;
				}break;
			}
		}
		UG_CATCH_THROW("LoadGridFromSMESH: Failed to interprete data in line " << curLine
						<< " of file '" << filename << "'");
	}

	return false;
}


////////////////////////////////////////////////////////////////////////
//	ExportGridToTETGEN
bool ExportGridToTETGEN(Grid& grid, const char* nodesFilename,
						const char* facesFilename, const char* elemsFilename,
						AVector3& aPos, std::vector<AFloat>* pvNodeAttributes,
						AInt* paNodeBoundaryMarker,
						AInt* paFaceBoundaryMarker,
						AInt* paElementAttribute)
{
	if(!grid.has_vertex_attachment(aPos))
			return false;

//	write nodes
	Grid::VertexAttachmentAccessor<AInt> aaIntVRT;
	{
		ofstream out(nodesFilename);

		if(!out)
			return false;

	//	write points. store indices with each point for later array-like access.
		AInt aInt;
		grid.attach_to_vertices(aInt);
		Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
		aaIntVRT.access(grid, aInt);


		int numAttribs = 0;
		int numBoundaryMarkers = 0;

	//	attachment-accessors for the nodes attributes
		vector<Grid::VertexAttachmentAccessor<AFloat> > vaaFloatVRT;
		if(pvNodeAttributes != NULL)
		{
			numAttribs = pvNodeAttributes->size();
			vaaFloatVRT.resize(pvNodeAttributes->size());
			for(uint i = 0; i < pvNodeAttributes->size(); ++i)
				vaaFloatVRT[i].access(grid, (*pvNodeAttributes)[i]);
		}

	//	attachment-accessor for the nodes boundary-marker
		Grid::VertexAttachmentAccessor<AInt> aaBMVRT;
		if(paNodeBoundaryMarker != NULL)
		{
			numBoundaryMarkers = 1;
			aaBMVRT.access(grid, *paNodeBoundaryMarker);
		}

	//	write number of nodes, dimension, number of attributes, boundary markers (0 or 1)
		out << grid.num_vertices() << " 3 " << numAttribs << " " << numBoundaryMarkers << endl;

		int counter = 0;
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++, counter++)
		{
			aaIntVRT[*iter] = counter;
			out << counter << " " <<	aaPos[*iter].x() << " " <<
										aaPos[*iter].y() << " " <<
										aaPos[*iter].z();

		//	write attributes:
			for(uint i = 0; i < vaaFloatVRT.size(); ++i)
			{
				out << " " << (vaaFloatVRT[i])[*iter];
			}

		//	write boundary markers:
			if(paNodeBoundaryMarker != NULL)
				out << " " << aaBMVRT[*iter];

			out << endl;
		}

		out.close();
	}

//	write facets
	if(facesFilename != NULL)
	{
		ofstream out(facesFilename);
		if(out)
		{
			Grid::FaceAttachmentAccessor<AInt> aaBMFACE;
			if(paFaceBoundaryMarker != NULL)
			{
				out << grid.num_faces() << " 1" << endl;
				aaBMFACE.access(grid, *paFaceBoundaryMarker);
			}
			else
				out << grid.num_faces() << " 0" << endl;

			for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter)
			{
				Face* f = *iter;
				out << f->num_vertices();

				for(uint i = 0; i < f->num_vertices(); ++i)
					out << " " << aaIntVRT[f->vertex(i)];

				if(paFaceBoundaryMarker != NULL)
					out << " " << aaBMFACE[f];

				out << endl;
			}
		}
		out.close();
	}

//	write tetrahedrons
	if(elemsFilename != NULL)
	{
		ofstream out(elemsFilename);
		if(out)
		{
			Grid::VolumeAttachmentAccessor<AInt> aaElementAttributeVOL;
			if(paElementAttribute != NULL)
			{
				aaElementAttributeVOL.access(grid, *paElementAttribute);
				out << grid.num<Tetrahedron>() << " 4 1" << endl;
			}
			else
				out << grid.num<Tetrahedron>() << " 4 0" << endl;

			int counter = 0;
			for(TetrahedronIterator iter = grid.begin<Tetrahedron>();
									iter != grid.end<Tetrahedron>(); iter++, counter++)
			{
				Tetrahedron* tet = *iter;
				out << counter << " " <<	aaIntVRT[tet->vertex(0)] << " " <<
											aaIntVRT[tet->vertex(1)] << " " <<
											aaIntVRT[tet->vertex(2)] << " " <<
											aaIntVRT[tet->vertex(3)];

				if(paElementAttribute != NULL)
					out << " " << aaElementAttributeVOL[tet];

				out << endl;
			}
		}

		out.close();
	}

	return true;
}

}//	end of namespace
