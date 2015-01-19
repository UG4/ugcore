//	created by Sebastian Reiter
//	y10 m11 d11
//	s.b.reiter@googlemail.com

#include <fstream>
#include <algorithm>
#include "file_io_stl.h"
#include "common/util/loader/loader_util.h"
#include "common/util/string_util.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"

using namespace std;

namespace ug
{

bool STLFileHasASCIIFormat(const char* filename)
{
	ifstream in(filename);
	UG_COND_THROW(!in, "Couldn't open file " << filename);

	string firstWord;
	in >> firstWord;
	transform(firstWord.begin(), firstWord.end(), firstWord.begin(), ::tolower);

	return firstWord.compare("solid") == 0;
}


bool LoadGridFromSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH,
					AVector3& aPos,
					AVector3& aNormFACE)
{
	if(STLFileHasASCIIFormat(filename)){
	//	load ascii-stl
		return LoadGridFromSTL_ASCII(grid, filename, pSH, aPos, aNormFACE);
	}
	else{
	//	load binary-stl
		return LoadGridFromSTL_BINARY(grid, filename, pSH, aPos, aNormFACE);
	}
}

bool LoadGridFromSTL_ASCII(Grid& grid, const char* filename,
						   ISubsetHandler* pSH,
						   AVector3& aPos,
						   AVector3& aNormFACE)
{
	UG_LOG("Loading STL with ASCII format\n");

//	open the file
	ifstream in(filename);
	UG_COND_THROW(!in, "Couldn't open file " << filename);

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	
	Grid::FaceAttachmentAccessor<AVector3> aaNormFACE;
	if(grid.has_face_attachment(aNormFACE))
		aaNormFACE.access(grid, aNormFACE);

	string buffer;
	vector<string> paramVec;
	bool bSuccess = true;
	int lineCount = 1;
	vector3 n(0.0, 0.0, 0.0);
	vector<vector3> positions;
	int newSubsetIndex = -1;
		
	while(!in.eof())
	{
	//	read the line
		getline(in, buffer);

	//	split the parameters
	//	todo: this causes problems with tabs or newlines which are used in between
	//		  coordinates or keywords. Use a more robust approach that searches
	//		  for keywords instead.
		TokenizeTrimString(buffer, paramVec, ' ');

	//	check if there are some at all
		if(!paramVec.empty())
		{
		//	check if the line is a comment
			string& str = paramVec[0];
			if(str.compare("vertex") == 0){
				if(paramVec.size() != 4){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  vertex not specified correctly in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
				
				vector3 v;
			//	read the position
				v.x() = atof(paramVec[1].c_str());
				v.y() = atof(paramVec[2].c_str());
				v.z() = atof(paramVec[3].c_str());
				
				positions.push_back(v);
			}
			else if(str.compare("facet") == 0)
			{
				if(paramVec.size() != 5){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  triangle not specified correctly in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
				
				if(paramVec[1].compare("normal") != 0){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  Missing 'normal' specifier in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
				
			//	read the normal
				n.x() = atof(paramVec[2].c_str());
				n.y() = atof(paramVec[3].c_str());
				n.z() = atof(paramVec[4].c_str());
			}
			else if(str.compare("outer") == 0){
				bool ok = false;
				if(paramVec.size() > 1){
					if(paramVec[1].compare("loop") == 0)
						ok = true;
				}
				
				if(!ok){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  expecting 'outer loop' in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
			}
			else if(str.compare("endfacet") == 0){
				if(positions.size() != 3){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  bad number of vertices specified for face in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
			//	create the vertices
				Vertex* v[3];
				for(size_t i = 0; i < 3; ++i){
					v[i] = *grid.create<RegularVertex>();
					aaPos[v[i]] = positions[i];
				}
			
			//	clear positions array for next iteration
				positions.clear();
				
			//	create the face
				Triangle* t = *grid.create<Triangle>(
								TriangleDescriptor(v[0], v[1], v[2]));
				
			//	assign the normal
				if(aaNormFACE.valid())
					aaNormFACE[t] = n;
					
			//	assign the subset
				if(pSH){
					pSH->assign_subset(t, newSubsetIndex);
				}	
			}
			else if(str.compare("solid") == 0){
			//	add a new subset
				if(pSH){
					newSubsetIndex = pSH->num_subsets();
					if(paramVec.size() > 1)
						pSH->subset_info(newSubsetIndex).name = paramVec[1];
				}
					
			}
		}
		lineCount++;
	}

	return bSuccess;
}

bool LoadGridFromSTL_BINARY(Grid& grid, const char* filename,
							ISubsetHandler* pSH,
							AVector3& aPos,
							AVector3& aNormFACE)
{
	UG_LOG("Loading STL with BINARY format\n");

	ifstream in(filename, ios::binary);
	UG_COND_THROW(!in, "Couldn't open file " << filename);

	char stl_header[80];
	in.read(stl_header, 80);
	UG_COND_THROW(!in, "Error while parsing binary stl header in file " << filename);

	uint32 numTris = 0;
	in.read((char*)&numTris, 4);
	UG_COND_THROW(!in, "Couldn't determine number of triangles in binary stl file " << filename);

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	
	Grid::FaceAttachmentAccessor<AVector3> aaNormFACE;
	if(grid.has_face_attachment(aNormFACE))
		aaNormFACE.access(grid, aNormFACE);

	int si = pSH->num_subsets();

	for(uint32 tri = 0; tri < numTris; ++tri){
		float d[12];
		in.read((char*)d, 12 * 4);
		UG_COND_THROW(!in, "Error while parsing trianlge in binary stl file " << filename);

		Vertex* v[3];
		for(int i = 0; i < 3; ++i){
			v[i] = *grid.create<RegularVertex>();
			float* di = d + (i + 1) * 3;
			aaPos[v[i]] = vector3(di[0], di[1], di[2]);
		}

		Triangle* f = *grid.create<Triangle>(TriangleDescriptor(v[0], v[1], v[2]));

		if(pSH)
			pSH->assign_subset(f, si);

		if(aaNormFACE.valid())
			aaNormFACE[f] = vector3(d[0], d[1], d[2]);

	//	read additional data
		char addData[2];
		in.read(addData, 2);
		UG_COND_THROW(!in, "Error while parsing additional triangle data in binary stl file " << filename);
	}

	return true;
}

bool SaveGridToSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH,
					AVector3& aPos)
{
	ofstream out(filename);
	if(!out){
		UG_LOG("Couldn't open file " << filename << " for writing\n");
		return false;
	}

	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("Specified vertex-position attachment missing!\n");
	}

	if(grid.num<Quadrilateral>() > 0){
		UG_LOG("WARNING: The specified grid contains quadrilaterals. "
				"Those won't be written to the stl file!\n");
	}

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

	const char* solidName = "stlFromUG4";
	out << "solid " << solidName << endl;

	for(Grid::traits<Triangle>::iterator iter = grid.begin<Triangle>();
		iter != grid.end<Triangle>(); ++iter)
	{
		Triangle* f = *iter;
		vector3 n;
		CalculateNormal(n, f, aaPos);
		out << "facet normal " << n.x() << " " << n.y() << " " << n.z() << endl;
		out << "outer loop" << endl;
		for(size_t i = 0; i < 3; ++i){
			vector3 v = aaPos[f->vertex(i)];
			out << "vertex " << v.x() << " " << v.y() << " " << v.z() << endl;
		}
		out << "endloop" << endl;
		out << "endfacet" << endl;
	}

	out << "endsolid " << solidName << endl;
	return true;
}

}
