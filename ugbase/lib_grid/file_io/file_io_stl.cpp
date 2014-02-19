//	created by Sebastian Reiter
//	y10 m11 d11
//	s.b.reiter@googlemail.com

#include <fstream>
#include "file_io_stl.h"
#include "common/util/loader/loader_util.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"

using namespace std;

namespace ug
{

bool LoadGridFromSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH,
					AVector3& aPos,
					AVector3& aNormFACE)
{
//	open the file
	ifstream in(filename);

	if(!in)
	{
		LOG(" file not found: " << filename << "\n");
		return false;
	}

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	
	Grid::FaceAttachmentAccessor<AVector3> aaNormFACE;
	if(grid.has_face_attachment(aNormFACE))
		aaNormFACE.access(grid, aNormFACE);

	char* BUFFER = new char[512];
	vector<string> paramVec;
	bool bSuccess = true;
	int lineCount = 1;
	vector3 n(0.0, 0.0, 0.0);
	vector<vector3> positions;
	int newSubsetIndex = -1;
		
	while(!in.eof())
	{
	//	read the line
		in.getline(BUFFER, 512);
		BUFFER[511] = 0;

	//	split the parameters
		split_parameters(paramVec, BUFFER, " \r");

	//	check if there are some at all
		if(!paramVec.empty())
		{
		//	check if the line is a comment
			string& str = paramVec[0];
			
			if(str.find("vertex") == 0){
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
			else if(str.find("facet") == 0)
			{
				if(paramVec.size() != 5){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  triangle not specified correctly in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
				
				if(paramVec[1].find("normal") != 0){
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
			else if(str.find("outer") == 0){
				bool ok = false;
				if(paramVec.size() > 1){
					if(paramVec[1].find("loop") == 0)
						ok = true;
				}
				
				if(!ok){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  expecting 'outer loop' in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
			}
			else if(str.find("endfacet") == 0){
				if(positions.size() != 3){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  bad number of vertices specified for face in line " << lineCount << endl);
					bSuccess = false;
					break;
				}
			//	create the vertices
				VertexBase* v[3];
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
			else if(str.find("solid") == 0){
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

	delete[] BUFFER;
	
	return bSuccess;
}

bool SaveGridToSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH,
					AVector3& aPos)
{
	ofstream out(filename);
	if(!out){
		UG_LOG("Couldn't load file " << filename << "\n");
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
