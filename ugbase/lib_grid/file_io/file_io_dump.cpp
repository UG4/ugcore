//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#include <fstream>
#include "file_io_dump.h"
#include "common/util/loader/loader_util.h"
#include "lib_grid/lg_base.h"

using namespace std;

namespace ug
{

static bool ReadTriangles(Grid& grid, ifstream& in,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int& lineCount)
{
	char* BUFFER = new char[512];
	vector<string> paramVec;
	bool bSuccess = true;
	
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
		//	check if the line is a commentary
			string& str = paramVec[0];
			if(str.find('#') == 0){
			//	theres nothing to do in here. Proceed by reading the next line.
			}
		//	if the first parameter contains "end", we're done in here.
			else if(str.find("end") == 0){
				break;
			}
		//	read the triangle
			else if(paramVec.size() == 9){
			//	this line contains the coordinates of a triangle.
			//	create vertices and assign the coordinates.
				Vertex* v[3];
				for(int i = 0; i < 3; ++i)
				{
					v[i] = *grid.create<Vertex>();
					aaPos[v[i]].x = atof(paramVec[i*3].c_str());
					aaPos[v[i]].y = atof(paramVec[i*3 + 1].c_str());
					aaPos[v[i]].z = atof(paramVec[i*3 + 2].c_str());
				}
			//	create a triangle
				grid.create<Triangle>(TriangleDescriptor(v[0], v[1], v[2]));
			}
			else
			{
			//	this line can't be interpreted correctly:
				bSuccess = false;
				break;
			}
		}
		lineCount++;
	}

	delete[] BUFFER;
	
	return bSuccess;
}

static bool ReadTetrahedrons(Grid& grid, ifstream& in,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int& lineCount)
{
	char* BUFFER = new char[512];
	vector<string> paramVec;
	bool bSuccess = true;
	
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
		//	check if the line is a commentary
			string& str = paramVec[0];
			if(str.find('#') == 0){
			//	theres nothing to do in here. Proceed by reading the next line.
			}
		//	if the first parameter contains "end", we're done in here.
			else if(str.find("end") == 0){
				break;
			}
		//	read the triangle
			else if(paramVec.size() == 12){
			//	this line contains the coordinates of a triangle.
			//	create vertices and assign the coordinates.
				Vertex* v[4];
				for(int i = 0; i < 4; ++i)
				{
					v[i] = *grid.create<Vertex>();
					aaPos[v[i]].x = atof(paramVec[i*3].c_str());
					aaPos[v[i]].y = atof(paramVec[i*3 + 1].c_str());
					aaPos[v[i]].z = atof(paramVec[i*3 + 2].c_str());
				}
			//	create a triangle
				grid.create<Tetrahedron>(TetrahedronDescriptor(v[0], v[1], v[2], v[3]));
			}
			else
			{
			//	this line can't be interpreted correctly:
				bSuccess = false;
				break;
			}
		}
		lineCount++;
	}

	delete[] BUFFER;
	
	return bSuccess;
}

bool LoadGridFromDUMP(Grid& grid, const char* filename,
					ISubsetHandler* pSH, AVector3& aPos)
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

	char* BUFFER = new char[512];
	vector<string> paramVec;
	bool bSuccess = true;
	int lineCount = 0;
	
	int oldDefSI = -1;
	
	if(pSH){
		oldDefSI = pSH->get_default_subset_index();
		pSH->set_default_subset_index(0);
	}
	
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
		//	check if the line is a commentary
			string& str = paramVec[0];
			if(str.find('#') == 0)
			{
			//	theres nothing to do in here. Proceed by reading the next line.
			}
		//	check whether a new block starts
			else if(str.find("begin") == 0)
			{
				if(paramVec.size() < 2){
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  missing block type specifier in line " << lineCount << endl);
					continue;
				}
				
				if(paramVec.size() > 2){
					if(pSH)
						pSH->set_default_subset_index(atoi(paramVec[2].c_str()));
				}
				
			//	check second paramenter
				string& type = paramVec[1];
				
				if(type.find("triangles") == 0){
					if(!ReadTriangles(grid, in, aaPos, lineCount)){
						LOG("  PROBLEM while reading from " << filename << ":\n");
						LOG("  ReadTriangles failed in line " << lineCount << endl);
						bSuccess = false;
					}
				}
				else if(type.find("tetrahedrons") == 0){
					if(!ReadTetrahedrons(grid, in ,aaPos, lineCount)){
						LOG("  PROBLEM while reading from " << filename << ":\n");
						LOG("  ReadTetrahedrons failed in line " << lineCount << endl);
						bSuccess = false;
					}
				}
				else{
					LOG("  PROBLEM while reading from " << filename << ":\n");
					LOG("  unknown block type specifier in line " << lineCount << endl);				
					bSuccess = false;
				}
			}
		}
		lineCount++;
	}

	delete[] BUFFER;

	if(pSH){
		pSH->set_default_subset_index(oldDefSI);
	}
	
	return bSuccess;
}

}
