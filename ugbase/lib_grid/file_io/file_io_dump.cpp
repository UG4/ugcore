//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#include <fstream>
#include "file_io_dump.h"
#include "common/util/loader/loader_util.h"

using namespace std;

namespace ug
{

bool LoadGridFromDUMP(Grid& grid, const char* filename, AVector3& aPos)
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

	int lineCount = 0;
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
			else if(paramVec.size() == 9)
			{
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
				LOG("PROBLEM while reading from " << filename << ":\n");
				LOG("  line " << lineCount << " contains wrong number of parameters:\n");
			}
		}
		lineCount++;
	}

	delete[] BUFFER;

	return true;
}

}
