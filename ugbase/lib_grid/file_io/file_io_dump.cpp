/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "file_io_dump.h"

#include <fstream>

#include "common/util/loader/loader_util.h"
#include "lib_grid/lg_base.h"

using namespace std;

namespace ug {

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
				RegularVertex* v[3];
				for(int i = 0; i < 3; ++i)
				{
					v[i] = *grid.create<RegularVertex>();
					aaPos[v[i]].x() = atof(paramVec[i*3].c_str());
					aaPos[v[i]].y() = atof(paramVec[i*3 + 1].c_str());
					aaPos[v[i]].z() = atof(paramVec[i*3 + 2].c_str());
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
				RegularVertex* v[4];
				for(int i = 0; i < 4; ++i)
				{
					v[i] = *grid.create<RegularVertex>();
					aaPos[v[i]].x() = atof(paramVec[i*3].c_str());
					aaPos[v[i]].y() = atof(paramVec[i*3 + 1].c_str());
					aaPos[v[i]].z() = atof(paramVec[i*3 + 2].c_str());
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
		UG_LOG(" file not found: " << filename << "\n");
		return false;
	}

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor aaPos(grid, aPos);

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
					UG_LOG("  PROBLEM while reading from " << filename << ":\n");
					UG_LOG("  missing block type specifier in line " << lineCount << endl);
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
						UG_LOG("  PROBLEM while reading from " << filename << ":\n");
						UG_LOG("  ReadTriangles failed in line " << lineCount << endl);
						bSuccess = false;
					}
				}
				else if(type.find("tetrahedrons") == 0){
					if(!ReadTetrahedrons(grid, in ,aaPos, lineCount)){
						UG_LOG("  PROBLEM while reading from " << filename << ":\n");
						UG_LOG("  ReadTetrahedrons failed in line " << lineCount << endl);
						bSuccess = false;
					}
				}
				else{
					UG_LOG("  PROBLEM while reading from " << filename << ":\n");
					UG_LOG("  unknown block type specifier in line " << lineCount << endl);
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
