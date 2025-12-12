/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include "file_io_msh.h"

#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>

#include "../lg_base.h"

using namespace std;

namespace ug {

bool LoadGridFromMSH(Grid& grid, const char* filename,
					 ISubsetHandler* psh, AVector3& aPos)
{
//	the position attachment	
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	open the file
	ifstream in(filename);
	if(!in){
		UG_LOG("File not found: " << filename << "\n");
		return false;
	}
//	this buffer is used to find the commands
	string buffer;	
	
//	here we'll store all nodes, so that we can access them by index
	vector<Vertex*> vrts;
	
//	iterate through the lines
	while(!in.eof())
	{
	//	check whether any commands can be found
		in >> buffer;
	//	make all entries upper case
		transform(buffer.begin(), buffer.end(), buffer.begin(), (int(*)(int)) toupper);
		//UG_LOG(buffer);
		
	//	compare with known commands
		if(buffer.find("$NODES") != string::npos){
		//	read the nodes array
			int numNodes;
			in >> numNodes;
			if(in.fail()){
				UG_LOG("LoadGridFromMSH: bad format in $NODES - numNodes\n");
				return false;
			}
			
		//	iterate through the nodes. Create a new vertex for each
		//	and assign the position
			int ind;
			vector3 p;
			for(int i = 0; i < numNodes; ++i){
				in >> ind >> p.x() >> p.y() >> p.z();
				
				if(in.fail()){
					UG_LOG("LoadGridFromMSH: bad format in $NODES\n");
					return false;
				}
				
				RegularVertex* vrt = *grid.create<RegularVertex>();
				aaPos[vrt] = p;
				vrts.push_back(vrt);
			}
		}
		else if(buffer.find("$ELEMENTS") != string::npos)
		{
			int numVrts = (int)vrts.size();
		//	read the number of elements
			int numElems;
			in >> numElems;
			if(in.fail()){
				UG_LOG("LoadGridFromMSH: bad format in $ELEMENTS - numElems\n");
				return false;
			}
			
		//	iterate through the entries and create the elements
			int ind, type, numTags;
			vector<int> tags;
			
			for(int i = 0; i < numElems; ++i){
				in >> ind >> type >> numTags;
				if(in.fail()){
					UG_LOG("LoadGridFromMSH: bad format in $ELEMENTS\n");
					return false;
				}
				
			//	read the tags (tags are ignored in the moment)
				tags.clear();
				for(int j = 0; j < numTags; ++j){
					int tag;
					in >> tag;
					tags.push_back(tag);
				}
				
				if(in.fail()){
					UG_LOG("LoadGridFromMSH: bad tags in $ELEMENTS\n");
					return false;
				}
				
				Face* f = nullptr;
				switch(type){
					case 2://triangles
					{
						int i1, i2, i3;
						in >> i1 >> i2 >> i3;
						if(in.fail()){
							UG_LOG("LoadGridFromMSH: bad format in triangle\n");
							return false;
						}
						
					//	make sure that indices are fine
						if((i1 < 0 || i1 >= numVrts) ||
						   (i2 < 0 || i2 >= numVrts) ||
						   (i3 < 0 || i3 >= numVrts))
						{
							UG_LOG("LoadGridFromMSH: bad vertex indices in triangle. ignoring triangle "
									<< grid.num<Triangle>() << ".\n");
							continue;
						}
						else{
							f = *grid.create<Triangle>(TriangleDescriptor(vrts[i1], vrts[i2], vrts[i3]));
						}
						
					}break;
					
					default:
					{
						UG_LOG("ERROR in LoadGridFromMSH: element type " << type << " not supported. aborting...\n");
						return false;
					}
				}//	end of switch
			//	assign subset index from tag
				if(psh){
					if(numTags < 3){
					// since there are no tags, we'll assign all elements to subset 0
						psh->assign_subset(f, 0);
					}
					else{
					//	tag 3 holds the mesh-partition. assign the index
						psh->assign_subset(f, tags[2]);
					}
				}
			}
		}
	}

//	done.
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
*/
}//	end of namespace
