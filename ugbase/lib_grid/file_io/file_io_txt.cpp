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

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "file_io_txt.h"
#include "../lg_base.h"

using namespace std;

namespace ug
{

bool LoadGridFromTXT(Grid& grid, const char* filename, AVector3& aPos)
{
	ifstream in(filename);

	if(!in)
		return false;

	//grid.clear();

	int numVrts, numElems;

	in >> numVrts;
	in >> numElems;

//	create points
//	store pointers to the vertices on the fly in a vector.
	vector<Vertex*>	vVrts;
	vector<int> vVrtIds;
	vVrts.resize(numVrts); vVrtIds.resize(numVrts);

	for(int i = 0; i < numVrts; ++i)
		vVrts[i] = *grid.create<RegularVertex>();

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	read the points
	{
		int i = 0;
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter, ++i)
		{
			in >> vVrtIds[i];
			in >> aaPos[*iter].x();
			in >> aaPos[*iter].y();
			in >> aaPos[*iter].z();
		}
		if (in.fail())
			UG_THROW ("LoadGridFromTXT: Failed to read vertices from '" << filename << "'");
	}

//	read the triangles and the quadrilaterals
	{
		string e_line;
		
		for(int i = 0; i < numElems; ++i)
		{
			int Index, vrt_id_1, vrt_id_2, vrt_id_3, vrt_id_4;
			vector<int>::iterator i1, i2, i3, i4;
			
			do
			{
				getline(in, e_line);
			}
			while (e_line.empty());
			stringstream ss(e_line, ios_base::in);
			
			ss >> Index;
			ss >> vrt_id_1;
			ss >> vrt_id_2;
			ss >> vrt_id_3;
			
			if (ss.fail())
				UG_THROW ("LoadGridFromTXT: Failed to read grid elements from '" << filename << "'");
			
			i1 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_1);
			i2 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_2);
			i3 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_3);
			
			if (i1 == vVrtIds.end() || i2 == vVrtIds.end() || i3 == vVrtIds.end())
				UG_THROW ("LoadGridFromTXT: Wrong index of a vertex of element " << Index << "in '" << filename << "'");
			
			ss >> vrt_id_4;
			if (ss.fail ())
				grid.create<Triangle>
					(TriangleDescriptor
						(vVrts[i1 - vVrtIds.begin()], vVrts[i2 - vVrtIds.begin()], vVrts[i3 - vVrtIds.begin()]));
			else
			{
				i4 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_4);
				if (i4 == vVrtIds.end())
					UG_THROW ("LoadGridFromTXT: Wrong index of the last vertex of quadrilateral " << Index << "in '" << filename << "'");
				grid.create<Quadrilateral>
					(QuadrilateralDescriptor
						(vVrts[i1 - vVrtIds.begin()], vVrts[i2 - vVrtIds.begin()], vVrts[i3 - vVrtIds.begin()], vVrts[i4 - vVrtIds.begin()]));
			}
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
	out << grid.num_vertices() << " " << grid.num<Triangle>() + grid.num<Quadrilateral>() << endl;

//	write the vertices
//	store in each vertex at which position it has been written to the file.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);

	{
		int counter = 0;

		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++)
		{
			out << counter << " " << 	aaPos[*iter].x() << " " <<
										aaPos[*iter].y() << " " <<
										aaPos[*iter].z() << endl;

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

		for(QuadrilateralIterator iter = grid.begin<Quadrilateral>(); iter != grid.end<Quadrilateral>(); ++iter, ++counter)
		{
			Quadrilateral* t = *iter;
			out << counter << " " 	<< aaInt[t->vertex(0)] << " "
									<< aaInt[t->vertex(1)] << " "
									<< aaInt[t->vertex(2)] << " "
									<< aaInt[t->vertex(3)] << endl;
		}
	}

	grid.detach_from_vertices(aInt);
	return true;
}

}//	end of namespace
