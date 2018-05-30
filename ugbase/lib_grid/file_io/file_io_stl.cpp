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

#include <fstream>
#include <algorithm>
#include "file_io_stl.h"
#include "common/util/loader/loader_util.h"
#include "common/util/string_util.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"

#include "stl_reader.h"

using namespace std;

namespace ug
{

bool LoadGridFromSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH,
					AVector3& aPos,
					AVector3& aNormFACE)
{
	vector<number> coords, normals;
	vector<unsigned int> tris, solids;

	if(!ReadStlFile(filename, coords, normals, tris, solids))
		return false;

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	
	const size_t numVrts = coords.size() / 3;
	grid.reserve<Vertex>(grid.num<Vertex>() + numVrts);
	vector<Vertex*> vrts(numVrts);

	for(size_t i = 0; i < numVrts; ++i){
		const size_t ci = i * 3;
		Vertex* v = *grid.create<RegularVertex>();
		vector3& p = aaPos[v];
		for(size_t j = 0; j < 3; ++j)
			p[j] = coords[ci + j];
		vrts[i] = v;
	}

	const size_t numTris = tris.size() / 3;
	grid.reserve<Face>(grid.num<Face>() + numTris);

	Grid::FaceAttachmentAccessor<AVector3> aaNormFACE;
	if(grid.has_face_attachment(aNormFACE))
		aaNormFACE.access(grid, aNormFACE);

	size_t curSolid = 0;
	for(size_t i = 0; i < numTris; ++i){
		const size_t ti = i*3;
		const unsigned int* t = &tris[ti];
		Face* f = *grid.create<Triangle>(
						TriangleDescriptor(vrts[t[0]], vrts[t[1]], vrts[t[2]]));
	
		if(pSH){	
			while((curSolid + 1 < solids.size()) && (i >= solids[curSolid + 1]))
				++curSolid;
			pSH->assign_subset(f, (int)curSolid);
		}
	
		if(aaNormFACE.valid()){
			vector3& n = aaNormFACE[f];
			for(size_t j = 0; j < 3; ++j)
				n[j] = normals[ti + j];
		}
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
