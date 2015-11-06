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
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "file_io_ncdf.h"

using namespace std;

namespace ug
{

bool SaveGridToNCDF(Grid& grid, const char* filename,
					ISubsetHandler* pSH,
					APosition aPos)
{
//	open the file to which we'll write
	ofstream out(filename);
	if(!out)
		return false;

//	access subset-handler
	if(!pSH){
		UG_LOG("ERROR: SubsetHandler required for NCDF (Exodus) export.\n");
		return false;
	}

	ISubsetHandler& sh = *pSH;

//	access the position attachment
	if(!grid.has_vertex_attachment(aPos))
		return false;

	if(grid.num_volumes() != grid.num<Tetrahedron>()){
		UG_LOG("Saving grid to EXODUS-format.\n"
				<< "  WARNING: only Tetrahedrons exported in the moment.\n");
	}

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	assign vertex indices
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
	AssignIndices(grid.vertices_begin(), grid.vertices_end(), aaInt, 1);

//	write exodus header
	int var4 = 4;
	out << "netcdf object {" << endl;
	out << "dimensions:" << endl;
	out << "\tlen_string = 33 ;" << endl;
	out << "\tlen_line = 81 ;" << endl;
	out << "\tfour = " << var4 << " ;" << endl;
	out << "\ttime_step = UNLIMITED ;" << endl;
	out << "\tnum_dim = 3 ;" << endl;
	out << "\tnum_nodes = " << grid.num_vertices() << " ;" << endl;
	out << "\tnum_elem = " << grid.num<Tetrahedron>() << " ;" << endl;
	out << "\tnum_el_blk = " << sh.num_subsets() << " ;" << endl;

	for(int i = 0; i < sh.num_subsets(); ++i){
		GridObjectCollection goc = sh.get_grid_objects_in_subset(i);
		out << "\tnum_el_in_blk" << i+1 << " = " << goc.num<Tetrahedron>() << " ;" << endl;
		out << "\tnum_nod_per_el" << i+1 << " = 4 ;" << endl;
	}

	out << "\tnum_qa_rec = 1 ;" << endl;

//	write variables
	out << endl;
	out << "variables:" << endl;
	out << "\tdouble time_whole(time_step) ;" << endl;
	out << "\tint eb_status(num_el_blk) ;" << endl;
	out << "\tint eb_prop1(num_el_blk) ;" << endl;
	out << "\t\teb_prop1:name = \"ID\" ;" << endl;

	for(int i = 0; i < sh.num_subsets(); ++i){
		out << "\tint connect" << i+1
			<< "(num_el_in_blk" << i+1
			<< ", num_nod_per_el" << i+1 << ") ;" << endl;
		out << "\t\tconnect" << i+1 << ":elem_type = \"TETRA4\" ;" << endl;
	}

	out << "\tdouble coord(num_dim, num_nodes) ;" << endl;
	out << "\tchar qa_records(num_qa_rec, four, len_string) ;" << endl;
	out << "\tchar coor_names(num_dim, len_string) ;" << endl;
	out << "\tint elem_map(num_elem) ;" << endl;
	out << "\tint elem_num_map(num_elem) ;" << endl;
	out << "\tint node_num_map(num_nodes) ;" << endl;

	out << "// global attributes:" << endl;
	out << "\t:api_version = 4.01f ;" << endl;
	out << "\t:version = 3.01f ;" << endl;
	out << "\t:floating_point_word_size = " << sizeof(number) << " ;" << endl;
	out << "\t:file_size = 0 ;" << endl;
	out << "\t:title = \"Exported from lib_grid.\" ;" << endl;

//	write data
	out << endl;
	out << "data:" << endl;

	out << "eb_status = ";
	for(int i = 0; i < sh.num_subsets(); ++i){
		out << "1";
		if(i + 1 < sh.num_subsets())
			out << ", ";
	}
	out << " ;" << endl << endl;

	out << "eb_prop1 = ";
	for(int i = 0; i < sh.num_subsets(); ++i){
		out << i+1;
		if(i + 1 < sh.num_subsets())
			out << ", ";
	}
	out << " ;" << endl << endl;

//	write elements
	for(int i = 0; i < sh.num_subsets(); ++i){
	//	get the goc for this subset
		GridObjectCollection goc = sh.get_grid_objects_in_subset(i);

		out << "connect" << i+1 << " =" << endl;
		for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
			for(TetrahedronIterator iter = goc.begin<Tetrahedron>(lvl);
				iter != goc.end<Tetrahedron>(lvl); ++iter)
			{
			//	last comma in each row
				if(iter != goc.begin<Tetrahedron>(lvl))
					out << "," << endl;

				Tetrahedron* tet = *iter;
				out << "  ";
				out << aaInt[tet->vertex(0)] << ", ";
				out << aaInt[tet->vertex(1)] << ", ";
				out << aaInt[tet->vertex(2)] << ", ";
				out << aaInt[tet->vertex(3)];
			}
		}
		out << " ;" << endl << endl;
	}

//	write coords
//	x
	out << "coord =" << endl << " ";
	size_t endlCounter = 1;
	for(VertexIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter, ++endlCounter)
	{
		if(endlCounter > 5){
			endlCounter = 1;
			out << endl << " ";
		}

		out << " " << aaPos[*iter].x() << ",";
	}

//	y
	for(VertexIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter, ++endlCounter)
	{
		if(endlCounter > 5){
			endlCounter = 1;
			out << endl << " ";
		}
		out << " " << aaPos[*iter].y() << ",";
	}

//	z
	for(VertexIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter, ++endlCounter)
	{
		if(iter != grid.vertices_begin()){
			out << ",";
			if(endlCounter > 5){
				endlCounter = 1;
				out << endl << " ";
			}
		}

		out << " " << aaPos[*iter].z();
	}
	out << " ;" << endl << endl;

//	write qa_records
	out << "qa_records =" << endl;
	out << "  \"lib_grid\"," << endl;
	out << "  \"...\"," << endl;
	out << "  \"...\"," << endl;
	out << "  \"...\" ;" << endl;

//	write coor_names
	out << endl;
	out << "coor_names =" << endl;
	out << "  \"x\"," << endl;
	out << "  \"y\"," << endl;
	out << "  \"z\" ;" << endl;

//	write elem_map
	out << endl;
	out << "elem_map = ";
	endlCounter = 1;
	for(size_t i = 0; i < grid.num<Tetrahedron>(); ++i, ++endlCounter)
	{
		out << i+1;
		if(i+1 < grid.num<Tetrahedron>()){
			out << ", ";

			if(endlCounter > 4){
				endlCounter = 0;
				out << endl << "  ";
			}
		}
	}
	out << " ;" << endl;

//	write elem_num_map
	out << endl;
	out << "elem_num_map = ";
	endlCounter = 1;
	for(size_t i = 0; i < grid.num<Tetrahedron>(); ++i, ++endlCounter)
	{
		out << i+1;
		if(i+1 < grid.num<Tetrahedron>()){
			out << ", ";

			if(endlCounter > 4){
				endlCounter = 0;
				out << endl << "  ";
			}
		}
	}
	out << " ;" << endl;

//	write node_num_map
	out << endl;
	out << "node_num_map = ";
	endlCounter = 1;
	for(size_t i = 0; i < grid.num<Vertex>(); ++i, ++endlCounter)
	{
		out << i+1;
		if(i+1 < grid.num<Vertex>()){
			out << ", ";

			if(endlCounter > 4){
				endlCounter = 0;
				out << endl << "  ";
			}
		}
	}
	out << " ;" << endl;

//	close object
	out << "}" << endl;
//	remove vertex indices
	grid.detach_from_vertices(aInt);

//	done
	return true;
}

}//	end of namespace
