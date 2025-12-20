/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

#include "file_io_tikz.h"
#include "common/util/string_util.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/bounding_box_util.h"
//ø #include "lib_grid/iterators/lg_for_each.h"

using namespace std;

namespace ug {

TikzExportDesc::
TikzExportDesc() :
	vrtRadius(0.1),
	vrtRimWidth(0.02),
	vrtColor(0.65, 0.65, 0.65),
	edgeWidth(0.025),
	edgeColor(0, 0, 0),
	faceRimWidth(0.015),
	faceColor(0.85, 0.85, 0.85),
	smallestVal(0.001)
{}

static vector2 trunk(const vector2& v, number smallestVal)
{
	vector2 t;
	t.x() = smallestVal * round(v.x() / smallestVal);
	t.y() = smallestVal * round(v.y() / smallestVal);
	return t;
}


struct TIKZElem{
	GridObject* elem;
	float zmin;
	float zmax;
	int elemId;
	int subsetId;

	template <typename vector_t>
	TIKZElem(GridObject* _e, int _si, const AABox<vector_t>& bbox) :
		elem(_e),
		elemId(_e->base_object_id()),
		subsetId(_si)
		{
			if(vector_t::Size == 3){
				zmin = bbox.min[2];
				zmax = bbox.max[2];
			}
			else
				zmin = zmax = 0;
		}

//	if x < y, x will be rendered before y
	bool operator <(const TIKZElem& e) const{
		if(zmin - e.zmin > SMALL)			return false;
		else if(e.zmin - zmin > SMALL)	return true;
		
		if(e.elemId > elemId)		return false;
		else if(e.elemId < elemId)	return true;

		if(zmax - e.zmax > SMALL)			return false;
		else if(e.zmax - zmax > SMALL)	return true;

		if(e.subsetId < subsetId)		return false;
		else if(e.subsetId > subsetId)	return true;

		return false;	
	}
};



bool ExportGridToTIKZ(Grid& grid, const char* filename, const ISubsetHandler* psh,
					  APosition aPos, TikzExportDesc desc)	
{
	UG_COND_THROW(!psh, "A subset handler is required to write a tikz file!\n");

	ofstream out(filename);
	UG_COND_THROW(!out, "Couldn't load file " << filename << "\n");


	vector<string>	subsetIdentifyer;
	vector<bool>	subsetNameIsDuplicate;
	for(int si = 0; si < psh->num_subsets(); ++si){
		string sname = psh->get_subset_name(si);
		for(size_t i = 0; i < sname.size(); ++i){
			char c = sname[i];
			if(!(isalnum(c) || c == '_'))
				sname[i] = '_';
		}

	//	check if the string is already contained in subsetIdentifyers and mark it
	//	as a duplicate if that's the case
		bool isDuplicate = false;
		for(size_t i = 0; i < subsetIdentifyer.size(); ++i){
			if(subsetIdentifyer[i] == sname){
				isDuplicate = true;
				break;
			}
		}
		subsetIdentifyer.push_back(sname);
		subsetNameIsDuplicate.push_back(isDuplicate);
	}

	number sml = max(SMALL, desc.smallestVal);

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

	out << "% This tex-file was exported from ProMesh (www.promesh3d.com)" << endl << endl;
	out << "% Call 'pdflatex' with this script to generate a .pdf file from it." << endl << endl;
	out << "% By placing a file 'custom_promesh_style.tex' in the same folder as this" << endl;
	out << "% script, you may override all local style definitions with custom definitions." << endl << endl;

	out << "\\documentclass[tikz]{standalone}" << endl;
	out << "\\begin{document}" << endl;
	out << "\\begin{tikzpicture}" << endl;

	out << endl;
	out << "\\definecolor{vertexColor}{rgb}{" << desc.vrtColor.x() << ", "
		<< desc.vrtColor.y() << ", " << desc.vrtColor.z() << "}" << endl;
	out << "\\definecolor{edgeColor}{rgb}{" << desc.edgeColor.x() << ", "
		<< desc.edgeColor.y() << ", " << desc.edgeColor.z() << "}" << endl;
	out << "\\definecolor{faceColor}{rgb}{" << desc.faceColor.x() << ", "
		<< desc.faceColor.y() << ", " << desc.faceColor.z() << "}" << endl;
	out << endl;

	out << "\\tikzset{vertexBase/.style={circle,draw=black,fill=vertexColor,line width=" << desc.vrtRimWidth << "cm,\n"
		<< "                        inner sep=0,minimum size=" << desc.vrtRadius << "cm}}\n";
	out << "\\tikzset{edgeBase/.style={draw=edgeColor, line width=" << desc.edgeWidth << "cm}}\n";
	out << "\\tikzset{faceBase/.style={draw=black, fill=faceColor, line width=" << desc.faceRimWidth <<"cm}}\n";
	
	out << endl;

	for(int si = 0; si < psh->num_subsets(); ++si){
		if(!subsetNameIsDuplicate[si]){
			const std::string& subsetName = subsetIdentifyer[si];
			out << "% '" << subsetName << "' subset styles" <<endl;
			out << "\\tikzset{vertex_" << subsetName << "/.style={vertexBase}}\n";
			out << "\\tikzset{edge_" << subsetName << "/.style={edgeBase}}\n";
			out << "\\tikzset{face_" << subsetName << "/.style={faceBase}}\n";
		}
	}
	
	out << endl;
	out << "\\InputIfFileExists{./custom_promesh_style.tex}{}{}" << endl;
	out << endl;


//	create a vector which contains all objects that shall be rendered and sort it
	vector<TIKZElem>	elems;
	for(auto _feI = grid.begin<Vertex>(); _feI != grid.end<Vertex>(); ++_feI){
		Vertex* v = *_feI;
		if(psh->get_subset_index(v) != -1)
			elems.emplace_back(v, psh->get_subset_index(v),
									 CalculateBoundingBox(v, aaPos));
	}

	for(auto _feI = grid.begin<Edge>(); _feI != grid.end<Edge>(); ++_feI){
		Edge* e = *_feI;
		if(psh->get_subset_index(e) != -1)
			elems.emplace_back(e, psh->get_subset_index(e),
									 CalculateBoundingBox(e, aaPos));
	}

	for(auto _feI = grid.begin<Face>(); _feI != grid.end<Face>(); ++_feI){
		Face* f = *_feI;
		if(psh->get_subset_index(f) != -1)
			elems.emplace_back(f, psh->get_subset_index(f),
									 CalculateBoundingBox(f, aaPos));
	}

	sort(elems.begin(), elems.end());

	out << "\\begin{scope}" << endl;
	for(size_t _vfeI = 0; _vfeI < elems.size(); ++_vfeI){
		TIKZElem& tikzElem = elems[_vfeI];
		const std::string& subsetName = subsetIdentifyer[tikzElem.subsetId];
		switch(tikzElem.elemId){
			case GridBaseObjectId::VERTEX:{
				auto vrt = static_cast<Vertex*>(tikzElem.elem);
				vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
				out << "\\node[vertex_" << subsetName << "] at (" << p.x() << "cm, " << p.y() << "cm) {};" << endl;
			}break;

			case GridBaseObjectId::EDGE:{
				auto* e = static_cast<Edge*>(tikzElem.elem);
				out << "\\draw[edge_" << subsetName << "]";
				for(size_t i = 0; i < e->num_vertices(); ++i){
					Vertex* vrt = e->vertex(i);
					vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
					out << " (" << p.x() << "cm, " << p.y() << "cm)";
					if(i+1 != e->num_vertices())
						out << " --";
				}
				out << ";" << endl;
			}break;

			case GridBaseObjectId::FACE:{
				auto f = static_cast<Face*>(tikzElem.elem);
				out << "\\filldraw[face_" << subsetName << "]";
				for(size_t i = 0; i < f->num_vertices(); ++i){
					Vertex* vrt = f->vertex(i);
					vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
					out << " (" << p.x() << "cm, " << p.y() << "cm) --";
				}
				out << " cycle;" << endl;
			}break;

		}

	}
	

	// for(int si = 0; si < psh->num_subsets(); ++si){
	// //	draw all faces which are in a subset
	// 	for(FaceIterator iter = grid.begin<Face>();
	// 		iter != grid.end<Face>(); ++iter)
	// 	{
	// 		Face* f = *iter;
	// 		if(psh->get_subset_index(f) != si)	// THIS IS NASTY! One should iterate over subset-elements directly instead
	// 			continue;
	// 		out << "\\filldraw[face" << si << "]";
	// 		for(size_t i = 0; i < f->num_vertices(); ++i){
	// 			Vertex* vrt = f->vertex(i);
	// 			vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
	// 			out << " (" << p.x() << "cm, " << p.y() << "cm) --";
	// 		}
	// 		out << " cycle;" << endl;
	// 	}
	// }

	// for(int si = 0; si < psh->num_subsets(); ++si){
	// //	draw all edges which are in a subset
	// 	for(EdgeIterator iter = grid.begin<Edge>();
	// 		iter != grid.end<Edge>(); ++iter)
	// 	{
	// 		Edge* e = *iter;
	// 		if(psh->get_subset_index(e) != si)	// THIS IS NASTY! One should iterate over subset-elements directly instead
	// 			continue;
	// 		out << "\\draw[edge" << si << "]";
	// 		for(size_t i = 0; i < e->num_vertices(); ++i){
	// 			Vertex* vrt = e->vertex(i);
	// 			vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
	// 			out << " (" << p.x() << "cm, " << p.y() << "cm)";
	// 			if(i+1 != e->num_vertices())
	// 				out << " --";
	// 		}
	// 		out << ";" << endl;
	// 	}
	// }

	// for(int si = 0; si < psh->num_subsets(); ++si){
	// //	draw all nodes which are in a subset
	// 	for(VertexIterator iter = grid.begin<Vertex>();
	// 		iter != grid.end<Vertex>(); ++iter)
	// 	{
	// 		Vertex* vrt = *iter;
	// 		if(psh->get_subset_index(vrt) != si)	// THIS IS NASTY! One should iterate over subset-elements directly instead
	// 			continue;
	// 		vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
	// 		out << "\\node[vertex" << si << "] at (" << p.x() << "cm, " << p.y() << "cm) {};" << endl;
	// 	}
	// }
	out << "\\end{scope}" << endl;

	out << "\\end{tikzpicture}" << endl;
	out << "\\end{document}" << endl;

	out.close();
	return true;
}

}//	end of namespace
