// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include <fstream>
#include <cmath>
#include "file_io_tikz.h"

using namespace std;

namespace ug{

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

bool ExportGridToTIKZ(Grid& grid, const char* filename, const ISubsetHandler* psh,
					  APosition aPos, TikzExportDesc desc)	
{
	ofstream out(filename);
	if(!out){
		UG_LOG("Couldn't load file " << filename << "\n");
		return false;
	}

	if(!psh){
		UG_LOG("A subset handler is required to write a tikz file!\n");
		return false;
	}

	number sml = max(SMALL, desc.smallestVal);

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

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
		out << "%subset styles for subset " << si << " (" << psh->subset_info(si).name << ")" << endl;
		out << "\\tikzset{vertex" << si << "/.style={vertexBase}}\n";
		out << "\\tikzset{edge" << si << "/.style={edgeBase}}\n";
		out << "\\tikzset{face" << si << "/.style={faceBase}}\n";
	}
	
	out << endl;

	out << "\\begin{scope}" << endl;
	for(int si = 0; si < psh->num_subsets(); ++si){
	//	draw all faces which are in a subset
		for(FaceIterator iter = grid.begin<Face>();
			iter != grid.end<Face>(); ++iter)
		{
			Face* f = *iter;
			if(psh->get_subset_index(f) != si)	// THIS IS NASTY! One should iterate over subset-elements directly instead
				continue;
			out << "\\filldraw[face" << si << "]";
			for(size_t i = 0; i < f->num_vertices(); ++i){
				Vertex* vrt = f->vertex(i);
				vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
				out << " (" << p.x() << "cm, " << p.y() << "cm) --";
			}
			out << " cycle;" << endl;
		}
	}

	for(int si = 0; si < psh->num_subsets(); ++si){
	//	draw all edges which are in a subset
		for(EdgeIterator iter = grid.begin<Edge>();
			iter != grid.end<Edge>(); ++iter)
		{
			Edge* e = *iter;
			if(psh->get_subset_index(e) != si)	// THIS IS NASTY! One should iterate over subset-elements directly instead
				continue;
			out << "\\draw[edge" << si << "]";
			for(size_t i = 0; i < e->num_vertices(); ++i){
				Vertex* vrt = e->vertex(i);
				vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
				out << " (" << p.x() << "cm, " << p.y() << "cm)";
				if(i+1 != e->num_vertices())
					out << " --";
			}
			out << ";" << endl;
		}
	}

	for(int si = 0; si < psh->num_subsets(); ++si){
	//	draw all nodes which are in a subset
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter)
		{
			Vertex* vrt = *iter;
			if(psh->get_subset_index(vrt) != si)	// THIS IS NASTY! One should iterate over subset-elements directly instead
				continue;
			vector2 p = trunk(vector2(aaPos[vrt].x(), aaPos[vrt].y()), sml);
			out << "\\node[vertex" << si << "] at (" << p.x() << "cm, " << p.y() << "cm) {};" << endl;
		}
	}
	out << "\\end{scope}" << endl;

	out << "\\end{tikzpicture}" << endl;
	out << "\\end{document}" << endl;

	out.close();
	return true;
}

}//	end of namespace
