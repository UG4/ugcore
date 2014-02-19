// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)
 
#include "grid_objects_1d.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	RegularEdge
bool RegularEdge::refine(std::vector<Edge*>& vNewEdgesOut, Vertex* newVertex,
				Vertex** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new RegularEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new RegularEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new RegularEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new RegularEdge(newVertex, vertex(1)));
	}
	return true;
}

bool RegularEdge::refine(std::vector<RegularEdge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<Edge*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge
bool ConstrainedEdge::refine(std::vector<Edge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new ConstrainedEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new ConstrainedEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, vertex(1)));
	}
	return true;
}

bool ConstrainedEdge::refine(std::vector<ConstrainedEdge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<Edge*>&>(vNewEdgesOut),
				newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge
bool ConstrainingEdge::refine(std::vector<Edge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new ConstrainingEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new ConstrainingEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, vertex(1)));
	}
	return true;
}

bool ConstrainingEdge::refine(std::vector<ConstrainingEdge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<Edge*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

template <> size_t
ConstrainingEdge::
num_constrained<Vertex>() const
{
	return num_constrained_vertices();
}

template <> size_t
ConstrainingEdge::
num_constrained<Edge>() const
{
	return num_constrained_edges();
}

template <> Vertex*
ConstrainingEdge::
constrained<Vertex>(size_t ind) const
{
	return constrained_vertex(ind);
}

template <> Edge*
ConstrainingEdge::
constrained<Edge>(size_t ind) const
{
	return constrained_edge(ind);
}

}// end of namespace
