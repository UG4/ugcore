// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)
 
#include "grid_objects_1d.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	Edge
bool Edge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex,
				VertexBase** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new Edge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new Edge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new Edge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new Edge(newVertex, vertex(1)));
	}
	return true;
}

bool Edge::refine(std::vector<Edge*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge
bool ConstrainedEdge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
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

bool ConstrainedEdge::refine(std::vector<ConstrainedEdge*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut),
				newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge
bool ConstrainingEdge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
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

bool ConstrainingEdge::refine(std::vector<ConstrainingEdge*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

template <> size_t
ConstrainingEdge::
num_constrained<VertexBase>() const
{
	return num_constrained_vertices();
}

template <> size_t
ConstrainingEdge::
num_constrained<EdgeBase>() const
{
	return num_constrained_edges();
}

template <> VertexBase*
ConstrainingEdge::
constrained<VertexBase>(size_t ind) const
{
	return constrained_vertex(ind);
}

template <> EdgeBase*
ConstrainingEdge::
constrained<EdgeBase>(size_t ind) const
{
	return constrained_edge(ind);
}

}// end of namespace
