// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)
 
#include "geometric_objects_1d.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	Edge::refine
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
//	ConstrainedEdge::refine
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
//	ConstrainingEdge::refine
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

}// end of namespace
