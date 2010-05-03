// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d19

#include <sstream>
#include "lib_grid/lib_grid.h"
#include "tmp_lib_grid_methods.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	TEMPORARY METHODS

////////////////////////////////////////////////////////////////////////
//	SaveSelectedEdgesToObj
bool SaveMarkedEdgesToObj(Grid& grid, const char* filename,
						 SubsetHandler& sh, APosition& aPos)
{
	ofstream out(filename);

	if(out)
	{
	//	write the header
		out << "# exported from libGrid" << endl;

		int indexDimension = 1;

	//	store indices in the vertices
		AInt aInt;
		grid.attach_to_vertices_dv(aInt, -1);
		Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);

	//	write vertex data
		{
			Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
			int counter = 1;
		//	iterate through all edges (from all subsets) and assign indices and output the vertices
			for(int si = 0; si < sh.num_subsets(); ++si)
			{
				for(EdgeBaseIterator iter = sh.begin<EdgeBase>(si);
					iter != sh.end<EdgeBase>(si); ++iter)
				{
					for(int i = 0; i < 2; ++i)
					{
						VertexBase* v = (*iter)->vertex(i);
						if(aaInt[v] == -1)
						{
							out << "v " << aaPos[v].x << " " << aaPos[v].y << " " <<
											aaPos[v].z << endl;
							aaInt[v] = counter++;
						}
					}
				}
			}
		}

	//	write objects
		for(int si = 0; si < sh.num_subsets(); ++si)
		{
			if(sh.subset_info(si).name.size() > 0)
				out << "o " << sh.subset_info(si).name << endl;
			else
				out << "o " << "sub_" << si << endl;

		//	write material reference
			out << "usemtl (null)" << endl;


		//	write the edges of this subset
			for(EdgeBaseIterator iter = sh.begin<EdgeBase>(si);
				iter != sh.end<EdgeBase>(si); ++iter)
			{
				EdgeBase* e = *iter;
				out << "f";
				for(uint i = 0; i < 2; ++i)
				{
					out << " " << aaInt[e->vertex(i)];
					for(int j = 1; j < indexDimension; ++j)
						out << '/' << aaInt[e->vertex(i)];
				}
				out << endl;
			}
		}

	//	clean up
		grid.detach_from_vertices(aInt);
	//	.obj done
		out.close();
		return true;
	}
	return false;
}
/*
template<typename Vector>
number TriangleValueMin(Vector& v1, Vector& v2, Vector& v3)
{
	Vector e1, e2, e3;
	VecSubtract(e1, v2, v1);
	VecSubtract(e2, v3, v1);
	VecNormalize(e1, e1);
	VecNormalize(e2, e2);
	
	number a1 = VecDot(e1, e2);
	
	VecSubtract(e3, v3, v2);
	VecNormalize(e3, e3);
	
	number a2 = VecDot(e3, e2);
	
	VecScale(e3, e3, -1.0f);
	number a3 = VecDot(e3, e1);
	
	return std::min(1.0-fabs(a1), std::min(1.0-fabs(a2), 1.0-fabs(a3)));
}

////////////////////////////////////////////////////////////////////////
template<typename Iterator, typename TVertexPositionAccessor>
number CalculateAreaValue(Iterator trisBegin, Iterator trisEnd, TVertexPositionAccessor& aaPos)
{
	number value = 1.0;
	for(Iterator iter = trisBegin; iter != trisEnd; iter++)
		value = std::min(value, TriangleValueMin(aaPos[(*iter)->vertex(0)], aaPos[(*iter)->vertex(1)], aaPos[(*iter)->vertex(2)]));

	return value;
}

////////////////////////////////////////////////////////////////////////
template<typename TVertexPositionAccessor>
inline number CalculateAreaValue(Grid& grid, VertexBase* v, 
								TVertexPositionAccessor& aaPos)
{
	return CalculateAreaValue(grid.associated_faces_begin(v),
							grid.associated_faces_end(v), aaPos);
}
*/
/*
////////////////////////////////////////////////////////////////////////
template <class TVertexPositionAccessor>
inline number EdgeLengthSq(EdgeBase* e, TVertexPositionAccessor& aaPos)
{
	return VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
}

////////////////////////////////////////////////////////////////////////
template <class TVertexPositionAccessor>
number GetMinEdgeLenSq(EdgeBaseIterator iterBegin, EdgeBaseIterator iterEnd,
					   TVertexPositionAccessor& aaPos)
{
	if(iterBegin == iterEnd)
		return 0;

	number minLen = EdgeLengthSq(*iterBegin, aaPos);
	++iterBegin;

	while(iterBegin != iterEnd){
		number elen = EdgeLengthSq(*iterBegin, aaPos);
		if(elen < minLen)
			minLen = elen;
		++iterBegin;
	}
	return minLen;
}

////////////////////////////////////////////////////////////////////////
template <class TVertexPositionAccessor>
inline number GetMinConnectedEdgeLenSq(Grid& grid, VertexBase* v,
									   TVertexPositionAccessor& aaPos)
{
	return GetMinEdgeLenSq(grid.associated_edges_begin(v),
							grid.associated_edges_end(v), aaPos);
}
*/


}//	end of namespace
