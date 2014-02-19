//	created by Martin Stepniewski
//	mastep@gmx.de
//	y08 m12 d07

#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "subdivision_loop.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	ProjectToLimitLoop
/// projects surface vertices to their limit subdivision surface position
bool ProjectToLimitLoop(Grid& grid, APosition& aProjPos)
{
//	grid management
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid,aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaProjPos(grid, aProjPos);

//	needed variables
	double const pi = 3.14159265;
	double x = 0;
	double y = 0;
	double z = 0;
	int valence = 0;
	const int numPrecalculated = 10;
	double beta[numPrecalculated];
	double b = 0;
	double chi = 0;

//	calculate weights for subdivision mask
	for(int i = 1; i < numPrecalculated; ++i)
	{
		double tmp = 0.375 + 0.25 * cos( (2.0 * pi) / (float)i );
		beta[i] = ( 0.625 - tmp * tmp ) / (float)i ;
	}

	beta[0] = 0;
	beta[6] = 0.0625;

//	iterate through all vertices, evaluate their limit positions and save them in their projection attachment
	for(VertexIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
	{
		Vertex* v = *vIter;
		valence = 0;
		x = 0;
		y = 0;
		z = 0;

		for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(v); eIter != grid.associated_edges_end(v); ++eIter)
		{
			Edge* e = *eIter;
			valence++;

			if(valence >= numPrecalculated)
			{
				double tmp = 0.375 + 0.25 * cos( (2.0*pi) / (float)valence );
				b = (0.625 - tmp*tmp) / (float)valence;
			}

			else
				b = beta[valence];

			chi = 1.0 / (0.375 / b + valence);

			if(aaPos[v].x() == aaPos[e->vertex(0)].x() && aaPos[v].y() == aaPos[e->vertex(0)].y() && aaPos[v].z() == aaPos[e->vertex(0)].z())
			{
				x += aaPos[e->vertex(1)].x();
				y += aaPos[e->vertex(1)].y();
				z += aaPos[e->vertex(1)].z();
			}

			else
			{
				x += aaPos[e->vertex(0)].x();
				y += aaPos[e->vertex(0)].y();
				z += aaPos[e->vertex(0)].z();
			}
		}

		x*=chi;
		y*=chi;
		z*=chi;

		aaProjPos[v].x() = aaPos[v].x() * (1.0 - (float)valence * chi);
		aaProjPos[v].y() = aaPos[v].y() * (1.0 - (float)valence * chi);
		aaProjPos[v].z() = aaPos[v].z() * (1.0 - (float)valence * chi);

		aaProjPos[v].x() += x;
		aaProjPos[v].y() += y;
		aaProjPos[v].z() += z;
	}

	return true;
}

}//	end of namespace

