// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 01.06.2011 (m,d,y)

#ifndef __H__UG__rule_util__
#define __H__UG__rule_util__

namespace ug{
namespace shared_rules{

const int MAX_NUM_INDS_OUT = 256;

int RecursiveRefine(int* newIndsOut, int* newEdgeVrts,
					const int faceVrtInds[][4],
					const int faceEdgeInds[][4],
					int numVrts, int numEdges, int numFaces);

}//	end of namespace
}//	end of namespace

#endif
