#ifndef __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
#define __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__

#include "../graph.h"
#include "rsamg_nodeinfo.h"

namespace ug
{

template<typename Matrix_type>
void
CreateStrongConnectionGraph(const Matrix_type &A, cgraph &graph, double theta=0.25);

void CreateMeasureOfImportancePQ(cgraph &strong, cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes);

void CreateAggressiveCoarseningGraph(cgraph &graph, cgraph &graph2, AMGNodes &nodes,
		int nrOfPaths, int *posInConnections);


void CreateMeasureOfImportanceAggressiveCoarseningPQ(cgraph &graphAC, nodeinfo_pq_type &PQ, AMGNodes &nodes);

int Coarsen(cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes);

void PreventFFConnections(cgraph &graphS, cgraph &graphST, AMGNodes &nodes);

} // namespace ug

#endif // __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
