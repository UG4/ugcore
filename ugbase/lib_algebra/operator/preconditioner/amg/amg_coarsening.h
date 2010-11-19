#ifndef __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
#define __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__

#include "graph.h"
#include "amg_nodeinfo.h"

namespace ug
{

template<typename Matrix_type>
void
CreateStrongConnectionGraph(const Matrix_type &A, cgraph &graph, double theta=0.25);

void CreateMeasureOfImportancePQ(cgraph &strong, cgraph &strongT, nodeinfo_pq_type &PQ, int &unassigned, amg_nodeinfo *nodes);

void CreateAggressiveCoarseningGraph(cgraph &graph, cgraph &graph2, amg_nodeinfo *nodes,
		int nrOfPaths, int *posInConnections);


void CreateMeasureOfImportanceAggressiveCoarseningPQ(cgraph &graphAC, nodeinfo_pq_type &PQ, int &unassigned, int &iNrOfCoarse, int *newIndex, amg_nodeinfo *nodes);

int Coarsen(cgraph &graph, nodeinfo_pq_type &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, amg_nodeinfo *nodes);

void PreventFFConnections(cgraph &graphS, cgraph &graphST, amg_nodeinfo *nodes, int *newIndex, int &nrOfCoarse);

} // namespace ug

#endif // __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
