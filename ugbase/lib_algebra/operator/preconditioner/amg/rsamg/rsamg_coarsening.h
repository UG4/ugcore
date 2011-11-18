#ifndef __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
#define __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__

#include "../graph.h"
#include "rsamg_nodeinfo.h"

namespace ug
{

template<typename Matrix_type>
void
CreateStrongConnectionGraph(const Matrix_type &A, cgraph &graph, double theta=0.25);

void CreateMeasureOfImportancePQ(const cgraph &strong, const cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes);

void CreateAggressiveCoarseningGraph(const cgraph &graph, cgraph &graph2, const AMGNodes &nodes,
		int nrOfPaths, int *posInConnections);


void CreateMeasureOfImportanceAggressiveCoarseningPQ(const cgraph &graphAC, nodeinfo_pq_type &PQ, AMGNodes &nodes);

int Coarsen(const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bAggressiveCoarsening);

void PreventFFConnections(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes);


void RemoveUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i);
void MarkUnassignedNeighborsFine(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, bool bMarkAsFineIndirect);
void ChangeRatingOfUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, int change);


} // namespace ug

#endif // __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
