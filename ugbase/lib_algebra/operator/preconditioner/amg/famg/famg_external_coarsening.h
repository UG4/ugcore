/*
 * famg_external_coarsening.h
 *
 *  Created on: 21.11.2011
 *      Author: mrupp
 */

#include "../rsamg/rsamg.h"
#include "../rsamg/rsamg_coarsening.h"


#ifndef FAMG_EXTERNAL_COARSENING_H_
#define FAMG_EXTERNAL_COARSENING_H_

namespace ug{

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::external_coarsening_calculate_prolongation()
{
	AMG_PROFILE_FUNC();
	stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "calculate prolongation... "); if(bTiming) SW.start();
	size_t N = A.num_rows();
	for(size_t i=0; i<N; i++)
	{
		if(rating[i].is_fine_direct() && rating.i_must_assign(i))
			calculator.get_all_neighbors_interpolation(i, PoldIndices, rating);
	}
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
}

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::rs_amg_external_coarsening()
{
	AMG_PROFILE_FUNC();
	stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "external coarsening... " ); if(bTiming) SW.start();
	nodeinfo_pq_type PQ;
	size_t N = A_OL2.num_rows();

	// use RSAMG's methods to do standard rs coarsening
#ifdef UG_PARALLEL
	AMGNodes nodes(N, PN);
	UG_DLOG(LIB_ALG_AMG, 0, m_famg.m_pParallelCoarsening->tostring() << "\n");
#else
	AMGNodes nodes(N);
#endif
	cgraph graphS, graphST;
	CreateStrongConnectionGraph(A_OL2, graphS, m_famg.m_dStrongConnectionExternal, nodes);

	graphST.set_as_transpose_of(graphS);
	CreateMeasureOfImportancePQ(graphS, graphST, PQ, nodes);
#ifdef UG_PARALLEL
	m_famg.m_pParallelCoarsening->coarsen(PN, masterLayouts, slaveLayouts, graphS, graphST, PQ, nodes, true);
#else
	Coarsen(graphS, graphST, PQ, nodes, true, false);
	PreventFFConnections(graphS, graphST, nodes);
#endif

	// aggressive coarsening
	if(0 && m_famg.get_aggressive_coarsening() == true && level == 0)
	{
		UG_DLOG(LIB_ALG_AMG, 2, std::endl << "building graph2... ");
		UG_DLOG(LIB_ALG_AMG, 2, "unassigned = " << nodes.get_unassigned() << "\n");

		//unassigned = 0; ??
		cgraph graphAC(N);
		size_t m_iAggressiveCoarseningNrOfPaths = 2;
		stdvector<int> posInConnections; posInConnections.resize(N, -1);
		CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, m_iAggressiveCoarseningNrOfPaths, &posInConnections[0]);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, nodes);
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
		// coarsen 2
		//------------------

		if(nodes.get_unassigned() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 2, std::endl << "skipping coarsening2: no unassigned nodes.");
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 2, std::endl << "coarsening2... ");
			Coarsen(graphAC, graphAC, PQ, nodes, true, true);
			//PreventFFConnections(graphS, graphST, nodes);
		}
	}

	// coarsening done. transfer AMG nodeinfo -> FAMG nodeinfo
	for(size_t i=0; i<N; i++)
	{
		if(graphS.is_isolated(i))
			rating.set_fine(i);
		else if(nodes[i].is_coarse())
			rating.external_set_coarse(i);
		else if(nodes[i].is_fine_direct() || nodes[i].is_fine())
			rating.set_fine(i);
		else if(nodes[i].is_dirichlet())
			rating.set_dirichlet(i);
		else
		{
			//UG_LOG("nodes " << i << " = " << nodes[i] << "\n");
			rating.set_fine(i);
			//UG_LOG("rating " << i << " = " << rating[i] << "\n");
		}
	}

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
	//if(nodes[i].is_unassigned_fine_indirect())
	//	rating.set_aggressive_fine(i);
}

}

#endif /* FAMG_EXTERNAL_COARSENING_H_ */
