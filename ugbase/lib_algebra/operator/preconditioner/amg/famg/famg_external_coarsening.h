/*
 * famg_external_coarsening.h
 *
 *  Created on: 21.11.2011
 *      Author: mrupp
 */

#include <fstream>

#include "../rsamg/rsamg.h"
#include "../rsamg/rsamg_coarsening.h"


#ifndef FAMG_EXTERNAL_COARSENING_H_
#define FAMG_EXTERNAL_COARSENING_H_

namespace ug{

template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::external_coarsening_calculate_prolongation()
{
	AMG_PROFILE_FUNC();
	Stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "calculate prolongation... "); if(bTiming) SW.start();
	size_t N = A.num_rows();
	for(size_t i=0; i<N; i++)
	{
		if(rating[i].is_fine_direct() && rating.i_must_assign(i))
			calculator.get_all_neighbors_interpolation(i, PoldIndices, rating);
	}
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
}

template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::rs_amg_external_coarsening()
{
	AMG_PROFILE_FUNC();
	Stopwatch SW;
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



template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::write_coarsening()
{
	PROFILE_FUNC_GROUP("debug");

	stdvector<MathVector<3> > &vec = m_famg.m_amghelper.positions[level];

	std::ios_base::openmode mode;
	if(pcl::GetProcRank()==0)
	{
		CreateDirectory((std::string(m_famg.m_writeMatrixPath) + "/coarsening").c_str(), 0777);
		mode = std::ios::out | ios::binary;
	}
	else mode = std::ios::out | std::ios::app | ios::binary;

	for(int i=0; i<pcl::GetProcRank(); i++)
		pcl::AllProcsTrue(true);
	{
		std::fstream file((std::string(m_famg.m_writeMatrixPath) + "/coarsening/rating" + ToString(level)).c_str(), mode);
		for(size_t i=0; i<vec.size(); i++)
		{
			Serialize(file, vec[i]);
			Serialize(file, rating[i].rating);
		}
	}
	for(int i=pcl::GetProcRank(); i<pcl::GetNumProcesses(); i++)
		pcl::AllProcsTrue(true);
}

template<int Ti>
struct MathVectorComp {
  inline bool operator() (const MathVector<Ti> &a, const MathVector<Ti> &b) const
  {
	  for(size_t i=0; i<Ti; i++)
		  if(a[i] != b[i]) return a[i] < b[i];
	  return false;
  }
};

template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::read_coarsening()
{
	PROFILE_FUNC_GROUP("debug");

	stdvector<MathVector<3> > &vec = m_famg.m_amghelper.positions[level];
	typedef std::map<MathVector<3>, size_t, MathVectorComp<3> > MathVecMap;
	MathVecMap m;
	for(size_t i=0; i<vec.size(); i++)
		m[vec[i]] = i;

	std::fstream file((std::string(m_famg.m_writeMatrixPath) + "/coarsening/rating" + ToString(level)).c_str(), std::ios::in | std::ios::binary);
	for(size_t i=0; i<vec.size(); i++)
	{
		MathVector<3> pos;
		double r;
		Deserialize(file, pos);
		Deserialize(file, r);

		MathVecMap::iterator it = m.find(pos);
		if(it != m.end())
		{
			rating[it->second].rating = r;
			UG_LOG("Setting node " << it->second << " to " << rating[it->second]);
		}
	}
}
}

#endif /* FAMG_EXTERNAL_COARSENING_H_ */
