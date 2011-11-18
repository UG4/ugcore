/*
 * rsamg_parallel_coarsening.h
 *
 *  Created on: 12.10.2011
 *      Author: mrupp
 */

#ifdef UG_PARALLEL

#ifndef RSAMG_PARALLEL_COARSENING_H_
#define RSAMG_PARALLEL_COARSENING_H_

#include "../graph.h"
#include "rsamg_nodeinfo.h"

#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"


namespace ug
{


class IParallelCoarsening
{
public:
	virtual ~IParallelCoarsening() {}
	virtual void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
				const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric) = 0;
	virtual int overlap_depth_master() = 0;
	virtual int overlap_depth_slave() = 0;
	virtual const char *tostring() = 0;
};


IParallelCoarsening *GetFullSubdomainBlockingCoarsening();
IParallelCoarsening *GetColorCoarsening();
IParallelCoarsening *GetRS3Coarsening();
IParallelCoarsening *GetCLJPCoarsening();
IParallelCoarsening *GetFalgoutCoarsening();
IParallelCoarsening *GetMinimumSubdomainBlockingCoarsening();
IParallelCoarsening *GetCoarseGridClassificationCoarsening();
IParallelCoarsening *GetSimpleParallelCoarsening();
}
#endif /* RSAMG_PARALLEL_COARSENING_H_ */

#endif /* UG_PARALLEL */
