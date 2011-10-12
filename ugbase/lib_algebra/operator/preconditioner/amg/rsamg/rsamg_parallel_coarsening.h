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


void FullSubdomainBlocking(const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes);

void ColoringCoarsen(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &OLCoarseningSendLayout, IndexLayout &OLCoarseningReceiveLayout,
		const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes);

}
#endif /* RSAMG_PARALLEL_COARSENING_H_ */

#endif /* UG_PARALLEL */
