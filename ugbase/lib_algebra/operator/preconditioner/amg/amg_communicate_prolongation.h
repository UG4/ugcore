/*
 * amg_communicate_prolongation.h
 *
 *  Created on: 27.04.2011
 *      Author: mrupp
 */

#ifndef AMG_COMMUNICATE_PROLONGATION_H_
#define AMG_COMMUNICATE_PROLONGATION_H_
#include "pcl/pcl.h"
#include "row_sender.h"

namespace ug
{

/**
 *
 *
 * send prolongation from masters to slaves
 *
 */

// FAMGLevelCalculator::communicate_prolongation
//---------------------------------------------------------------------------
/** send prolongation from Master nodes to Slave nodes.
 */
template<typename matrix_type, typename TNodes>
void communicate_prolongation(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &masterLayout, IndexLayout &slaveLayout,
		IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout,
		ParallelNodes &PN, matrix_type &P, TNodes &nodes)
{
	//pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();

	UG_DLOG(LIB_ALG_AMG, 4, "\n********** communicate_prolongation **************\n\n")
	AMG_PROFILE_FUNC();

	//UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelCommunicateProlongation);

	//UG_DLOG(LIB_ALG_AMG, 1, "\n\n\nOL2Layout:\n")
	//PrintLayout(A_OL2.get_communicator(), OL2MasterLayout, OL2SlaveLayout);

	// 1. get global IDs (we would only need IDs on OL1)
	/*std::vector<AlgebraID> localToGlobal;
	GenerateGlobalAlgebraIDs(localToGlobal, A_OL2.num_rows(), OL2MasterLayout, OL2SlaveLayout);
	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		for(size_t i=0; i<localToGlobal.size(); i++)
		{	UG_DLOG(LIB_ALG_AMG, 4, "local " << i << " = " << localToGlobal[i] << "\n");	}
	}
	NewNodesNummerator globalToLocal(localToGlobal);*/

	// 2. send from master interface of A to slaves rows of P (with global IDs)

	UG_DLOG(LIB_ALG_AMG, 4, "send prolongation\n");

	rowSender CRowSender(P, PN);
	rowSender.set_mode(RSM_OVERWRITE_MODE);
	rowSender.set_exchange_over_interfaces(masterLayout, slaveLayout);
	rowSender.set_original_interfaces(masterLayout, slaveLayout);
	rowSender.add_new_interfaces_to(newMasterLayout, newSlaveLayout);
	rowSender.communicate(communicator);


	PN.resize(P.num_rows());
	for(size_t i=0; i<P.num_rows(); i++)
	{
		if(P.num_connections(localRowIndex) == 1)
		{
			if(nodes[i].is_coarse() == false)
			{
				nodes.external_set_coarse(j);
				UG_DLOG(LIB_ALG_AMG, 4, " post-setted " << i << " coarse.\n");
			}
		}
		else
		{
			if(nodes[i].is_fine() == false)
			{
				UG_ASSERT(nodes.is_slave(i, 0), localRowIndex);
				nodes.set_fine(localRowIndex);
				UG_DLOG(LIB_ALG_AMG, 4, " post-setted " << localRowIndex << " fine.\n");
			}

			for(typename matrix_type::const_row_iterator conn = P.begin_row(i); conn != P.end_row(i); ++conn)
			{
				size_t j = conn.index();
				if(nodes[j].is_coarse() == false)
				{
					nodes.external_set_coarse(j);
					UG_DLOG(LIB_ALG_AMG, 4, " post-setted " << j << " coarse.\n");
				}
			}
		}
	}

}

template<typename T>
void AddCoarseToLayout(IndexLayout &src, IndexLayout &dest, const T &nodes)
{
	for(IndexLayout::iterator iter = src.begin(); iter != src.end(); ++iter)
	{
		IndexLayout::Interface &interface = src.interface(iter);
		int pid = src.proc_id(iter);

		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t localIndex = interface.get_element(iter);
			if(nodes[localIndex].is_coarse())
			{
				IndexLayout::Interface &destInterface = dest.interface(pid);
				destInterface.push_back(localIndex);
			}
		}
	}
}



}
#endif /* FAMG_COMMUNICATE_PROLONGATION_H_ */
