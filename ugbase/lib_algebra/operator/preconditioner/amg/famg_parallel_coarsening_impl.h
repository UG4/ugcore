/**
 * \file famg_parallel_coarsening_impl.h
 *
 * \author Martin Rupp
 *
 * \date 09.02.2011
 *
 * implementation file for famg parallel coarsening
 *
 * Goethe-Center for Scientific Computing 2011.
 *
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_PARALLEL_COARSENING_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_PARALLEL_COARSENING_IMPL_H__

#include "ug.h"
#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"

#include "stopwatch.h"
#include "common/assert.h"


namespace ug
{


template<typename matrix_type, typename prolongation_matrix_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type>::color_process_graph()
{
	stopwatch SW;

	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	UG_LOG("\ncoloring processor OL graph..."); if(bTiming) SW.start();
	// add processors of overlap 1 to pidsOL
	std::set<int> pidsOL;
	for(IndexLayout::iterator iter = OLSendingLayout.begin(); iter != OLSendingLayout.end(); ++iter)
		pidsOL.insert(OLSendingLayout.proc_id(iter));
	for(IndexLayout::iterator iter = OLReceivingLayout.begin(); iter != OLReceivingLayout.end(); ++iter)
		pidsOL.insert(OLReceivingLayout.proc_id(iter));

	m_myColor = ColorProcessorGraph(communicator, pidsOL, processesWithLowerColor, processesWithHigherColor);
	if(bTiming) UG_LOG("took " << SW.ms() << " ms");
	UG_LOG("my color is " << m_myColor);
}

template<typename matrix_type, typename prolongation_matrix_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type>::receive_coarsening_from_processes_with_lower_color()
{
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();

	// issue receive of coarsening data from processes with lower color
	if(processesWithLowerColor.size() > 0)
	{
		stopwatch SW;
		UG_LOG("\nWaiting for processes "); if(bTiming) SW.start();
		stdvector< stdvector<char> > coarseNodes;
		coarseNodes.resize(processesWithLowerColor.size());

		for(size_t i=0; i<processesWithLowerColor.size(); i++)
		{
			int pid = processesWithLowerColor[i];

			size_t s = OLReceivingLayout.interface(pid).size();
			coarseNodes[i].resize(s, -1);
			UG_LOG(pid << ", awaiting  " << s << " bytes.");
			communicator.receive_raw(pid, &coarseNodes[i][0], s);
		}
		UG_LOG("which have higher color to receive coarse nodes... ");
		communicator.communicate();
		UG_LOG("done. processing data...");

		// set nodes coarse
		for(size_t i=0; i<processesWithLowerColor.size(); i++)
		{
			int pid = processesWithLowerColor[i];
			UG_LOG("\nfrom processor " << pid << ": ");

			int j=0;
			IndexLayout::Interface &interface = OLReceivingLayout.interface(pid);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t index = interface.get_element(iter);
				//UG_LOG((int)coarseNodes[i][j] << " ");
				UG_LOG(index);
				if(rating.is_master(index))
					UG_LOG(" master");
				if(coarseNodes[i][j++])
				{
					rating.external_set_coarse(index);
					size_t newIndex = rating.newIndex[index];
					UG_ASSERT(newIndex != -1, "");

					if(rating.is_master(index))
						nextLevelMasterLayout.interface(pid).push_back(newIndex);
					else
						nextLevelSlaveLayout.interface(pid).push_back(newIndex);

					UG_LOG(" coarse ");
				}
				else
				{
					rating.set_fine(index);
					UG_LOG(" fine ");
				}
			}
		}
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");
	}
	else
		UG_LOG("\nno processes with lower color.")
}

template<typename matrix_type, typename prolongation_matrix_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type>::send_coarsening_data_to_processes_with_higher_color()
{
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();

	stopwatch SW;
	UG_LOG("\nsend coarsening data to processes "); if(bTiming) SW.start();
	for(size_t i=0; i<processesWithHigherColor.size(); i++)
	{
		int pid = processesWithHigherColor[i];
		UG_LOG(pid << ": ");
		BinaryStream s;

		IndexLayout::Interface &interface = OLSendingLayout.interface(pid);
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t index = interface.get_element(iter);
			char bCoarse = rating[i].is_coarse();
			UG_LOG(index << (bCoarse ? " is coarse " : " is fine "));
			Serialize(s, bCoarse);
		}

		UG_LOG("sending " << s.size() << " of data to pid " << pid << " ");
		communicator.send_raw(pid, s.buffer(), s.size(), true);
	}
	UG_LOG("with higher color...")

	communicator.communicate();
	UG_LOG("done.");
	if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
}


template<typename matrix_type, typename prolongation_matrix_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type>::create_OL2_matrix()
{
	stopwatch SW;

	// get the Overlap 2 matrix
	UG_LOG("\nGenerate Overlap 2..."); if(bTiming) SW.start();

	std::vector<IndexLayout> masterLayouts, slaveLayouts;
	IndexLayout totalMasterLayout, totalSlaveLayout;
	GenerateOverlap(A, A_OL2, totalMasterLayout, totalSlaveLayout, masterLayouts, slaveLayouts, 2);

	std::vector<MathVector<3> > vec2(&m_famg.amghelper.positions[0], &m_famg.amghelper.positions[m_famg.amghelper.size]);
	vec2.resize(A_OL2.num_rows());

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// get positions of newly created indices
	UG_LOG("\nGet positions of Overlap 2 nodes..."); if(bTiming) SW.start();
	ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	communicator.send_data(totalMasterLayout, copyPol);
	communicator.receive_data(totalSlaveLayout, copyPol);
	communicator.communicate();
	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// write overlap 2 matrix as debug output
	std::stringstream ss; ss << "A_OL2_" << pcl::GetProcRank() << ".mat";
	WriteMatrixToConnectionViewer(ss.str().c_str(), A_OL2, &vec2[0], 2);



	AddLayout(OLSendingLayout, masterLayouts[0]);
	AddLayout(OLSendingLayout, slaveLayouts[0]);

	AddLayout(OLReceivingLayout, slaveLayouts[0]);
	AddLayout(OLReceivingLayout, masterLayouts[0]);


	size_t N = A_OL2.num_rows();
	famg_nodes rating(N);

	// set overlap type of the nodes
	for(IndexLayout::iterator iter = masterLayouts[0].begin(); iter != masterLayouts[0].end(); ++iter)
	{
		IndexLayout::Interface &interface = masterLayouts[0].interface(iter);
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			size_t index = interface.get_element(iter2);
			rating.OLtype[index] = 1;
		}
	}
	for(size_t i=0; i < slaveLayouts.size(); i++)
	{
		for(IndexLayout::iterator iter = slaveLayouts[i].begin(); iter != slaveLayouts[i].end(); ++iter)
		{
			IndexLayout::Interface &interface = slaveLayouts[i].interface(iter);
			for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t index = interface.get_element(iter2);
				rating.OLtype[index] |= 1 << i+1;
			}
		}
	}

	rating.print_OL_types();

	// print some layouts
	UG_LOG("\n");
	for(size_t i=0; i<slaveLayouts.size(); i++)
	{
		UG_LOG("Slave OL Level " << i << ":\n");
		PrintLayout(slaveLayouts[i]);
	}
	for(size_t i=0; i<masterLayouts.size(); i++)
	{
		UG_LOG("Master OL Level " << i << ": ");
		PrintLayout(masterLayouts[i]);
	}
}

}

#endif // __H__LIB_ALGEBRA__AMG__FAMG_PARALLEL_COARSENING_IMPL_H__
