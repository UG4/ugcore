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

#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"

#include "stopwatch.h"


namespace ug
{


// FAMGLevelCalculator::color_process_graph
//-------------------------------------------
/** calculates coloring of processes for coarsening
 * calculates m_myColor so that each processor which is connected through
 * OLCoarseningSendLayout or OLCoarseningReceiveLayout gets a different color,
 * and gets arrays processesWithLowerColor and processesWithHigherColor
 * \sa ColorProcessorGraph
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	color_process_graph()
{
	stopwatch SW;

	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	UG_LOG("\ncoloring processor OL graph..."); if(bTiming) SW.start();
	// add processors of overlap 1 to pidsOL
	std::set<int> pidsOL;
	for(IndexLayout::iterator iter = OLCoarseningSendLayout.begin(); iter != OLCoarseningSendLayout.end(); ++iter)
		pidsOL.insert(OLCoarseningSendLayout.proc_id(iter));
	for(IndexLayout::iterator iter = OLCoarseningReceiveLayout.begin(); iter != OLCoarseningReceiveLayout.end(); ++iter)
		pidsOL.insert(OLCoarseningReceiveLayout.proc_id(iter));

	m_myColor = ColorProcessorGraph(communicator, pidsOL, processesWithLowerColor, processesWithHigherColor);
	if(bTiming) UG_LOG("took " << SW.ms() << " ms");
	UG_LOG("my color is " << m_myColor);
}

// FAMGLevelCalculator::receive_coarsening_from_processes_with_lower_color
//---------------------------------------------------------------------------
/** receives coarsening data from processes with lower color
 * we receive data from processes in processesWithLowerColor
 * they send us data about which nodes they have set coarse
 * for this, we use OLCoarseningReceiveLayout, so that nodes we have in common with processor pid1,
 * are in the interface OLCoarseningReceiveLayout.interface(pid1)
 *
 * If i is in the OLCoarseningReceiveLayout to processor pid1,
 * and it is coarse, we have to add it to an interface for the next level.
 * that is: nextMasterLayout.interface(pid1) if i is master on this processor,
 * or nextSlaveLayout.interface(pid1) if i is slave on this processor
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	receive_coarsening_from_processes_with_lower_color()
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

			size_t s = OLCoarseningReceiveLayout.interface(pid).size();
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
			UG_LOG("\nfrom processor " << pid << ":\n");

			int j=0;
			IndexLayout::Interface &interface = OLCoarseningReceiveLayout.interface(pid);
			IndexLayout::Interface &nextMasterInterface = nextLevelMasterLayout.interface(pid);
			IndexLayout::Interface &nextSlaveInterface = nextLevelSlaveLayout.interface(pid);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t index = interface.get_element(iter);
				//UG_LOG((int)coarseNodes[i][j] << " ");
				UG_LOG(" " << index);
				if(rating.is_master(index))

					UG_LOG(" master");
				if(coarseNodes[i][j++])
				{
					rating.external_set_coarse(index);
					int newIndex = rating.newIndex[index];
					UG_ASSERT(newIndex != -1, "");

					// add node to next Level Interface
					if(rating.is_master(index))
						nextMasterInterface.push_back(newIndex);
					else
						nextSlaveInterface.push_back(newIndex);

					UG_LOG(" coarse ");
				}
				else
				{
					rating[index].set_uninterpolateable();
					UG_LOG(" fine ");
				}
				UG_LOG("\n");
			}
		}
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");
	}
	else
		UG_LOG("\nno processes with lower color.")

	AH.set_communicator(communicator);
	AH.set_slave_layout(nextLevelSlaveLayout);
	AH.set_master_layout(nextLevelMasterLayout);
}

// FAMGLevelCalculator::send_coarsening_data_to_processes_with_higher_color
//---------------------------------------------------------------------------
/** sends coarsening data to processes with higher color
 * we send data to processes in processesWithHigherColor about which nodes we have set coarse.
 * for this, we use OLCoarseningSendLayout.
 * If i is in the OLCoarseningSendLayout to processor pid1,
 * and it is coarse, we have to add it to an interface for the next level.
 * that is: nextMasterLayout.interface(pid1) if i is master on this processor,
 * or nextSlaveLayout.interface(pid1) if i is slave on this processor
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	send_coarsening_data_to_processes_with_higher_color()
{
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();

	stopwatch SW;
	UG_LOG("\nsend coarsening data to processes "); if(bTiming) SW.start();
	for(size_t i=0; i<processesWithHigherColor.size(); i++)
	{
		int pid = processesWithHigherColor[i];
		UG_LOG(pid << ": ");
		BinaryStream s;

		IndexLayout::Interface &interface = OLCoarseningSendLayout.interface(pid);
		IndexLayout::Interface &nextMasterInterface = nextLevelMasterLayout.interface(pid);
		IndexLayout::Interface &nextSlaveInterface = nextLevelSlaveLayout.interface(pid);

		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t index = interface.get_element(iter);
			char bCoarse = rating[index].is_coarse();
			UG_LOG(index << (bCoarse ? " is coarse " : " is fine "));
			Serialize(s, bCoarse);

			if(bCoarse)
			{
				int newIndex = rating.newIndex[index];
				UG_ASSERT(newIndex != -1, "");

				// add node to next Level Interface
				if(rating.is_master(index))
					nextMasterInterface.push_back(newIndex);
				else
					nextSlaveInterface.push_back(newIndex);
			}
		}

		UG_LOG("sending " << s.size() << " of data to pid " << pid << " ");
		communicator.send_raw(pid, s.buffer(), s.size(), true);
	}
	UG_LOG("with higher color...")

	communicator.communicate();
	UG_LOG("done.");
	if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
}


// FAMGLevelCalculator::add_connections_between_slave_nodes
//---------------------------------------------------------------------------
/** adds connections between slave nodes to OLCoarseningSendLayout/OLCoarseningReceiveLayout
 * in the pcl, there are only connections between master to slave, and not
 * between all slaves. we add those connections, so that we can send coarsening
 * data from one node to all processes which have this node (as master or as slave).
  *
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	add_connections_between_slave_nodes(IndexLayout &masterLayout, IndexLayout slaveLayout)
{
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	UG_LOG("\n\nadd_connections_between_slave_nodes...\n")

	// 1. get for every master node where it is slave
	std::map<size_t, std::vector<int> > slaveOnProc;
	//slaveOnProc.resize(A_OL2.num_rows());

	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		IndexLayout::Interface &interface = masterLayout.interface(iter);
		int pid = masterLayout.proc_id(iter);
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			size_t index = interface.get_element(iter2);
			slaveOnProc[index].push_back(pid);
		}
	}

	// 2. send information to slaves
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		BinaryStream stream;

		IndexLayout::Interface &interface = masterLayout.interface(iter);
		int pid = masterLayout.proc_id(iter);
		size_t j=0;
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2, ++j)
		{
			size_t index = interface.get_element(iter2);
			size_t s = slaveOnProc[index].size()-1;
			if(s <= 0) continue;
			Serialize(stream, j);
			Serialize(stream, s);
			for(size_t i=0; i < s+1; i++)
			{
				int pid2 = slaveOnProc[index][i];
				if(pid2 != pid)
					Serialize(stream, pid2);
			}
		}
		UG_LOG("Sending " << stream.size() << " bytes of data to processor " << pid);
		communicator.send_raw(pid, stream.buffer(), stream.size(), false);
	}

	// 3. communicate
	StreamPack pack;
	std::vector<int> pids;

	for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		int pid = slaveLayout.proc_id(iter);
		pids.push_back(pid);
		communicator.receive_raw(pid, *pack.get_stream(pid));
	}
	communicator.communicate();

	/* sort pids. this is important!
	* imagine 4 processes. when 2 of them have 2 common slave nodes, they have to be inserted in the same order to the interfaces:
	* process 1: interface to 2: 5 (master). interface to 3: 5 (master)
	* process 2 interface to 1: 5 (slave). interface to 4: 10 (slave)
	* process 3 interface to 1: 5 (slave). interface to 4: 10 (slave)
	* process 4 interface to 2: 10 (master). interface to 3: 10 (master)
	* process 1 sends to process 2 / 3: node 0 of our interface is connected with processor 3 / 2
	* process 4 sends to process 2 / 3: node 0 of our interface is connected with processor 3 / 2
	* so process 2 starts with information from processor 1: inserting node 5 to his interface to processor 3. then information from processor 4: node 10 to pid 3.
	* so his interface is 5, 10 with pid 3.
	* if process 3 does in same order, he gets also interface 5, 10 with pid 2. */
	sort(pids.begin(), pids.end());

	// 4. process data

	for(size_t i=0; i<pids.size(); i++)
	{
		int pid = pids[i];
		BinaryStream &stream = *pack.get_stream(pid);
		UG_LOG("Received " << stream.size() << " bytes of data from processor " << pid << ":\n");

		std::vector<size_t> indices;
		IndexLayout::Interface &interface = slaveLayout.interface(pid);
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			indices.push_back(interface.get_element(iter2));

		while(stream.can_read_more())
		{
			size_t index, s;
			Deserialize(stream, index);
			index = indices[index];
			Deserialize(stream, s);
			UG_LOG(" got " << s << " other slave connections from node " << index << ": ")
			for(size_t i=0; i<s; i++)
			{
				int pid2;
				UG_ASSERT(stream.can_read_more(), "stream said to have " << s << " entries, but got only " << i);
				Deserialize(stream, pid2);
				UG_ASSERT(pid2 != pcl::GetProcRank(), "");
				OLCoarseningReceiveLayout.interface(pid2).push_back(index);
				OLCoarseningSendLayout.interface(pid2).push_back(index);
				UG_LOG(pid2 << " ");
			}
			UG_LOG("\n");
		}

	}

	UG_LOG("\n");
}


// FAMGLevelCalculator::create_OL2_matrix
//---------------------------------------------------------------------------
/**
 * creates overlap 2 matrix A_OL2,
 * \sa GnerateOverlap
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::create_OL2_matrix()
{
	stopwatch SW;

	// 1. get the Overlap 2 matrix
	//-------------------------------
	UG_LOG("\nGenerate Overlap 2..."); if(bTiming) SW.start();

	/*{
		UG_SET_DEBUG_LEVELS(4);
		std::vector<IndexLayout> masterLayouts, slaveLayouts;
		IndexLayout totalMasterLayout, totalSlaveLayout;
		matrix_type mat;
		GenerateOverlap2(A, mat, totalMasterLayout, totalSlaveLayout, masterLayouts, slaveLayouts, 1, 0, true, true);
		std::vector<MathVector<3> > vec2(&m_famg.m_amghelper.positions[0], &m_famg.m_amghelper.positions[m_famg.m_amghelper.size]);
		vec2.resize(A_OL2.num_rows());
		// 2. get positions of newly created indices
			//--------------------------------------------
		ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
		pcl::ParallelCommunicator<IndexLayout> &communicator = mat.get_communicator();
		communicator.send_data(totalMasterLayout, copyPol);
		communicator.receive_data(totalSlaveLayout, copyPol);
		communicator.communicate();

		// debug: write overlap 2 matrix as debug output
		WriteMatrixToConnectionViewer(GetProcFilename("A_fullRow", ".mat").c_str(), mat, &vec2[0], 2);
		UG_SET_DEBUG_LEVELS(0);
	}*/


	std::vector<IndexLayout> masterLayouts, slaveLayouts;
	IndexLayout totalMasterLayout, totalSlaveLayout;
	GenerateOverlap(A, A_OL2, totalMasterLayout, totalSlaveLayout, masterLayouts, slaveLayouts, 2);

	std::vector<MathVector<3> > vec2(&m_famg.m_amghelper.positions[0], &m_famg.m_amghelper.positions[m_famg.m_amghelper.size]);
	vec2.resize(A_OL2.num_rows());

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// 2. get data on newly created indices
	//--------------------------------------------
	// use ONE communicate

	UG_LOG("\nGet data on Overlap 2 nodes..."); if(bTiming) SW.start();

	ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	communicator.send_data(totalMasterLayout, copyPol);
	communicator.receive_data(totalSlaveLayout, copyPol);

	std::vector<ComPol_VecCopy< Vector<double> > > vecCopyPol;
	vecCopyPol.resize(m_testvectors.size());
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		vecCopyPol[i].set_vector(&m_testvectors[i]);
		communicator.send_data(totalMasterLayout, vecCopyPol[i]);
		communicator.receive_data(totalSlaveLayout, vecCopyPol[i]);
	}

	communicator.communicate();

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// debug: write overlap 2 matrix as debug output
	WriteMatrixToConnectionViewer(GetProcFilename("A_OL2", ".mat").c_str(), A_OL2, &vec2[0], 2);

	// create OLCoarseningSendLayout, OLCoarseningReceiveLayout
	//------------------------------------------------------------
	AddLayout(OLCoarseningSendLayout, masterLayouts[0]);
	AddLayout(OLCoarseningSendLayout, slaveLayouts[0]);

	AddLayout(OLCoarseningReceiveLayout, slaveLayouts[0]);
	AddLayout(OLCoarseningReceiveLayout, masterLayouts[0]);

	// add connections between slave nodes
	add_connections_between_slave_nodes(masterLayouts[0], slaveLayouts[0]);

	size_t N = A_OL2.num_rows();
	rating.create(N);

	// set overlap type of the nodes
	//-------------------------------
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

	//rating.print_OL_types();


	// debug: print some layouts
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
