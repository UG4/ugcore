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
#include "lib_algebra/parallelization/parallel_coloring.h"
#include "stopwatch.h"


namespace ug
{
void AddConnectionsBetweenSlaves(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &masterLayout, IndexLayout &slaveLayout, IndexLayout &allToAllSend,
		IndexLayout &allToAllReceive);

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
	AMG_PROFILE_FUNC();
	stopwatch SW;

	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	UG_DLOG(LIB_ALG_AMG, 1, "\ncoloring processor OL graph..."); if(bTiming) SW.start();
	// add processors of overlap 1 to pidsOL
	std::set<int> pidsOL;
	for(IndexLayout::iterator iter = OLCoarseningSendLayout.begin(); iter != OLCoarseningSendLayout.end(); ++iter)
	{
		if(OLCoarseningSendLayout.interface(iter).size())
			pidsOL.insert(OLCoarseningSendLayout.proc_id(iter));
	}
	for(IndexLayout::iterator iter = OLCoarseningReceiveLayout.begin(); iter != OLCoarseningReceiveLayout.end(); ++iter)
	{
		if(OLCoarseningReceiveLayout.interface(iter).size())
			pidsOL.insert(OLCoarseningReceiveLayout.proc_id(iter));
	}

	m_myColor = ColorProcessorGraph(communicator, pidsOL, processesWithLowerColor, processesWithHigherColor);
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
	UG_DLOG(LIB_ALG_AMG, 1, "\nmy color is " << m_myColor);
	UG_DLOG(LIB_ALG_AMG, 1, "i am connected to");
	for(std::set<int>::iterator it = pidsOL.begin(); it != pidsOL.end(); ++it)
		UG_DLOG(LIB_ALG_AMG, 1, (*it) << " ");
	UG_DLOG(LIB_ALG_AMG, 1, "\nprocesses with lower color:");
	for(size_t i=0; i<processesWithLowerColor.size(); i++)
		UG_DLOG(LIB_ALG_AMG, 1, processesWithLowerColor[i] << " ");
	UG_DLOG(LIB_ALG_AMG, 1, "\nprocesses with higher color:");
	for(size_t i=0; i<processesWithHigherColor.size(); i++)
		UG_DLOG(LIB_ALG_AMG, 1, processesWithHigherColor[i] << " ");
	UG_DLOG(LIB_ALG_AMG, 1, "\n");
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
	AMG_PROFILE_FUNC();
	UG_DLOG(LIB_ALG_AMG, 1, "\n*** receive coarsening data from processes with lower color ***\n");
	if(processesWithLowerColor.size() == 0)
	{
		UG_DLOG(LIB_ALG_AMG, 1, "no processes with lower color.");
		return;
	}
	pcl::ParallelCommunicator<IndexLayout> communicator = A_OL2.get_communicator();


	// issue receive of coarsening data from processes with lower color
	/*IF_DEBUG(LIB_ALG_AMG, 11)
	{
		pcl::ProcessCommunicator lowerPC = A_OL2.get_process_communicator().create_sub_communicator(processesWithLowerColor);
		communicator.enable_communication_debugging(lowerPC);
	}*/

	stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 1, "\nWaiting for processes "); if(bTiming) SW.start();
	stdvector< stdvector<char> > states;
	states.resize(processesWithLowerColor.size());

	for(size_t i=0; i<processesWithLowerColor.size(); i++)
	{
		int pid = processesWithLowerColor[i];

		size_t s = OLCoarseningReceiveLayout.interface(pid).size();
		states[i].resize(s, -1);
		UG_DLOG(LIB_ALG_AMG, 1, pid << ", awaiting  " << s << " bytes.");
		communicator.receive_raw(pid, &states[i][0], s);
	}
	UG_DLOG(LIB_ALG_AMG, 1, "which have lower color to receive coarse nodes...\n");
	AMG_PROFILE_BEGIN(FAMG_recv_coarsening_communicate);
	communicator.communicate();
	AMG_PROFILE_END();
	UG_DLOG(LIB_ALG_AMG, 1, "done. processing data...");

	// set nodes coarse
	for(size_t i=0; i<processesWithLowerColor.size(); i++)
	{
		int pid = processesWithLowerColor[i];
		UG_DLOG(LIB_ALG_AMG, 3, "\nfrom processor " << pid << ":\n");

		int j=0;
		IndexLayout::Interface &interface = OLCoarseningReceiveLayout.interface(pid);
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t index = interface.get_element(iter);

			int state = states[i][j++];

			if(state == FAMG_FINE_RATING || state == FAMG_UNCALCULATED_FINE_RATING || state == FAMG_AGGRESSIVE_FINE_RATING)
				rating.external_set_uncalculated_fine(index);
			else if(state == FAMG_COARSE_RATING)
				rating.external_set_coarse(index);
			else if(state == 0.0)
				; // nothing
			else
			{
				UG_ASSERT(0, "state is " << state << "?");
			}

			UG_DLOG(LIB_ALG_AMG, 3, index << ": got state " << state << " now " << rating[index] << ", " << rating.OL_type(index) << "\n");
		}
	}
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 3, "took " << SW.ms() << " ms");
}

// FAMGLevelCalculator::send_coarsening_data_to_processes_with_higher_color
//---------------------------------------------------------------------------
/** sends coarsening data to processes with higher color
 * we send data to processes in processesWithHigherColor about which nodes we have set coarse.
 * for this, we use OLCoarseningSendLayout.
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	send_coarsening_data_to_processes_with_higher_color()
{
	AMG_PROFILE_FUNC();
	stopwatch SW; if(bTiming) SW.start();
	UG_DLOG(LIB_ALG_AMG, 1, "\n*** send coarsening data to processes ***\n");
	if(processesWithHigherColor.size() == 0)
	{
		UG_DLOG(LIB_ALG_AMG, 1, "no processes with higher color.");
		return;
	}

	pcl::ParallelCommunicator<IndexLayout> communicator = A_OL2.get_communicator();

	/*IF_DEBUG(LIB_ALG_AMG, 11)
	{
		pcl::ProcessCommunicator higherPC = A_OL2.get_process_communicator().create_sub_communicator(processesWithHigherColor);
		communicator.enable_communication_debugging(higherPC);
	}*/

	for(size_t i=0; i<processesWithHigherColor.size(); i++)
	{
		int pid = processesWithHigherColor[i];
		UG_DLOG(LIB_ALG_AMG, 1, "Process " << pid << ":\n");
		BinaryStream s;

		IndexLayout::Interface &interface = OLCoarseningSendLayout.interface(pid);
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t index = interface.get_element(iter);
			char state = rating[index].get_state();
			UG_ASSERT(state <= 0 , "node " << index << " still has a valid rating");

			Serialize(s, state);

			UG_DLOG(LIB_ALG_AMG, 3, rating.info(index));
		}

		UG_DLOG(LIB_ALG_AMG, 1, "sending " << s.size() << " of data to pid " << pid << "\n");
		communicator.send_raw(pid, s.buffer(), s.size(), true);
	}
	UG_DLOG(LIB_ALG_AMG, 1, "with higher color...")

	AMG_PROFILE_BEGIN(FAMG__coarsening_communicate);
	communicator.communicate();
	AMG_PROFILE_END();
	UG_DLOG(LIB_ALG_AMG, 1, "done.");
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
}



// FAMGLevelCalculator::create_interface
//---------------------------------------------------------------------------
/**
 * this function copies the layout 'layout' to 'nextLevelLayout',
 * but skips fine nodes (which have newIndex[i] == -1), and adds the coarse nodes
 * with the new indices.
 * called once with slave and with master, we have our next level layouts.
 * \param layout
 * \param nextLevelLayout
 * \param newIndex
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	update_interface_with_newIndex(IndexLayout &, IndexLayout &nextLevelLayout, stdvector<int> &newIndex)
{
	AMG_PROFILE_FUNC();
	UG_DLOG(LIB_ALG_AMG, 3, "*** create_interface ***");
	//PrintLayout(layout);
	UG_DLOG(LIB_ALG_AMG, 3, "\n");
	for(IndexLayout::iterator iter = nextLevelLayout.begin(); iter != nextLevelLayout.end(); ++iter)
	{
		IndexLayout::Interface &interface = nextLevelLayout.interface(iter);
		UG_DLOG(LIB_ALG_AMG, 3, "to processor " << nextLevelLayout.proc_id(iter) << ": ");
		/*IndexLayout::Interface &nextLevelinterface = nextLevelLayout.interface(pid);*/
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			size_t &i = interface.get_element(iter2);
			UG_DLOG(LIB_ALG_AMG, 3, i << "(new " << newIndex[i] << ")  ");
			i = newIndex[i];
			/*if(newIndex[i] != -1)
			{
				UG_DLOG(LIB_ALG_AMG, 1, i << "(new " << newIndex[i] << ")  ");
				nextLevelinterface.push_back(newIndex[i]);
			}*/
		}
		UG_DLOG(LIB_ALG_AMG, 3, "\n");
	}

	//UG_DLOG(LIB_ALG_AMG, 1, "total layout:\n");
	//PrintLayout(nextLevelLayout);
}


/*
template<typename matrix_type>
bool GenerateOverlap(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout, size_t overlapDepthMaster, size_t overlapDepthSlave, bool masterDirichletLast, bool slaveDirichletLast)
{
	std::vector<IndexLayout> vMasterLayouts;
	std::vector<IndexLayout> vSlaveLayouts;
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts);
	c.m_overlapDepthMaster = overlapDepthMaster;
	c.m_overlapDepthSlave = overlapDepthSlave;
	c.m_masterDirichletLast = masterDirichletLast;
	c.m_slaveDirichletLast = slaveDirichletLast;
	bool b = c.calculate();
	newMat.set_layouts(newMasterLayout, newSlaveLayout);
	return b;
}*/

// FAMGLevelCalculator::create_OL2_matrix
//---------------------------------------------------------------------------
/**
 * creates overlap 2 matrix A_OL2,
 * \sa GnerateOverlap
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::create_OL2_matrix()
{
	AMG_PROFILE_FUNC();
	stopwatch SW;

	// get the Overlap 2 matrix
	//-------------------------------

	AMG_PROFILE_BEGIN(AMG_GenerateOverlap2)

	std::vector<IndexLayout> masterLayouts;
	std::vector<IndexLayout> slaveLayouts;
	//GenerateOverlap(A_OL1, A_OL2, OL2MasterLayout, OL2SlaveLayout, 1, 0, true, true);

	AddLayout(OL2MasterLayout, A.get_master_layout());
	AddLayout(OL2SlaveLayout, A.get_slave_layout());

	GenerateOverlap(A, A_OL2, OL2MasterLayout, OL2SlaveLayout, masterLayouts, slaveLayouts, overlapSize, 2);


	// OL1
	AddLayout(OL1MasterLayout, masterLayouts[0]);
	AddLayout(OL1MasterLayout, masterLayouts[1]);
	AddLayout(OL1SlaveLayout, slaveLayouts[0]);
	AddLayout(OL1SlaveLayout, slaveLayouts[1]);



	/*UG_DLOG(LIB_ALG_AMG, 1, "\n\n\nOL1Layout:\n")
	PrintLayout(A_OL2.get_communicator(), OL1MasterLayout, OL1SlaveLayout);

	UG_DLOG(LIB_ALG_AMG, 1, "\n\n\nOL2Layout:\n")
	PrintLayout(A_OL2.get_communicator(), OL2MasterLayout, OL2SlaveLayout);
	*/

	/*UG_LOG("\n\nA layout, level " << level << ":\n")
	PrintLayout(A_OL2.get_communicator(), A.get_master_layout(), A.get_slave_layout());
	UG_LOG("\n\nA1 layout, level " << level << ":\n")
	PrintLayout(A_OL2.get_communicator(), OL1MasterLayout, OL1SlaveLayout);
	UG_LOG("\n\nA2 layout, level " << level << ":\n")
	PrintLayout(A_OL2.get_communicator(), OL1MasterLayout, OL1SlaveLayout);*/

#ifdef UG_DEBUG
	UG_ASSERT(TestLayout(A_OL2.get_communicator(), A.get_master_layout(), A.get_slave_layout()) == true, "A layout wrong, level " << level);
	UG_ASSERT(TestLayout(A_OL2.get_communicator(), OL1MasterLayout, OL1SlaveLayout) == true, "A layout wrong, level " << level);
	UG_ASSERT(TestLayout(A_OL2.get_communicator(), OL1MasterLayout, OL1SlaveLayout) == true, "A layout wrong, level " << level);
#endif

	// get testvectors on newly created indices
	//--------------------------------------------
	// todo: use ONE communicate
	//std::vector<ComPol_VecCopy< Vector<double> > > vecCopyPol;
	//vecCopyPol.resize(m_testvectors.size());

	AMG_PROFILE_NEXT(AMG_get_testvecs_on_OL2);

	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		ComPol_VecCopy< Vector<double> > vecCopyPol;
		m_testvectors[i].resize(A_OL2.num_rows());
		vecCopyPol.set_vector(&m_testvectors[i]);
		pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
		communicator.send_data(OL2MasterLayout, vecCopyPol);
		communicator.receive_data(OL2SlaveLayout, vecCopyPol);
		communicator.communicate();
	}



	// create OLCoarseningSendLayout, OLCoarseningReceiveLayout
	//------------------------------------------------------------
	AMG_PROFILE_NEXT(AMG_create_coarsening_layout);

	// in the pcl, there are only connections between master to slave, and not
	// between all slaves. we add those connections, so that we can send coarsening
	// data from one node to all processes which have this node (as master or as slave).

	// master0 and slave0 nodes can be set coarse by anyone, so we need to send coarsening information
	// master -> slave
	AddLayout(OLCoarseningSendLayout, OL1MasterLayout);
	AddLayout(OLCoarseningReceiveLayout, OL1SlaveLayout);

	// slave -> master
	AddLayout(OLCoarseningSendLayout, OL1SlaveLayout);
	AddLayout(OLCoarseningReceiveLayout, OL1MasterLayout);

	AddConnectionsBetweenSlaves(A_OL2.get_communicator(), OL1MasterLayout, OL1SlaveLayout, OLCoarseningSendLayout, OLCoarseningReceiveLayout);

	//UG_DLOG(LIB_ALG_AMG, 1, "OLCoarseningLayout :\n")
	//PrintLayout(A_OL2.get_communicator(), OLCoarseningSendLayout, OLCoarseningSendLayout);


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
			int pid = slaveLayouts[i].proc_id(iter);
			for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t index = interface.get_element(iter2);
				rating.OLtype[index] |= (1 << (i+1));
				if(i==0)
					rating.m_masterOn[index] = pid;
			}
		}
	}


	// debug: write overlap 2 matrix as debug output

	// 2. get famg helper positions
	//-------------------------------
	AMG_PROFILE_NEXT(AMG_get_famg_helper_pos);

	if(m_famg.m_amghelper.has_positions())
	{
		std::vector<MathVector<3> > &vec2 = m_famg.m_amghelper.positions[level];
		vec2.resize(A_OL2.num_rows());

		ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
		pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
		communicator.send_data(OL2MasterLayout, copyPol);
		communicator.receive_data(OL2SlaveLayout, copyPol);
		communicator.communicate();

		AMG_PROFILE_NEXT(create_OL2_matrix_debug_output);
		if(m_famg.m_writeMatrices)
			WriteMatrixToConnectionViewer(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_A_OL2_L") + ToString(level), ".mat").c_str(),
					A_OL2, &vec2[0], 2);
	}

	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		AMG_PROFILE_NEXT(create_OL2_print_layouts);
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

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	calculate_uncalculated_fine_nodes()
{
	AMG_PROFILE_FUNC();
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(rating[i].is_uncalculated_fine() == false || rating.i_must_assign(i) == false)
			continue;
		if(prolongation_calculated[i])
		{
			calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
			prolongation_calculated[i] = true;
			rating.update_rating(i, possible_parents[i]);
		}

		neighborstruct2 &n = possible_parents[i][0];

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;
			UG_ASSERT(rating[node].is_coarse(), "parent node " << rating.get_original_index(node) << " of node " << rating.get_original_index(i) << " has to be coarse");
			PoldIndices(i, node) = n.parents[j].value;
		}
	}
}


}

#endif // __H__LIB_ALGEBRA__AMG__FAMG_PARALLEL_COARSENING_IMPL_H__
