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

#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif

#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#include "lib_algebra/parallelization/parallel_coloring.h"
#include "../stopwatch.h"
#include "lib_algebra/common/connection_viewer_output.h"
#include "../send_interface.h"
#include "../rsamg/rsamg_impl.h"

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
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelColoring)
	stopwatch SW;

	UG_DLOG(LIB_ALG_AMG, 0, "\ncoloring processor OL graph..."); if(bTiming) SW.start();
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

	m_myColor = ColorProcessorGraph(PN.get_communicator(), pidsOL, processesWithLowerColor, processesWithHigherColor);
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
	UG_DLOG(LIB_ALG_AMG, 0, "\nmy color is " << m_myColor);
	UG_DLOG(LIB_ALG_AMG, 1, "i am connected to");
	for(std::set<int>::iterator it = pidsOL.begin(); it != pidsOL.end(); ++it)
		UG_DLOG(LIB_ALG_AMG, 1, (*it) << " ");
	UG_DLOG(LIB_ALG_AMG, 1, "\nprocesses with lower color:");
	for(size_t i=0; i<processesWithLowerColor.size(); i++)
		UG_DLOG(LIB_ALG_AMG, 1, processesWithLowerColor[i] << " ");
	UG_DLOG(LIB_ALG_AMG, 1, "\nprocesses with higher color:");
	for(size_t i=0; i<processesWithHigherColor.size(); i++)
		UG_DLOG(LIB_ALG_AMG, 1, processesWithHigherColor[i] << " ");
	UG_DLOG(LIB_ALG_AMG, 0, "\n");
}




class FAMGCoarseningCommunicationScheme : public CommunicationScheme<FAMGCoarseningCommunicationScheme >
{
public:
	typedef char value_type;
	FAMGCoarseningCommunicationScheme(FAMGNodes &r) : rating(r)
	{

	}
	char send(int pid, size_t index) const
	{
		int state = rating[index].get_state();
		if(state == FAMG_FINE_RATING || state == FAMG_UNCALCULATED_FINE_RATING || state == FAMG_AGGRESSIVE_FINE_RATING)
			return 0;
		else if(state == FAMG_COARSE_RATING)
			return 1;
		else if(state == 0.0)
			return 2;
		else
		{
			UG_ASSERT(0, "state is " << state << "?");
			return 0;
		}
	}
	void receive(int pid, size_t index, char state)
	{
		if(state == 0)
		{
			//rating.external_set_uncalculated_fine(index);
			rating.set_fine(index);
		}
		else if(state == 1)
			rating.external_set_coarse(index);
		// else ; // nothing

	}

	inline int get_element_size() const
	{
		return sizeof(int);
	}
private:
	FAMGNodes &rating;
};


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
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelRecvCoarsening);

	stopwatch SW;
	if(bTiming) SW.start();

	UG_DLOG(LIB_ALG_AMG, 0, "\nreceive coarsening data from processes with lower color");
	if(processesWithLowerColor.size() == 0)
	{
		UG_DLOG(LIB_ALG_AMG, 0, "no processes with lower color.");
		return;
	}

	FAMGCoarseningCommunicationScheme scheme(rating);
	ReceiveOnInterfaces(PN.get_communicator(), processesWithLowerColor, OLCoarseningReceiveLayout, scheme);

	// issue receive of coarsening data from processes with lower color
	/*IF_DEBUG(LIB_ALG_AMG, 11)
	{
		pcl::ProcessCommunicator lowerPC = A_OL2.get_process_communicator().create_sub_communicator(processesWithLowerColor);
		communicator.enable_communication_debugging(lowerPC);
	}*/

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
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
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelSendCoarsening);

	stopwatch SW; if(bTiming) SW.start();
	UG_DLOG(LIB_ALG_AMG, 0, "\nsend coarsening data to processes...");
	if(processesWithHigherColor.size() == 0)
	{
		UG_DLOG(LIB_ALG_AMG, 0, "no processes with higher color.");
		return;
	}


	FAMGCoarseningCommunicationScheme scheme(rating);
	SendOnInterfaces(PN.get_communicator(), processesWithHigherColor, OLCoarseningSendLayout, scheme);

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");

	/*IF_DEBUG(LIB_ALG_AMG, 11)
	{
		pcl::ProcessCommunicator higherPC = A_OL2.get_process_communicator().create_sub_communicator(processesWithHigherColor);
		communicator.enable_communication_debugging(higherPC);
	}*/
}


// TODO: one "bug" remains: dirichlet nodes, which have only connection to themselfs = 1.0, get afterwards 2.0 (because rows are not additive there)
template<typename matrix_type>
bool GenerateOverlap2(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		IndexLayout &totalMasterLayout, IndexLayout &totalSlaveLayout, std::vector<IndexLayout> &vMasterLayouts, std::vector<IndexLayout> &vSlaveLayouts,
		size_t overlapDepthMaster, size_t overlapDepthSlave, std::vector<size_t> &overlapSize, bool masterDirichletLast, bool slaveDirichletLast,
		ParallelNodes &PN)
{
	PROFILE_FUNC();
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> >
		c(mat, newMat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts, PN);
	c.m_overlapDepthMaster = overlapDepthMaster;
	c.m_overlapDepthSlave = overlapDepthSlave;
	c.m_masterDirichletLast = masterDirichletLast;
	c.m_slaveDirichletLast = slaveDirichletLast;
	bool b = c.calculate();
	overlapSize = c.m_overlapSize;
	return b;
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
	UG_DLOG(LIB_ALG_AMG, 0, "\ncreate OL2 matrix...\n");
#ifdef UG_DEBUG
	//int iDebugLevelMatrixPre = GET_DEBUG_LEVEL(LIB_ALG_MATRIX);
	UG_SET_DEBUG_LEVEL(LIB_ALG_MATRIX, m_famg.iDebugLevelOverlapMatrix);
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelOverlapAMG);
#endif

	AMG_PROFILE_FUNC();
	stopwatch SW;

	// get the Overlap 2 matrix
	//-------------------------------

	AMG_PROFILE_BEGIN(AMG_GenerateOverlap2)

	//GenerateOverlap(A_OL1, A_OL2, OL2MasterLayout, OL2SlaveLayout, 1, 0, true, true);

	AddLayout(OL2MasterLayout, A.get_master_layout());
	AddLayout(OL2SlaveLayout, A.get_slave_layout());

	GenerateOverlap2(A, A_OL2, OL2MasterLayout, OL2SlaveLayout, masterLayouts, slaveLayouts,
			2, 1, overlapSize, false, false, PN);
	// PN.print();

	// OL1
	AddLayout(OL1MasterLayout, masterLayouts[0]);
	AddLayout(OL1SlaveLayout, slaveLayouts[0]);


	/*
	PRINTLAYOUT(PN.get_communicator(), A.get_master_layout(), A.get_slave_layout());
	PRINTLAYOUT(PN.get_communicator(), OL1MasterLayout, OL1SlaveLayout);
	PRINTLAYOUT(PN.get_communicator(), OL2MasterLayout, OL2SlaveLayout);
*/
/*	for(int i=0; i<3; i++)
 	{
		UG_LOG("\n\n\noverlap[" << i << "] layout, level " << level << ":\n");
		PRINTLAYOUT(PN.get_communicator(), masterLayouts[i], slaveLayouts[i]);
	}*/

#ifdef UG_DEBUG
	UG_ASSERT(TestLayout(A.get_process_communicator(), PN.get_communicator(), A.get_master_layout(), A.get_slave_layout()) == true, "A layout wrong, level " << level);
	UG_ASSERT(TestLayout(A.get_process_communicator(), PN.get_communicator(), OL1MasterLayout, OL1SlaveLayout) == true, "A layout wrong, level " << level);
	UG_ASSERT(TestLayout(A.get_process_communicator(), PN.get_communicator(), OL1MasterLayout, OL1SlaveLayout) == true, "A layout wrong, level " << level);
#endif

	// create OLCoarseningSendLayout, OLCoarseningReceiveLayout
	//------------------------------------------------------------
	AMG_PROFILE_NEXT(AMG_create_coarsening_layout);

	// in the pcl, there are only connections between master to slave, and not
	// between all slaves. we add those connections, so that we can send coarsening
	// data from one node to all processes which have this node (as master or as slave).

	// master0 and slave0 nodes can be set coarse by anyone, so we need to send coarsening information
	// master -> slave
	/*AddLayout(OLCoarseningSendLayout, OL1MasterLayout);
	AddLayout(OLCoarseningReceiveLayout, OL1SlaveLayout);

	// slave -> master
	AddLayout(OLCoarseningSendLayout, OL1SlaveLayout);
	AddLayout(OLCoarseningReceiveLayout, OL1MasterLayout);*/


	AddLayout(OLCoarseningSendLayout, masterLayouts[0]);
	AddLayout(OLCoarseningSendLayout, slaveLayouts[0]);
	AddLayout(OLCoarseningSendLayout, masterLayouts[1]);
	AddLayout(OLCoarseningSendLayout, slaveLayouts[1]);

	AddLayout(OLCoarseningReceiveLayout, slaveLayouts[0]);
	AddLayout(OLCoarseningReceiveLayout, masterLayouts[0]);
	AddLayout(OLCoarseningReceiveLayout, slaveLayouts[1]);
	AddLayout(OLCoarseningReceiveLayout, masterLayouts[1]);

	AddConnectionsBetweenSlaves(PN.get_communicator(), OL1MasterLayout, OL1SlaveLayout, OLCoarseningSendLayout, OLCoarseningReceiveLayout);

	//UG_DLOG(LIB_ALG_AMG, 1, "OLCoarseningLayout :\n")
	//PrintLayout(PN.get_communicator(), OLCoarseningSendLayout, OLCoarseningSendLayout);


	size_t N = A_OL2.num_rows();
	rating.create(N);

	/*
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
	}*/

	// debug: write overlap 2 matrix as debug output

	// 2. get famg helper positions
	//-------------------------------
	AMG_PROFILE_NEXT(AMG_get_famg_helper_pos);

	if(pcl::GetNumProcesses() > 1 && m_famg.m_amghelper.has_positions())
	{
		m_famg.m_amghelper.update_overlap_positions(level, PN.get_communicator(), OL2MasterLayout, OL2SlaveLayout, A_OL2.num_rows());

		AMG_PROFILE_NEXT(create_OL2_matrix_debug_output);
		if(m_famg.m_writeMatrices)
				WriteMatrixToConnectionViewer((m_famg.m_writeMatrixPath + std::string("AMG_A_OL2_L") + ToString(level) + ".mat").c_str(),
						A_OL2, &m_famg.m_amghelper.positions[level][0], 2);
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

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");

//	UG_SET_DEBUG_LEVEL(LIB_ALG_MATRIX, iDebugLevelMatrixPre);
}

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	calculate_uncalculated_fine_nodes()
{
	AMG_PROFILE_FUNC();
	size_t iNrOfUncalculated=0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(rating[i].is_uncalculated_fine() == false || rating.i_must_assign(i) == false)
			continue;
		iNrOfUncalculated ++;
		if(prolongation_calculated[i] == false)
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
	UG_LOG("calculated " << iNrOfUncalculated << " of uncalculated nodes.\n");
}


}

#endif // __H__LIB_ALGEBRA__AMG__FAMG_PARALLEL_COARSENING_IMPL_H__
