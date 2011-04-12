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
		pidsOL.insert(OLCoarseningSendLayout.proc_id(iter));
	for(IndexLayout::iterator iter = OLCoarseningReceiveLayout.begin(); iter != OLCoarseningReceiveLayout.end(); ++iter)
		pidsOL.insert(OLCoarseningReceiveLayout.proc_id(iter));

	m_myColor = ColorProcessorGraph(communicator, pidsOL, processesWithLowerColor, processesWithHigherColor);
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
	UG_DLOG(LIB_ALG_AMG, 1, "my color is " << m_myColor);
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
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();

	// issue receive of coarsening data from processes with lower color
	if(processesWithLowerColor.size() > 0)
	{
		stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 1, "\nWaiting for processes "); if(bTiming) SW.start();
		stdvector< stdvector<char> > coarseNodes;
		coarseNodes.resize(processesWithLowerColor.size());

		for(size_t i=0; i<processesWithLowerColor.size(); i++)
		{
			int pid = processesWithLowerColor[i];

			size_t s = OLCoarseningReceiveLayout.interface(pid).size();
			coarseNodes[i].resize(s, -1);
			UG_DLOG(LIB_ALG_AMG, 1, pid << ", awaiting  " << s << " bytes.");
			communicator.receive_raw(pid, &coarseNodes[i][0], s);
		}
		UG_DLOG(LIB_ALG_AMG, 1, "which have higher color to receive coarse nodes... ");
		communicator.communicate();
		UG_DLOG(LIB_ALG_AMG, 1, "done. processing data...");

		// set nodes coarse
		for(size_t i=0; i<processesWithLowerColor.size(); i++)
		{
			int pid = processesWithLowerColor[i];
			UG_DLOG(LIB_ALG_AMG, 3, "\nfrom processor " << pid << ":\n");

			int j=0;
			IndexLayout::Interface &interface = OLCoarseningReceiveLayout.interface(pid);
			IndexLayout::Interface &nextMasterInterface = nextLevelMasterLayout.interface(pid);
			IndexLayout::Interface &nextSlaveInterface = nextLevelSlaveLayout.interface(pid);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t index = interface.get_element(iter);
				//UG_DLOG(LIB_ALG_AMG, 3, (int)coarseNodes[i][j] << " ");
				UG_DLOG(LIB_ALG_AMG, 3, " " << index);
				if(rating.is_master(index))
					UG_DLOG(LIB_ALG_AMG, 3, " master");

				if(coarseNodes[i][j++])
				{
					rating.external_set_coarse(index);
					int newIndex = rating.newIndex[index];
					UG_ASSERT(newIndex != -1, "");

					// add node to next Level Interface
					if(rating.is_master(index))
						nextMasterInterface.push_back(newIndex);
					else if(rating.is_master_on(index, pid))
						nextSlaveInterface.push_back(newIndex);

					UG_DLOG(LIB_ALG_AMG, 3, " coarse ");
				}
				else
				{
					// todo: das ist falsch, da macht der nacher coarse nodes draus, und das ist doof
					// man mŸsste vielleicht hier direkt die interpolation ausrechnen. dann mŸsste man sie aber
					// mitschicken?!? nochmal genau Ÿberlegen, wie die parallelisierung funktionieren soll
					rating.set_uninterpolateable(index);
					UG_DLOG(LIB_ALG_AMG, 3, " fine ");
				}
				UG_DLOG(LIB_ALG_AMG, 3, "\n");
			}
		}
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 3, "took " << SW.ms() << " ms");
	}
	else
		UG_DLOG(LIB_ALG_AMG, 3, "\nno processes with lower color.")

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
	AMG_PROFILE_FUNC();
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();

	stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 1, "\nsend coarsening data to processes\n"); if(bTiming) SW.start();
	for(size_t i=0; i<processesWithHigherColor.size(); i++)
	{
		int pid = processesWithHigherColor[i];
		UG_DLOG(LIB_ALG_AMG, 1, "Process " << pid << ":\n");
		BinaryStream s;

		IndexLayout::Interface &interface = OLCoarseningSendLayout.interface(pid);
		IndexLayout::Interface &nextMasterInterface = nextLevelMasterLayout.interface(pid);
		IndexLayout::Interface &nextSlaveInterface = nextLevelSlaveLayout.interface(pid);

		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t index = interface.get_element(iter);
			char bCoarse = rating[index].is_coarse();
			UG_DLOG(LIB_ALG_AMG, 3, index << (bCoarse ? " is coarse " : " is fine "));
			Serialize(s, bCoarse);
			if(rating.is_master(index))
				UG_DLOG(LIB_ALG_AMG, 3, "master");
			UG_DLOG(LIB_ALG_AMG, 3, "\n");
			if(bCoarse)
			{
				int newIndex = rating.newIndex[index];
				UG_ASSERT(newIndex != -1, "");

				// add node to next Level Interface
				if(rating.is_master(index))
					nextMasterInterface.push_back(newIndex);
				else if(rating.is_master_on(index, pid))
					nextSlaveInterface.push_back(newIndex);
			}
		}

		UG_DLOG(LIB_ALG_AMG, 1, "sending " << s.size() << " of data to pid " << pid << "\n");
		communicator.send_raw(pid, s.buffer(), s.size(), true);
	}
	UG_DLOG(LIB_ALG_AMG, 1, "with higher color...")

	communicator.communicate();
	UG_DLOG(LIB_ALG_AMG, 1, "done.");
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
}


// FAMGLevelCalculator::add_connections_between_slave_nodes
//---------------------------------------------------------------------------
/** adds connections between slave nodes to OLCoarseningSendLayout/OLCoarseningReceiveLayout
 * in the pcl, there are only connections between master to slave, and not
 * between all slaves. we add those connections, so that we can send coarsening
 * data from one node to all processes which have this node (as master or as slave).
  *
 */



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

	// 1. get the Overlap 2 matrix
	//-------------------------------
	UG_DLOG(LIB_ALG_AMG, 1, "\nGenerate Overlap 2..."); if(bTiming) SW.start();

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

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// 2. get data on newly created indices
	//--------------------------------------------
	// use ONE communicate

	UG_DLOG(LIB_ALG_AMG, 1, "\nGet data on Overlap 2 nodes..."); if(bTiming) SW.start();

	ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
	pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
	communicator.send_data(totalMasterLayout, copyPol);
	communicator.receive_data(totalSlaveLayout, copyPol);
	communicator.communicate();

	//std::vector<ComPol_VecCopy< Vector<double> > > vecCopyPol;
	//vecCopyPol.resize(m_testvectors.size());
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		ComPol_VecCopy< Vector<double> > vecCopyPol;
		m_testvectors[i].resize(A_OL2.num_rows());
		vecCopyPol.set_vector(&m_testvectors[i]);
		communicator.send_data(totalMasterLayout, vecCopyPol);
		communicator.receive_data(totalSlaveLayout, vecCopyPol);
		communicator.communicate();
	}



	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// create OLCoarseningSendLayout, OLCoarseningReceiveLayout
	//------------------------------------------------------------

	// in the pcl, there are only connections between master to slave, and not
	// between all slaves. we add those connections, so that we can send coarsening
	// data from one node to all processes which have this node (as master or as slave).

	// master -> master
	AddLayout(OLCoarseningSendLayout, masterLayouts[0]);
	AddLayout(OLCoarseningReceiveLayout, slaveLayouts[0]);

	// slave -> master
	AddLayout(OLCoarseningSendLayout, slaveLayouts[0]);
	AddLayout(OLCoarseningReceiveLayout, masterLayouts[0]);

	AddConnectionsBetweenSlaves(A_OL2.get_communicator(), masterLayouts[0], slaveLayouts[0], OLCoarseningSendLayout, OLCoarseningReceiveLayout);


	UG_DLOG(LIB_ALG_AMG, 1, "OLCoarseningLayout :\n")
	//PrintLayout(A_OL2.get_communicator(), OLCoarseningSendLayout, OLCoarseningSendLayout);

/*	UG_DLOG(LIB_ALG_AMG, 1, "OLCoarseningSendLayout :\n");
	PrintLayout(OLCoarseningSendLayout);

	UG_DLOG(LIB_ALG_AMG, 1, "OLCoarseningRecvLayout :\n");
	PrintLayout(OLCoarseningSendLayout);*/

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
				rating.OLtype[index] |= (1 << (i+1));
				if(i==0)
					rating.m_masterOn[index] = slaveLayouts[i].proc_id(iter);
			}
		}
	}

	//rating.print_OL_types();


	PROFILE_BEGIN(create_OL2_matrix_debug_output)

	// debug: write overlap 2 matrix as debug output
	if(m_famg.m_writeMatrices)
		WriteMatrixToConnectionViewer(GetProcFilename(m_famg.m_writeMatrixPath, "A_OL2", ".mat").c_str(),
			A_OL2, &vec2[0], 2);

	IF_DEBUG(LIB_ALG_AMG, 4)
	{
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
	PROFILE_END();
}

}

#endif // __H__LIB_ALGEBRA__AMG__FAMG_PARALLEL_COARSENING_IMPL_H__
