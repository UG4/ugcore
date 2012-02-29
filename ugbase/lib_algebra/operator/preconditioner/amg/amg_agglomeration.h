#include "agglomeration.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/global_layout.h"
#endif

namespace ug
{

#ifdef UG_PARALLEL
template<typename TAlgebra>
bool AMGBase<TAlgebra>::gather_vertical(vector_type &vec, vector_type &collectedVec, size_t level,
		ParallelStorageType type)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	pcl::ParallelCommunicator<IndexLayout> &com = A.communicator();

	if(!levels[level]->bHasBeenMerged)
	{
		if(&vec == &collectedVec) return true;
		else
		{
			for(size_t i=0; i<vec.size(); i++)
				collectedVec[i] = vec[i]; // don't use that
		}
	}
	else
	{
		if(m_agglomerateLevel == level)
		{
			ComPol_VecAdd<vector_type > compolAdd(&collectedVec, &vec);
			com.send_data(agglomerateSlaveLayout, compolAdd);
			com.communicate();
		}
		else
		{

			UG_ASSERT(&vec != &collectedVec, "");
			collectedVec.set(0.0);
			for(size_t i=0; i<vec.size(); i++)
				collectedVec[i] = vec[i];

			UG_ASSERT(vec.has_storage_type(type), "");
			if(type == PST_ADDITIVE)
			{
				ComPol_VecAdd<vector_type > compolAdd(&collectedVec, &vec);
				com.receive_data(L.agglomerateMasterLayout, compolAdd);
				com.communicate();
				collectedVec.set_storage_type(PST_ADDITIVE);
			}
			else if(type == PST_CONSISTENT)
			{
				ComPol_VecCopy<vector_type > compolCopy(&collectedVec, &vec);
				com.receive_data(L.agglomerateMasterLayout, compolCopy);
				com.communicate();
				collectedVec.set_storage_type(PST_CONSISTENT);
			}
			else { UG_ASSERT(0, "unsupported."); }
		}
	}
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::broadcast_vertical(vector_type &vec, vector_type &collectedVec, size_t level,
		ParallelStorageType type)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	pcl::ParallelCommunicator<IndexLayout> &com = A.communicator();

	if(!levels[level]->bHasBeenMerged)
	{
		if(&vec == &collectedVec) return true;
		else
		{
			for(size_t i=0; i<vec.size(); i++)
				vec[i] = collectedVec[i]; // don't use that
		}
	}
	else
	{
		if(m_agglomerateLevel == level)
		{
			ComPol_VecCopy<vector_type> compolCopy(&vec, &collectedVec);
			com.receive_data(agglomerateSlaveLayout, compolCopy);
			com.communicate();
			vec.set_storage_type(type);
			if(type == PST_ADDITIVE)
			{
				vec.set_storage_type(PST_ADDITIVE);
				SetLayoutValues(&vec, agglomerateSlaveLayout, 0.0); //!!!
			}
			else if(type == PST_CONSISTENT) {	}
			else { UG_ASSERT(0, "unsupported."); }
		}
		else
		{
			UG_ASSERT(&vec != &collectedVec, "");
			for(size_t i=0; i<vec.size(); i++)
				vec[i] = collectedVec[i];
			UG_ASSERT(collectedVec.has_storage_type(type), "");
			vec.set_storage_type(type);

			ComPol_VecAdd<vector_type > compolCopy(&vec, &collectedVec);
			com.send_data(L.agglomerateMasterLayout, compolCopy);
			com.communicate();
		}
	}
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::isMergingSlave(size_t level)
{
	return levels[level]->bHasBeenMerged && m_agglomerateLevel == level;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::isMergingMaster(size_t level)
{
	return levels[level]->bHasBeenMerged && m_agglomerateLevel != level;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::isNotMerging(size_t level)
{
	return levels[level]->bHasBeenMerged && m_agglomerateLevel != level;
}
#endif


template<typename TAlgebra>
bool AMGBase<TAlgebra>::add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
#ifdef UG_PARALLEL
	if(levels[level]->bLevelHasMergers)
	{
		//matrix_operator_type &collA = *L.pAgglomeratedA;

		gather_vertical(d, L.collD, level, PST_ADDITIVE);

		bool b=true;
		if(!isMergingSlave(level))
		{
			L.collC.set(0.0);
			b = add_correction_and_update_defect2(L.collC, L.collD, *L.pAgglomeratedA, level);
		}
		broadcast_vertical(c, L.collC, level, PST_CONSISTENT);
		broadcast_vertical(d, L.collD, level, PST_ADDITIVE);
		return b;
	}
	else
		return add_correction_and_update_defect2(c, d, A, level);
#else
	return add_correction_and_update_defect2(c, d, A, level);
#endif

}

/*
 * überall brauche ich
 * - int iMergeLevel
 *
 * std::vector<ProcessCommunicator> m_vProcessCommunicator;
 * std::vector<IndexLayout> m_vAgglomerateMasterLayout
 * std::vector<vector_type> collC, collD;
 * std::vector<matrix_type> collA;
 * std::vector<bool> bHasBeenMerged;
 *
 * IndexLayout agglomerateSlaveLayout;
 * size_t m_agglomerateLevel;
 *
 * auf jedem level, bei dem ich mitmach
 * - einen process_communicator (für participate usw)
 * - ein IndexLayout m_vAgglomerateMasterLayout
 * - zusätzliche Vektoren collC, collD
 * - zusätzliche Matrix collA
 * - ein flag bHasBeenMerged
 *
 * nur auf dem gröbsten, auf dem ich mitmache (level == iMergeLevel)
 * - ein IndexLayout agglomerateSlaveLayout
 * - einen process_communicator (für participate usw)
 *
 * A[0] - 1024 entries. coarsen
 * A[1] - 512 entries. Merge with processor 2 -> A[1] has now 1024 entries
 * A[2] - 512 entries. Merge to processor 3. finished.
 *
 *
 *
 * agglomerate()
 *  ...
 *  // (ProcessCommunicator L.agglomeratedPC)
 *  if bMergesWithAnotherProcessor
 *  	L.agglomeratedPC = A.get_process_communicator().create_sub_communicator(false)
 *  else
 *  	L.agglomeratedPC = A.get_process_communicator().create_sub_communicator(true)
 *
 *
 *	d additiv, c consistent
 *
 * MGC(c, d, level)
 * if level == m_aggloLevel
 * 		sendtoMaster(c, d, level)
 * else
 * 		if bHasMerged[level] == true
 * 			getFromSlaves(c, collC[level], d, collD[level], level)
 * 			c = collC[level]; d = collD[level];
 * 		endif
 * 		smooth(c, d, A[level])
 *
 *
 * !!! AUFPASSEN AUFPASSEN AUFPASSEN !!!
 *  mit additiv und konsistent usw!!!!!
 * !!! AUFPASSEN AUFPASSEN AUFPASSEN !!!
 * */
#ifdef UG_PARALLEL
template<typename TAlgebra>
bool AMGBase<TAlgebra>::agglomerate(size_t level)
{

	//UG_ASSERT(0, "not working");
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	if(L.m_levelInformation.get_nr_of_nodes_min() > m_minNodesOnOneProcessor
			|| A.process_communicator().size() == 1)
	{
		L.bHasBeenMerged = false;
		L.pAgglomeratedA = L.pA;
		L.bLevelHasMergers = false;
		return true;
	}



	L.bLevelHasMergers=true;

	// 1. create global Algebra IDs

	pcl::ParallelCommunicator<IndexLayout> &communicator = A.communicator();
	pcl::ProcessCommunicator &pc = A.process_communicator();

	ParallelNodes PN(A.communicator(), A.master_layout(), A.slave_layout(), A.num_rows());

	// 1. ein prozessor erhält alle daten
	// 2. ein maß definiert, wie leicht man agglomerieren kann:
	// - reduktion des interfaces

	BinaryBuffer rstream;


	/*UG_LOG("ProcessCommunicator is ");
	for(size_t i=0; i<pc.size(); i++) UG_LOG(pc.get_proc_id(i) << " ");
	UG_LOG("\n");*/

	if(pcl::GetProcRank() == pc.get_proc_id(0))
	{
		//UG_LOG("Receiving interface information.\n");

		// receive interface information and
		typedef std::map<int, BinaryBuffer> BufferMap;
		std::vector<BinaryBuffer> receivepack;
		receivepack.resize(pc.size());
		for(size_t i=1; i<pc.size(); i++)
			communicator.receive_raw(pc.get_proc_id(i), receivepack[i]);

		communicator.communicate();

		// the participating processors can be something like 5, 3, 22, 7
		// so we need to map 5 -> 0, 3 -> 1, 22 -> 2 and 7 -> 3
		// global -> "local"
		// so that agglomeration algorithms are easier to implement.

		std::map<int, int> globalToLocalPID;
		for(size_t i=0; i<pc.size(); i++)
		{
			globalToLocalPID[pc.get_proc_id(i)] = i;
			//UG_LOG("global " << pc.get_proc_id(i) << " local " << i << "\n");
		}

		std::vector<size_t> sizes;
		sizes.resize(pc.size());
		std::vector<std::map<int, size_t> > connections;
		connections.resize(pc.size());

		sizes[0] = A.num_rows();
		for(IndexLayout::iterator iter = A.master_layout().begin(); iter != A.master_layout().end(); ++iter)
		{
			size_t s = A.master_layout().interface(iter).size();
			int pid = A.master_layout().proc_id(iter);
			int localID = globalToLocalPID[pid];
			connections[0][localID] += s;
			connections[localID][0] += s;
		}

		for(size_t i=1; i<pc.size(); i++)
		{
			BinaryBuffer &stream = receivepack[i];

			size_t s; int pid;
			Deserialize(stream, s);
			sizes[i] = s;
			while(!stream.eof())
			{
				Deserialize(stream, pid);
				Deserialize(stream, s);
				int localID = globalToLocalPID[pid];
				connections[localID][i] += s;
				connections[i][localID] += s;
			}
		}

		// output info
		for(size_t i=0; i<pc.size(); i++)
		{
			/*UG_LOG("Processor nr. " << i << " has PID " << pc.get_proc_id(i) << ", " << sizes[i] << " nodes and connections to ");
			for(std::map<int, size_t>::iterator it = connections[i].begin(); it != connections[i].end(); ++it)
				UG_LOG("processor nr. " << (*it).first << " PID " << pc.get_proc_id((*it).first) << " (" << (*it).second << " connections) ")
			UG_LOG("\n");*/
		}

		std::vector<std::vector<int> > mergeWith(pc.size());
		EasyAgglomeration(sizes, connections, mergeWith, m_minNodesOnOneProcessor, m_preferredNodesOnOneProcessor);


		for(size_t i=0; i<pc.size(); i++)
		{
			// send each processor
			// size
			// if size > 0
			// 		size = nr of processors which merge with you
			//		processors [1..n]
			// else
			//		master processor to merge to
			// endif
			// connected neighbors which merge with others

			BinaryBuffer stream;

			std::vector<int> mergeWithGlobal;
			for(size_t j=0; j<mergeWith[i].size(); j++)
				mergeWithGlobal.push_back(pc.get_proc_id(mergeWith[i][j]));
			Serialize(stream, mergeWithGlobal);

			std::map<int, int> merge; // maps old -> new

			for(std::map<int, size_t>::iterator it = connections[i].begin(); it != connections[i].end(); ++it)
			{
				int j = it->first;
				if(mergeWith[j][0] != j)
					merge[pc.get_proc_id(j)] = pc.get_proc_id(mergeWith[j][0]);
			}

			Serialize(stream, merge);

			if(i == 0)
				rstream = stream;
			else
				communicator.send_raw(pc.get_proc_id(i),
					stream.buffer(), stream.write_pos(), false);

		}
		communicator.communicate();
	}
	else
	{
		//UG_LOG("Sending interface information.\n");
		BinaryBuffer sstream;
		Serialize(sstream, A.num_rows());
		for(IndexLayout::iterator iter = A.master_layout().begin(); iter != A.master_layout().end(); ++iter)
		{
			int pid = A.master_layout().proc_id(iter);
			size_t s = A.master_layout().interface(iter).size();
			Serialize(sstream, pid);
			Serialize(sstream, s);
		}
		communicator.send_raw(pc.get_proc_id(0), sstream.buffer(), sstream.write_pos(), false);
		communicator.communicate();

		communicator.receive_raw(pc.get_proc_id(0), rstream);
		communicator.communicate();
	}

	// receive informations
	// 1. - should i merge and i am master? -> list of merging processors
	//    - should i merge and i am slave? -> master processor
	//    - i shouldnt merge
	// 2. which processors are active afterwards

	matrix_operator_type M;

	std::vector<int> mergeWith;
	Deserialize(rstream, mergeWith);

	std::map<int, int> merge; // maps old -> new
	Deserialize(rstream, merge);

	/*UG_LOG("Neighbor merge information: ");
	for(std::map<int, int>::iterator it = merge.begin(); it != merge.end(); ++it) UG_LOG(it->first << " -> " << it->second);
	UG_LOG("\n");*/

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.master_layout(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.slave_layout(), PN);

	if(mergeWith.size() > 1)
	{
		L.bHasBeenMerged = true;
		UG_LOG("I am the merging father for pids ");
		for(size_t i=1; i<mergeWith.size(); i++) UG_LOG(mergeWith[i] << " ");
		UG_LOG("\n");
		mergeWith.erase(mergeWith.begin());

		PrintGlobalLayout(globalMasterLayout, "globalmasterLayout before");
		PrintGlobalLayout(globalSlaveLayout, "globalSlaveLayout before");

		ReceiveMatrix(A, L.collectedA, L.agglomerateMasterLayout, mergeWith, PN);
		ReceiveGlobalLayout(communicator, mergeWith, globalMasterLayout, globalSlaveLayout);

		PrintGlobalLayout(globalMasterLayout, "received globalmasterLayout");
		PrintGlobalLayout(globalSlaveLayout, "received globalSlaveLayout");

		MergeGlobalLayout(globalMasterLayout, merge);
		MergeGlobalLayout(globalSlaveLayout, merge);

		PrintGlobalLayout(globalMasterLayout, "AM globalmasterLayout");
		PrintGlobalLayout(globalSlaveLayout, "AM globalSlaveLayout");

		CreateLayoutFromGlobalLayout(L.masterLayout2, globalMasterLayout, PN);
		CreateLayoutFromGlobalLayout(L.slaveLayout2, globalSlaveLayout, PN);

		L.agglomeratedPC = A.process_communicator().create_sub_communicator(true);

		L.collectedA.set_layouts(L.masterLayout2, L.slaveLayout2);
		L.collectedA.set_process_communicator(L.agglomeratedPC);
		L.pAgglomeratedA  = &L.collectedA;

		L.collC.resize(L.collectedA.num_rows());
		SetParallelVectorAsMatrix(L.collC, L.collectedA, PST_CONSISTENT);
		L.collD.resize(L.collectedA.num_rows());
		SetParallelVectorAsMatrix(L.collD, L.collectedA, PST_ADDITIVE);

		/*PRINTPC(L.agglomeratedPC);
		UG_LOG("L.masterLayout2\n");
		PrintLayout(L.masterLayout2);
		UG_LOG("L.slaveLayout2\n");
		PrintLayout(L.slaveLayout2);
		UG_LOG("level = " << level << ". old size of matrix at this level: " <<L.collectedA.num_rows() << "\n");*/

		if(m_amghelper.has_positions())
			m_amghelper.receive_agglomerate_positions(level, A.communicator(), L.agglomerateMasterLayout,
				L.collectedA.num_rows());
		//merging_as_master(L);
	}
	else
	{
		int pid = mergeWith[0];
		if(pid == pcl::GetProcRank())
		{
			L.bHasBeenMerged = false;
			//UG_LOG("Not merging.\n");
			L.agglomeratedPC = A.process_communicator().create_sub_communicator(true);

			MergeGlobalLayout(globalMasterLayout, merge);
			MergeGlobalLayout(globalSlaveLayout, merge);
			CreateLayoutFromGlobalLayout(L.masterLayout2, globalMasterLayout, PN);
			CreateLayoutFromGlobalLayout(L.slaveLayout2, globalSlaveLayout, PN);

			PrintGlobalLayout(globalMasterLayout, "AM globalmasterLayout");
					PrintGlobalLayout(globalSlaveLayout, "AM globalSlaveLayout");

			//PRINTPC(L.agglomeratedPC);
			/*UG_LOG("others merged, might have to change Layouts:");
			UG_LOG("masterLayout2\n");
			PrintLayout(L.masterLayout2);
			UG_LOG("slaveLayout2\n");
			PrintLayout(L.masterLayout2);*/

			// todo: this is some waste because we're just changing layouts.
			// change vecs tooo....
			L.collectedA = A;
			L.collectedA.set_layouts(L.masterLayout2, L.slaveLayout2);
			L.collectedA.set_process_communicator(L.agglomeratedPC);
			L.pAgglomeratedA  = &L.collectedA;

			L.collC.resize(L.collectedA.num_rows());
			SetParallelVectorAsMatrix(L.collC, L.collectedA, PST_CONSISTENT);
			L.collD.resize(L.collectedA.num_rows());
			SetParallelVectorAsMatrix(L.collD, L.collectedA, PST_ADDITIVE);
		}
		else
		{
			L.bHasBeenMerged = true;
			UG_LOG("MERGING to " << pid << "\n");

			SendMatrix(A, agglomerateSlaveLayout, pid, PN);

			MergeGlobalLayout(globalMasterLayout, merge);
			MergeGlobalLayout(globalSlaveLayout, merge);

			// remove interfaces to pid, so we do not unnecessary send them
			globalMasterLayout.erase(pid);
			globalSlaveLayout.erase(pid);

			PrintGlobalLayout(globalMasterLayout, "send globalmasterLayout");
			PrintGlobalLayout(globalSlaveLayout, "send globalSlaveLayout");
			SendGlobalLayout(communicator, globalMasterLayout, globalSlaveLayout, pid);

			L.agglomeratedPC = A.process_communicator().create_sub_communicator(false);
			/*PRINTPC(L.agglomeratedPC);
			UG_LOG("agglomerateSlaveLayout\n");
			PrintLayout(agglomerateSlaveLayout);
			UG_LOG("MERGING done \n");*/

			/*L.masterLayout.clear();
			L.slaveLayout.clear();
			A.set_layouts(L.masterLayout, L.slaveLayout);
			A.set_process_communicator(L.agglomeratedPC);*/
			m_agglomerateLevel = level;

			L.pAgglomeratedA  = NULL;
			if(m_amghelper.has_positions())
				m_amghelper.send_agglomerate_positions(level, A.communicator(), agglomerateSlaveLayout);
			//merging_as_slave(L);
		}
	}
	return true;
}
#endif

}
