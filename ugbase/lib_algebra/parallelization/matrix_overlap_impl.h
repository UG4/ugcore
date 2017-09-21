/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_matrix_overlap_impl
#define __H__UG_matrix_overlap_impl

#include "algebra_layouts.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"
#include "parallelization_util.h"

namespace ug{

template <class TLayout>
void GetLayoutTargetProcs(std::vector<int>& procsOut, const TLayout& layout)
{	
	procsOut.clear();
	typedef typename TLayout::const_iterator const_iter_t;
	for(const_iter_t iter = layout.begin(); iter != layout.end(); ++iter)
	{
		procsOut.push_back(layout.interface(iter).get_target_proc());
	}
}



/// Highly specialized communication policy for matrix overlap creation
/** \note	This policy is only intended to be used for slave->master communication.
 * \note 	This policy is only used for internal implementation for overlap creation,
 * 			e.g. in CreateOverlap.
 *
 * After communicating from slave->master, call post_process() to create the actual
 * overlap. GlobalIDs will also be updated for new entries.
 *
 * \warning	Please make sure that the matrix and the array of globalIDs passed to
 *			the constructor exist until the instance of this class is destroyed.
 * \sa CreateOverlap*/
template <class TMatrix>
class ComPol_MatCreateOverlap
	: public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
		typedef typename TMatrix::value_type block_t;
	///	Constructor setting the vector
	/**
	 * vGlobalID must have size >= mat.num_rows()
	 */
		ComPol_MatCreateOverlap(TMatrix& rMat, AlgebraIDVec& vGlobalID)
			: m_mat(rMat), m_globalIDs(vGlobalID)
		{
			UG_COND_THROW(vGlobalID.size() < m_mat.num_rows(), "Not enough GlobalIDs");
		//	fill the map global->local
			GenerateAlgebraIDHashList(m_algIDHash, m_globalIDs);
		}

	///	writes the interface values into a buffer that will be sent
		virtual bool collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_collect, "algebra parallelization");
			typedef typename TMatrix::row_iterator row_iterator;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	count number of row entries
				const row_iterator rowEnd = m_mat.end_row(index);
				size_t numRowEntry = 0;
				for(row_iterator it_k = m_mat.begin_row(index); it_k != rowEnd; ++it_k)
					numRowEntry++;

			//	write number of row entries to stream
				Serialize(buff, numRowEntry);

			//	write entries and global id to stream
				for(row_iterator it_k = m_mat.begin_row(index); it_k != rowEnd; ++it_k)
				{
					const size_t k = it_k.index();
					block_t& a_ik = it_k.value();

				//	write global entry to buffer
					Serialize(buff, m_globalIDs[k]);

				//	write entry into buffer
					Serialize(buff, a_ik);
				}
			}

		///	done
			return true;
		}


	///	writes values from a buffer into the interface
		virtual bool extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_extract, "algebra parallelization");

		//	we'll read global ids into this variable
			AlgebraID gID;

		//	we'll read blocks into this var
			block_t block;

			const int targetProc = interface.get_target_proc ();

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	read the number of connections
				size_t numConnections = 0;
				Deserialize(buff, numConnections);

			//	read each connection
				for(size_t i_conn = 0; i_conn < numConnections; ++i_conn){
					Deserialize(buff, gID);
					Deserialize(buff, block);

				//	if gID exists on this process, then set the connection to
				//	the new value.
					size_t conInd;
					if(m_algIDHash.get_entry(conInd, gID)){
						m_mat(index, conInd) += block;
					}
					else {
						const ExtCon ec(index, m_globalIDs[index], gID, targetProc);
						m_recvNewCons [ec] = block;
						m_recvNewIDs.insert (gID);
					}
				}
			}

		///	done
			return true;
		}


		/// After communication is done, this method should be called to create the overlap
		virtual void post_process()
		{
		//	create new master overlap
			using namespace std;

		//	we create a new layout, since the old one may be shared between many
		//	different vectors and matrices. H-Master and H-Slave layouts are still
		//	the same.
			SmartPtr <AlgebraLayouts> newLayout = make_sp(new AlgebraLayouts);
			*newLayout = *m_mat.layouts();
			newLayout->enable_overlap(true);
			m_mat.set_layouts (newLayout);

			const size_t oldSize = m_mat.num_rows();
			const size_t numNewInds = m_recvNewIDs.size();
			const size_t newSize = oldSize + numNewInds;

		//	add new entries to the algebra hash and to the globalID array
			{
				m_globalIDs.reserve(newSize);
				size_t i = oldSize;
				for(set<AlgebraID>::iterator iter = m_recvNewIDs.begin();
				    iter != m_recvNewIDs.end(); ++iter, ++i)
				{
					m_algIDHash.insert(*iter, i);
					m_globalIDs.push_back(*iter);
				}
			}

			if(newSize != oldSize) {
			//	For each new DoF we set a DirichletRow
				m_mat.resize_and_keep_values (newSize, newSize);
				for(size_t i = oldSize; i < newSize; ++i){
					SetDirichletRow(m_mat, i);
				}
			}


			vector<int> slaveProcs; // process ranks of processes with associated slave interfaces
			GetLayoutTargetProcs(slaveProcs, newLayout->master());

		//	the following buffers and vectors are used to collect globalIDs of
		//	newly created entries sorted by the slave process from which they
		//	were received.
			BinaryBuffer sendBuf;
			// size of the message for the i-th process in slaveProcs
			vector<int> msgSizeForSlaveProcs(slaveProcs.size(), 0);
			int slaveInd = -1;

		//	create new master-overlap and add connections to matrix
			{
				IndexLayout::Interface* itfc = NULL;
				vector<bool> added;// keeps track of whether an index was already added.

				for (typename map<ExtCon, block_t>::iterator iter = m_recvNewCons.begin();
				     iter != m_recvNewCons.end(); ++iter)
				{
					const ExtCon& curExtCon = iter->first;
					const int targetProc = curExtCon.conProc;
					if (!itfc || targetProc != itfc->get_target_proc()){
						itfc = &newLayout->master_overlap().interface(targetProc);
						added.clear();
						added.resize(numNewInds, false);

					//	find the index of the corresponding slave proc
						slaveInd = -1;
						for(size_t islave = 0; islave < slaveProcs.size(); ++islave){
							if(slaveProcs[islave] == targetProc){
								slaveInd = int(islave);
								break;
							}
						}
						UG_COND_THROW(slaveInd == -1, "slaveProcs doas not contain referenced slave rank");
					}

					size_t toInd;
					UG_COND_THROW(!m_algIDHash.get_entry(toInd, curExtCon.toID),
					              "Expected AlgebraID not found in Hash");

					if(!added[toInd - oldSize]){
						itfc->push_back(toInd);
						const size_t oldWritePos = sendBuf.write_pos();
						Serialize(sendBuf, curExtCon.toID);
						msgSizeForSlaveProcs[slaveInd] += int(sendBuf.write_pos() - oldWritePos);
					}
					added[toInd - oldSize] = true;

					m_mat(curExtCon.fromInd, toInd) += iter->second;
				}
			}

			if(!msgSizeForSlaveProcs.empty()){
				UG_LOG("msgSizeForSlaveProcs: " << msgSizeForSlaveProcs[0] << endl);
			}

		//	master processing done!
		//	now find all processes which contain master interfaces to local slave interfaces
			vector<int> masterProcs;
			GetLayoutTargetProcs(masterProcs, newLayout->slave());
			
			UG_LOG("masterProcs: " << masterProcs.size() << endl);
			if(masterProcs.size())
				UG_LOG("masterProcs[0]: " << masterProcs[0] << endl);

			vector<int> recvSizes (masterProcs.size());
			BinaryBuffer recvBuf;

			newLayout->proc_comm().distribute_data(
					recvBuf, &recvSizes.front(), &masterProcs.front(),
					(int)masterProcs.size(),
					sendBuf.buffer(), &msgSizeForSlaveProcs.front(),
					&slaveProcs.front(), (int)slaveProcs.size(), 7553173);

			if(!recvSizes.empty()){
				UG_LOG("recvSizes: " << recvSizes[0] << endl);
			}

		//	Extract content of recvBuf and create new slave-overlap
			for(size_t i = 0; i < masterProcs.size(); ++i)
			{
				IndexLayout::Interface& itfc = newLayout->slave_overlap()
													.interface(masterProcs[i]);

				size_t oldReadPos = recvBuf.read_pos();
				while(recvBuf.read_pos() < oldReadPos + recvSizes[i]){
					AlgebraID globID;
					Deserialize(recvBuf, globID);
					size_t locID;
					if(m_algIDHash.get_entry(locID, globID)){
						itfc.push_back(locID);
					}
					else{
						UG_THROW("GlobalID " << globID << " expected on this process "
						         "but not found.");
					}
				}
			}

		//	set h-slave rows to dirichlet rows
			// IndexLayout& slaveLayout = newLayout->slave();
			// for(size_t lvl = 0; lvl < slaveLayout.num_levels(); ++lvl){
			// 	for(IndexLayout::const_iterator interfaceIter = slaveLayout.begin(lvl);
			// 		interfaceIter != slaveLayout.end(lvl); ++interfaceIter)
			// 	{
			// 		const IndexLayout::Interface& interface = slaveLayout.interface(interfaceIter);
			// 		for(IndexLayout::Interface::const_iterator iter = interface.begin();
			// 			iter != interface.end(); ++iter)
			// 		{
			// 			SetDirichletRow(m_mat, interface.get_element(iter));
			// 		}
			// 	}
			// }
		}

	private:
		TMatrix& m_mat;

	///	map localID->globalID
		AlgebraIDVec&		m_globalIDs;

	///	map globalID->localID
		AlgebraIDHashList	m_algIDHash;

		struct ExtCon {
			ExtCon (size_t _fromInd, const AlgebraID& _fromID,
			        const AlgebraID& _toID, int _conProc) :
				fromInd (_fromInd),
				fromID (_fromID),
				toID (_toID),
				conProc (_conProc)
			{}

			size_t		fromInd;
			AlgebraID	fromID;
			AlgebraID	toID;
			int			conProc;

			bool operator < (const ExtCon& ec) const
			{
				if(conProc < ec.conProc) return true;
			 	else if(conProc > ec.conProc) return false;
			 	if(toID < ec.toID) return true;
			 	else if(toID > ec.toID) return false;
			 	return (fromID < ec.fromID);
			}
		};

	///	New connections received from other processes.
		std::map<ExtCon, block_t>	m_recvNewCons;
	///	Here we just collect the new IDs which are added to the local matrix
		std::set<AlgebraID>			m_recvNewIDs;
};



template <class TMatrix>
void CreateOverlap (TMatrix& mat)
{
	using namespace std;
	PROFILE_FUNC_GROUP("algebra parallelization");

	vector<AlgebraID> globalIDs;

	GenerateGlobalAlgebraIDs(mat.layouts()->comm(),
	                         globalIDs,
	                         mat.num_rows(),
	                         mat.layouts()->master(),
	                         mat.layouts()->slave());

	ComPol_MatCreateOverlap<TMatrix> comPolOverlap(mat, globalIDs);
	mat.layouts()->comm().send_data(mat.layouts()->slave(), comPolOverlap);
	mat.layouts()->comm().receive_data(mat.layouts()->master(), comPolOverlap);
	mat.layouts()->comm().communicate();

	comPolOverlap.post_process();

//	todo: 	Remove this once Overlap creation is stable. Check for redistributed grids
//			which have h-masters and h-slaves on one process!
	TestHorizontalAlgebraLayouts(mat, NULL, false);

//	Make matrix partialy consistent on slave interfaces
	ComPol_MatCopyRowsOverlap0<TMatrix> comPolMatCopy(mat, globalIDs);
	mat.layouts()->comm().send_data(mat.layouts()->master(), comPolMatCopy);
	mat.layouts()->comm().receive_data(mat.layouts()->slave(), comPolMatCopy);
	mat.layouts()->comm().communicate();
}
	
}//	end of namespace

#endif	//__H__UG_matrix_overlap_impl
