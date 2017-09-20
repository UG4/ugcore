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

/// Extends the matrix with connections to other processes
/** This policy is only intended to be used for slave->master communication.
 */
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
			: m_mat(rMat), m_vGlobalID(vGlobalID)
		{
			UG_ASSERT(vGlobalID.size() >= m_mat.num_rows(), "too few Global ids");
			m_tmpConProcs.resize(vGlobalID.size(), -1);

		//	we create a new layout, since the old one may be shared between many
		//	different vectors and matrices. H-Master and H-Slave layouts are still
		//	the same.
			m_newLayout = make_sp(new AlgebraLayouts);
			*m_newLayout = *m_mat.layouts();
			m_newLayout->enable_overlap(true);
			m_mat.set_layouts (m_newLayout);
		}

	///	writes the interface values into a buffer that will be sent
		virtual bool collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_collect, "algebra parallelization");
			typedef typename TMatrix::row_iterator row_iterator;

		//	mark all entries in this interface
			const int targetProc = interface.get_target_proc ();

			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				m_tmpConProcs [interface.get_element(iter)] = targetProc;
			}


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
					Serialize(buff, m_vGlobalID[k]);

				//	write entry into buffer
					Serialize(buff, a_ik);

				//	store as sent new connection if the connection did not yet
				//	exist on the target proc
					if (m_tmpConProcs [k] != targetProc) {
						m_sentNewIDs.insert(SentID(k, m_vGlobalID[k], targetProc));
					}
				}
			}

		///	done
			return true;
		}

		virtual bool
		end_layout_collection(const Layout* pLayout)	
		{
			using namespace std;

		//	create new slave-overlap
			{
				IndexLayout::Interface* itfc = NULL;

				for (typename set<SentID>::iterator iter = m_sentNewIDs.begin();
				     iter != m_sentNewIDs.end(); ++iter)
				{
					const SentID& curSentID = *iter;
					const int targetProc = curSentID.conProc;
					if (!itfc || targetProc != itfc->get_target_proc())
						itfc = &m_newLayout->slave_overlap().interface(targetProc);

					itfc->push_back(curSentID.ind);
				}
			}

		//	set h-slave rows to dirichlet rows
			IndexLayout& slaveLayout = m_newLayout->slave();
			for(size_t lvl = 0; lvl < slaveLayout.num_levels(); ++lvl){
				for(IndexLayout::const_iterator interfaceIter = slaveLayout.begin(lvl);
					interfaceIter != slaveLayout.end(lvl); ++interfaceIter)
				{
					const IndexLayout::Interface& interface = slaveLayout.interface(interfaceIter);
					for(IndexLayout::Interface::const_iterator iter = interface.begin();
						iter != interface.end(); ++iter)
					{
						SetDirichletRow(m_mat, interface.get_element(iter));
					}
				}
			}

			return true;
		}



		virtual bool
		begin_layout_extraction(const Layout* pLayout)
		{
		//	fill the map global->local
			GenerateAlgebraIDHashList(m_algIDHash, m_vGlobalID);
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
						const ExtCon ec(index, m_vGlobalID[index], gID, targetProc);
						m_recvNewCons [ec] = block;
						m_recvNewIDs.insert (gID);
					}
				}
			}

		///	done
			return true;
		}

		virtual bool
		end_layout_extraction(const Layout* pLayout)
		{
		//	create new master overlap
			using namespace std;

			const size_t oldSize = m_mat.num_rows();
			const size_t numNewInds = m_recvNewIDs.size();
			const size_t newSize = oldSize + numNewInds;

			if(newSize == oldSize)
				return true;

		//	add new entries to the algebra hash
			{
				size_t i = oldSize;
				for(set<AlgebraID>::iterator iter = m_recvNewIDs.begin();
				    iter != m_recvNewIDs.end(); ++iter, ++i)
				{
					m_algIDHash.insert(*iter, i);
				}
			}


		//	For each new DoF we set a DirichletRow
			m_mat.resize_and_keep_values (newSize, newSize);
			for(size_t i = oldSize; i < newSize; ++i){
				SetDirichletRow(m_mat, i);
			}

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
						itfc = &m_newLayout->master_overlap().interface(targetProc);
						added.clear();
						added.resize(numNewInds, false);
					}

					size_t toInd;
					UG_COND_THROW(!m_algIDHash.get_entry(toInd, curExtCon.toID),
					              "Expected AlgebraID not found in Hash");

					if(!added[toInd - oldSize])
						itfc->push_back(toInd);
					added[toInd - oldSize] = true;

					m_mat(curExtCon.fromInd, toInd) += iter->second;
				}
			}

			return true;
		}

	private:
		TMatrix& m_mat;

	///	map localID->globalID
		AlgebraIDVec&		m_vGlobalID;

	///	map globalID->localID
		AlgebraIDHashList	m_algIDHash;

	///	helper array to identify whether a coefficient is contained in an interface to a specific process
		std::vector<int>	m_tmpConProcs;

	///	In this fresh layout we'll combine old and new interfaces
		SmartPtr <AlgebraLayouts> m_newLayout;


		struct SentID {
			SentID (size_t _ind, const AlgebraID& _ID, int _conProc) :
				ind (_ind),
				ID (_ID),
				conProc (_conProc)
			{}

			size_t		ind;
			AlgebraID	ID;
			int			conProc;

			bool operator < (const SentID& ec) const
			{
				if(conProc < ec.conProc) return true;
			 	else if(conProc > ec.conProc) return false;
			 	return (ID < ec.ID);
			}
		};

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

	///	Internal DoF indices to whom connections were sent to other processes.
		std::set<SentID>	m_sentNewIDs;
	///	New connections received from other processes.
		std::map<ExtCon, block_t>	m_recvNewCons;
	///	Here we just collect the new IDs which are added to the local matrix
		std::set<AlgebraID>	m_recvNewIDs;
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

//	todo: 	Remove this once Overlap creation is stable. Check for redistributed grids
//			which have h-masters and h-slaves on one process!
	TestHorizontalAlgebraLayouts(mat);
}
	
}//	end of namespace

#endif	//__H__UG_matrix_overlap_impl
