/*
 * parallel_matrix.h
 *
 *  Created on: 19.10.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__

#include "pcl/pcl.h"
#include "parallel_index_layout.h"
#include "parallelization_util.h"
#include "parallel_storage_type.h"

namespace ug
{

/**
 * A ParallelMatrix is a wrapper around a sequential vector to make it usable in parallel.
 * It has all the function a sequential matrix supports, since it is publically derived from it.
 * Furthermore the ParallelStorageType is remembered and can be switched.
 * In addition some functions of the sequential matrix are overwritten to adapted the functionality
 * to parallel (e.g. set)
 */
template <typename TMatrix>
class ParallelMatrix : public TMatrix
{
	private:
	// 	disallow copy constructor
	ParallelMatrix(const ParallelMatrix&);

	public:
		typedef ParallelMatrix<TMatrix> this_type;

	public:
		ParallelMatrix()
		: TMatrix(), m_type(PST_UNDEFINED),
			m_pSlaveLayout(NULL), m_pMasterLayout(NULL),
			m_pVerticalSlaveLayout(NULL), m_pVerticalMasterLayout(NULL), m_pCommunicator(NULL)
		{}

		ParallelMatrix(	IndexLayout& slaveLayout, IndexLayout masterLayout,
						IndexLayout& verticalSlaveLayout, IndexLayout& verticalMasterLayout)
		: TMatrix(), m_type(PST_UNDEFINED),
			m_pSlaveLayout(&slaveLayout), m_pMasterLayout(&masterLayout),
			m_pVerticalSlaveLayout(&verticalSlaveLayout), m_pVerticalMasterLayout(&verticalMasterLayout), m_pCommunicator(NULL)
		{}

		/////////////////////////
		// Layouts and communicator
		/////////////////////////

		inline void set_slave_layout(IndexLayout& layout)	{m_pSlaveLayout = &layout;}
		inline void set_master_layout(IndexLayout& layout)	{m_pMasterLayout = &layout;}
		inline void set_vertical_slave_layout(IndexLayout& layout)	{m_pVerticalSlaveLayout = &layout;}
		inline void set_vertical_master_layout(IndexLayout& layout) {m_pVerticalMasterLayout = &layout;}

		inline IndexLayout& get_slave_layout()	{return *m_pSlaveLayout;}
		inline IndexLayout& get_master_layout()	{return *m_pMasterLayout;}
		inline IndexLayout& get_vertical_slave_layout()	{return *m_pVerticalSlaveLayout;}
		inline IndexLayout& get_vertical_master_layout() {return *m_pVerticalMasterLayout;}

		inline void set_communicator(pcl::ParallelCommunicator<IndexLayout>& pc) {m_pCommunicator = &pc;}
		inline pcl::ParallelCommunicator<IndexLayout>& get_communicator() {return *m_pCommunicator;}

		inline void set_process_communicator(const pcl::ProcessCommunicator& pc)	{m_processCommunicator = pc;}
		inline pcl::ProcessCommunicator& get_process_communicator()					{return m_processCommunicator;}

		/////////////////////////
		// Storage type handling
		/////////////////////////

		// sets the storage type
		void set_storage_type(ParallelStorageType type) {m_type = type;}

		// adds a storage type
		void add_storage_type(ParallelStorageType type) {m_type |= type;}

		// removes a storage type
		void remove_storage_type(ParallelStorageType type) {m_type &= ~type;}

		// changes to the requested storage type if possible
		bool change_storage_type(ParallelStorageType type);

		// returns if the current storage type has a given representation
		bool has_storage_type(ParallelStorageType type) const {return (bool)(m_type & type);}

		// returns storage type mask
		ParallelStorageType get_storage_mask() const { return (ParallelStorageType) m_type; }

		// copies the storage type from another vector
		void copy_storage_type(const this_type& v) {m_type = v.m_type;}


	private:
		// type of storage  (i.e. consistent, additiv, additiv unique)
		int m_type;

		// index layout for slave dofs
		IndexLayout* m_pSlaveLayout;

		// index layout for master dofs
		IndexLayout* m_pMasterLayout;

		// index layout for vertical slave dofs
		IndexLayout* m_pVerticalSlaveLayout;

		// index layout for vertical master dofs
		IndexLayout* m_pVerticalMasterLayout;

		// communicator for direct neighbor communication
		pcl::ParallelCommunicator<IndexLayout>* m_pCommunicator;

		// process communicator (world by default)
		pcl::ProcessCommunicator m_processCommunicator;
};

} // end namespace ug

#include "parallel_matrix_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__ */
