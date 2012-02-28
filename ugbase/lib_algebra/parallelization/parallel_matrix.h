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
#include "lib_algebra/common/operations.h"

namespace ug
{

///\ingroup lib_algebra_parallelization

///\brief Wrapper for sequential matrices to handle them in parallel
/**
 * A ParallelMatrix is a wrapper around a sequential vector to make it usable
 * in parallel. It has all the function a sequential matrix supports, since it
 * is derived from it. Furthermore the ParallelStorageType is remembered.
 * Currently only additive storage type is implemented. In addition some
 * functions of the sequential matrix are overwritten to adapted the functionality
 * to parallel (e.g. set)
 *
 * Please Note:
 * The Implementation is not yet finished.
 *
 *\tparam		TMatrix		Sequential Matrix Type
 */
template <typename TMatrix>
class ParallelMatrix : public TMatrix
{
	public:
		enum {rows_sorted=TMatrix::rows_sorted};


	private:
	// 	disallow copy constructor
		ParallelMatrix(const ParallelMatrix&);

	public:
	///	own type
		typedef ParallelMatrix<TMatrix> this_type;

	public:
	///	Default Constructor
		ParallelMatrix()
		: TMatrix(), m_type(PST_UNDEFINED),
		  	m_pSlaveLayout(NULL), m_pMasterLayout(NULL),
			m_pCommunicator(NULL)
		{}

	///	Constructor setting the layouts
		ParallelMatrix(	IndexLayout& slaveLayout, IndexLayout masterLayout)
		: TMatrix(), m_type(PST_UNDEFINED),
			m_pCommunicator(NULL)
		{
			set_layouts(masterLayout, slaveLayout);
		}

		//////////////////////////////
		// Layouts and communicator
		//////////////////////////////

	///	sets the domain decomposition level to be used
		void set_layouts(IndexLayout& masterLayout, IndexLayout& slaveLayout)
		{
		//	set layout
			m_pMasterLayout = &masterLayout;
			m_pSlaveLayout = &slaveLayout;
		}

	///	sets slave layout
		void set_slave_layout(IndexLayout& slaveLayout)
		{
			m_pSlaveLayout = &slaveLayout;
		}

	///	sets slave layout
		void set_master_layout(IndexLayout& masterLayout)
		{
			m_pMasterLayout = &masterLayout;
		}

	///	returns the slave layout
		IndexLayout& slave_layout() const
		{
			UG_ASSERT(m_pSlaveLayout != NULL, "No Slave Layout set, but requested.");
			return *m_pSlaveLayout;
		}

	///	returns the master layout
		IndexLayout& master_layout() const
		{
			UG_ASSERT(m_pMasterLayout != NULL, "No Slave Layout set, but requested.");
			return *m_pMasterLayout;
		}

	///	sets a communicator
		void set_communicator(pcl::ParallelCommunicator<IndexLayout>& pc)
		{
			m_pCommunicator = &pc;
		}

	///	returns the communicator
		pcl::ParallelCommunicator<IndexLayout>&
		communicator()
		{
			UG_ASSERT(m_pCommunicator != NULL, "No communicator set, but requested.");
			return *m_pCommunicator;
		}

		const pcl::ParallelCommunicator<IndexLayout>&
		communicator() const
		{
			UG_ASSERT(m_pCommunicator != NULL, "No communicator set, but requested.");
			return *m_pCommunicator;
		}

	///	sets a process communicator
		void set_process_communicator(const pcl::ProcessCommunicator& pc)
		{
			m_processCommunicator = pc;
		}

	///	returns the process communicator
		pcl::ProcessCommunicator&
		process_communicator() {return m_processCommunicator;}

		const pcl::ProcessCommunicator&
		process_communicator() const {return m_processCommunicator;}

		/////////////////////////
		// Storage type handling
		/////////////////////////

	/// sets the storage type
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		void set_storage_type(uint type) {m_type = type;}

	/// adds a storage type
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		void add_storage_type(uint type) {m_type |= type;}

	/// removes a storage type
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		void remove_storage_type(uint type) {m_type &= ~type;}

	/// changes to the requested storage type if possible
		bool change_storage_type(ParallelStorageType type);

	/// returns if the current storage type has a given representation
	/**	type may be any or-combination of constants enumerated in ug::ParallelStorageType.*/
		bool has_storage_type(uint type) const
			{return type == PST_UNDEFINED ? m_type == PST_UNDEFINED : (m_type & type) == type;}

	/// returns storage type mask
		uint get_storage_mask() const { return m_type; }

	/// copies the storage type from another vector
		void copy_storage_type(const this_type& v) {m_type = v.m_type;}

		/////////////////////////
		// OverWritten functions
		/////////////////////////

	/// calculate res = A x
		template<typename TPVector>
		bool apply(TPVector &res, const TPVector &x) const;

	/// calculate res = A.T x
		template<typename TPVector>
		bool apply_transposed(TPVector &res, const TPVector &x) const;

	/// calculate res -= A x
		template<typename TPVector>
		bool matmul_minus(TPVector &res, const TPVector &x) const;

	///	copy layouts from another parallel matrix
		void copy_layouts(const this_type &v)
		{
			m_pSlaveLayout = v.m_pSlaveLayout;
			m_pMasterLayout = v.m_pMasterLayout;

			m_pCommunicator = v.m_pCommunicator;
			m_processCommunicator = v.m_processCommunicator;
		}

		///	assignment
		this_type &operator =(const this_type &M)
		{
		//	forward to sequential matrices
			TMatrix::operator= (*dynamic_cast<const TMatrix*>(&M));

		//	copy storage type and layouts
			copy_storage_type(M);
			copy_layouts(M);

		//	we're done
			return *this;
		}

	private:
	//  type of storage  (i.e. consistent, additiv, additiv unique)
		uint m_type;

	// 	index layout for slave dofs (0 is process-wise (finest grained) partition)
		IndexLayout* m_pSlaveLayout;

	// 	index layout for master dofs
		IndexLayout* m_pMasterLayout;

	// 	communicator for direct neighbor communication
		pcl::ParallelCommunicator<IndexLayout>* m_pCommunicator;

	// 	process communicator (world by default)
		pcl::ProcessCommunicator m_processCommunicator;
};

//	predaclaration.
//	this type may already be declared somewhere else, which shouldn't hurt.
template<typename T>
struct matrix_algebra_type_traits;

template<typename T>
struct matrix_algebra_type_traits<ParallelMatrix<T> >
{
	static const int type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

} // end namespace ug

#include "parallel_matrix_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__ */
