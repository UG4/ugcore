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
		: TMatrix(), m_type(PST_UNDEFINED), m_ddlev(0)
		{
			m_vpSlaveLayout.clear();
			m_vpMasterLayout.clear();
			m_vCommunicator.clear();
			m_vProcessCommunicator.clear();
		}

	///	Constructor setting the layouts
		ParallelMatrix(	IndexLayout& slaveLayout, IndexLayout masterLayout)
		: TMatrix(), m_type(PST_UNDEFINED), m_ddlev(0)
		{
			set_slave_layout(slaveLayout);
			set_master_layout(masterLayout);
			m_vCommunicator.clear();
			m_vProcessCommunicator.clear();
		}

		//////////////////////////////
		// Layouts and communicator
		//////////////////////////////

	///	returns the number of domain decomposition level
		size_t num_layouts() const
		{
		//	check that level exists
			if(m_vpSlaveLayout.size() != m_vpMasterLayout.size())
			{
				UG_LOG("ERROR in 'use_layout':"
						" Different numbers of Master and Slave Layouts.\n");
				throw(UGFatalError("Cannot determine number of Domain Decomp"));
			}
			return m_vpSlaveLayout.size();
		}

	///	sets the domain decomposition level to be used
		bool use_layout(size_t ddlev)
		{
		//	check that level exists
			if(!(ddlev < num_layouts()))
			{
				UG_LOG("ERROR in 'use_layout':"
						" accessing level, that does not exist.\n");
				return false;
			}

		//	check that layouts are really set
			if(m_vpSlaveLayout.at(ddlev) == NULL || m_vpMasterLayout.at(ddlev) == NULL)
			{
				UG_LOG("ERROR in 'use_layout':"
						" Although level exist, layouts are not set.\n");
				return false;
			}

		//	set new level
			m_ddlev = ddlev;

		//	we're done
			return true;
		}

	///	sets the slave layout
		void set_slave_layout(IndexLayout& layout, size_t ddlev = 0)
		{
		//	resize if needed
			if(m_vpSlaveLayout.size() <= ddlev)
				m_vpSlaveLayout.resize(ddlev+1, NULL);

		//	remember layout
			m_vpSlaveLayout[ddlev] = &layout;
		}

	///	sets the slave layouts
		void set_slave_layouts(std::vector<IndexLayout>& layouts)
		{
		//	clear old
			m_vpSlaveLayout.clear();

		//	set new
			for(size_t i = 0; i < layouts.size(); ++i)
			{
				m_vpSlaveLayout.push_back(&layouts[i]);
			}
		}

	///	sets the master layout
		void set_master_layout(IndexLayout& layout, size_t ddlev = 0)
		{
		//	resize if needed
			if(m_vpMasterLayout.size() <= ddlev)
				m_vpMasterLayout.resize(ddlev+1, NULL);

		//	remember layout
			m_vpMasterLayout[ddlev] = &layout;
		}

	///	sets the master layouts
		void set_master_layouts(std::vector<IndexLayout>& layouts)
		{
		//	clear old
			m_vpMasterLayout.clear();

		//	set new
			for(size_t i = 0; i < layouts.size(); ++i)
			{
				m_vpMasterLayout.push_back(&layouts[i]);
			}
		}

	///	returns the slave layout
		IndexLayout& get_slave_layout(size_t ddlev = 0)	{return *(m_vpSlaveLayout.at(ddlev));}

	///	returns the master layout
		IndexLayout& get_master_layout(size_t ddlev = 0) {return *(m_vpMasterLayout.at(ddlev));}

	///	sets a communicator
		void set_communicator(pcl::ParallelCommunicator<IndexLayout>& pc, size_t ddlev = 0)
		{
		//	resize if needed
			if(m_vCommunicator.size() <= ddlev)
				m_vCommunicator.resize(ddlev+1, NULL);

		//	remember layout
			m_vCommunicator[ddlev] = &pc;
		}

	///	sets all communicators
		void set_communicators(std::vector<pcl::ParallelCommunicator<IndexLayout> >& pc)
		{
		//	clear old
			m_vCommunicator.clear();

		//	set new
			for(size_t i = 0; i < pc.size(); ++i)
			{
				m_vCommunicator.push_back(&pc[i]);
			}
		}

	///	returns the communicator
		pcl::ParallelCommunicator<IndexLayout>&
		get_communicator(size_t ddlev = 0) {return *(m_vCommunicator.at(ddlev));}

	///	sets a process communicator
		void set_process_communicator(const pcl::ProcessCommunicator& pc, size_t ddlev = 0)
		{
		//	resize if needed
		if(m_vProcessCommunicator.size() <= ddlev)
			m_vProcessCommunicator.resize(ddlev+1);

		//	remember layout
			m_vProcessCommunicator[ddlev] = pc;
		}

	///	sets all process communicators
		void set_process_communicators(std::vector<pcl::ProcessCommunicator>& pc)
		{
		//	clear old
			m_vProcessCommunicator.clear();

		//	set new
			m_vProcessCommunicator = pc;
		}

	///	returns the process communicator
		pcl::ProcessCommunicator&
		get_process_communicator(size_t ddlev = 0) {return m_vProcessCommunicator.at(ddlev);}

		/////////////////////////
		// Storage type handling
		/////////////////////////

	/// sets the storage type
		void set_storage_type(ParallelStorageType type) {m_type = type;}

	/// adds a storage type
		void add_storage_type(ParallelStorageType type) {m_type |= type;}

	/// removes a storage type
		void remove_storage_type(ParallelStorageType type) {m_type &= ~type;}

	/// changes to the requested storage type if possible
		bool change_storage_type(ParallelStorageType type);

	/// returns if the current storage type has a given representation
		bool has_storage_type(ParallelStorageType type) const {return (bool)(m_type & type);}

	/// returns storage type mask
		ParallelStorageType get_storage_mask() const { return (ParallelStorageType) m_type; }

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


	private:
	//  type of storage  (i.e. consistent, additiv, additiv unique)
		int m_type;

	//	current domain decomposition level
		int m_ddlev;

	//  index layout for slave dofs
		std::vector<IndexLayout*> m_vpSlaveLayout;

	//  index layout for master dofs
		std::vector<IndexLayout*> m_vpMasterLayout;

	//  communicator for direct neighbor communication
		std::vector<pcl::ParallelCommunicator<IndexLayout>*> m_vCommunicator;

	//  process communicator (world by default)
		std::vector<pcl::ProcessCommunicator> m_vProcessCommunicator;
};

//	predaclaration.
//	this type may already be declared somewhere else, which shouldn't hurt.
template<typename T>
struct matrix_algebra_type_traits;

template<typename T>
struct matrix_algebra_type_traits<ParallelMatrix<T> >
{
	static const int type = matrix_algebra_type_traits<T>::type;
};

} // end namespace ug

#include "parallel_matrix_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX__ */
