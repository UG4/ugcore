/*
 * parallel_vector.h
 *
 *  Created on: 3.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__

#include "pcl/pcl.h"
#include "parallel_index_layout.h"
#include "parallelization_util.h"
#include "parallel_storage_type.h"
#include "common/assert.h"

namespace ug
{

///\ingroup lib_algebra_parallelization

/// Parallelization Wrapper to make a Vector usable in parallel
/**
 * A ParallelVector is a wrapper around a sequential vector to make it usable
 * in parallel. It has all the function a sequential vector supports, since it
 * is derived from it. Furthermore the ParallelStorageType is remembered and
 * can be switched. In addition some functions of the sequential vector are
 * overwritten to adapted the functionality to parallel (e.g. norm, set).
 *
 * \tparam 	TVector		Sequential Vector type
 */
template <typename TVector>
class ParallelVector : public TVector
{
	public:
		typedef typename TVector::value_type value_type;
		typedef typename TVector::vector_type vector_type;

	private:
	// 	disallow copy constructor
		//ParallelVector(const ParallelVector&);

	public:
	///	own type
		typedef ParallelVector<TVector> this_type;

	public:
	///	Default constructor
		ParallelVector()
		: TVector(), m_type(PST_UNDEFINED),
		  	m_pSlaveLayout(NULL), m_pMasterLayout(NULL),
			m_pVerticalSlaveLayout(NULL), m_pVerticalMasterLayout(NULL),
			m_pCommunicator(NULL)
		{}

	/// Resizing constructor
		ParallelVector(size_t length)
		: TVector(length), m_type(PST_UNDEFINED),
			m_pSlaveLayout(NULL), m_pMasterLayout(NULL),
			m_pVerticalSlaveLayout(NULL), m_pVerticalMasterLayout(NULL),
			m_pCommunicator(NULL)
		{}

	///	Constructor setting the Layouts
		ParallelVector(	IndexLayout& slaveLayout, IndexLayout masterLayout,
						IndexLayout& verticalSlaveLayout,
						IndexLayout& verticalMasterLayout)
		: TVector(), m_type(PST_UNDEFINED),
		  	m_pSlaveLayout(NULL), m_pMasterLayout(NULL),
			m_pVerticalSlaveLayout(&verticalSlaveLayout),
			m_pVerticalMasterLayout(&verticalMasterLayout),
			m_pCommunicator(NULL)
		{
			set_layouts(masterLayout, slaveLayout);
		}

		/////////////////////////////////
		// Layouts and communicator
		/////////////////////////////////

	///	sets layouts
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

	///	sets the vertical slave layout
		void set_vertical_layouts(IndexLayout& masterLayout, IndexLayout& slaveLayout)
		{
			m_pVerticalSlaveLayout = &slaveLayout;
			m_pVerticalMasterLayout = &masterLayout;
		}

	///	returns the slave layout
		IndexLayout& slave_layout() const
		{
			UG_ASSERT(m_pSlaveLayout != NULL,
			          "No Horizontal Slave Layout set, but requested.");
			return *m_pSlaveLayout;
		}

	///	returns the master layout
		IndexLayout& master_layout() const
		{
			UG_ASSERT(m_pMasterLayout != NULL,
			          "No Horizontal Master Layout set, but requested.");
			return *m_pMasterLayout;
		}

	///	returns the vertical slave layout
		IndexLayout& vertical_slave_layout() const
		{
			UG_ASSERT(m_pVerticalSlaveLayout != NULL,
			          "No Vertical Slave Layout set, but requested.");
			return *m_pVerticalSlaveLayout;
		}

	///	returns the vertical slave layout
		IndexLayout& vertical_master_layout() const
		{
			UG_ASSERT(m_pVerticalMasterLayout != NULL,
			          "No Vertical Master Layout set, but requested.");
			return *m_pVerticalMasterLayout;
		}

	///	sets a communicator
		void set_communicator(pcl::InterfaceCommunicator<IndexLayout>& pc)
		{
			m_pCommunicator = &pc;
		}

	///	returns the communicator
		pcl::InterfaceCommunicator<IndexLayout>& communicator() const
		{
			UG_ASSERT(m_pCommunicator != NULL,
			          "No Parallel Communicator set, but requested.");
			return *m_pCommunicator;
		}

	///	sets a process communicator
		void set_process_communicator(pcl::ProcessCommunicator& pc)
		{
			m_processCommunicator = pc;
		}

	///	returns the process communicator
		pcl::ProcessCommunicator&
		process_communicator() {return m_processCommunicator;}

	///	returns the process communicator
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

		//////////////////////////////////////////////////
		// overwritten functions of sequential vector
		//////////////////////////////////////////////////

	/// two norm (overwrites TVector::norm())
	/**
	 * Returns the two norm of the Vector. First, the two norm of each process
	 * is computed using the norm() method of the sequential vector. Then,
	 * the norms are summarized over all processes.
	 */
		number norm() const;

	/// dotprod (overwrites TVector::dotprod())
	/**
	 * Returns the dot product of the vector. First, the dot prod of each process
	 * is computed using the dotprod() method of the sequential vector. Then,
	 * the results are summarized over all processes.
	 */
		inline number dotprod(const this_type& v);

	/// set all entries to value and the storage type
		bool set(number w, ParallelStorageType type);

	///	sets all entries to a value and the storage type to consistent
		bool set(number w){return set(w, PST_CONSISTENT);}

	/// set all entries to a random number (overwrites TVector::set_random(number from, number to))
		bool set_random(number from, number to, ParallelStorageType type);
		bool set_random(number from, number to)
		{
			return set_random(from, to, PST_CONSISTENT);
		}

	///	assignment
		this_type &operator =(const this_type &v)
		{
		//	forward to sequential vectors
			TVector::operator=(*dynamic_cast<const TVector*>(&v));

		//	copy storage type and layouts
			copy_storage_type(v); copy_layouts(v);

		//	we're done
			return *this;
		}

	///	subtract a vector
		this_type &operator -=(const this_type &v)
		{
		//	compute storage mask
			uint mask = get_storage_mask() & v.get_storage_mask();

		//	check mask
			UG_ASSERT(mask != 0, "ERROR in 'ParallelVector::operator-=': Storage"
					" types not compatible. This: "<< get_storage_mask() <<
					" and other: " << v.get_storage_mask() << ".");
			if(mask == 0)
				throw(UG_ERROR_IncompatibleParallelStorageType(
						get_storage_mask(), v.get_storage_mask()));

		//	set this vector to mask
			m_type = mask;

		//	forward
			TVector::operator-=(*dynamic_cast<const TVector*>(&v));

		//	we're done
			return *this;
		}

	///	add a vector
		this_type &operator +=(const this_type &v)
		{
		//	compute parallel storage mask
			uint mask = get_storage_mask() & v.get_storage_mask();

		//	check mask
			UG_ASSERT(mask != 0, "ERROR in 'ParallelVector::operator+=': Storage"
					" types not compatible. This: "<< get_storage_mask() <<
					" and other: " << v.get_storage_mask() << ".");
			if(mask == 0)
				throw(UG_ERROR_IncompatibleParallelStorageType(
						get_storage_mask(), v.get_storage_mask()));

		//	set to new mask
			m_type = mask;

		// 	forward to sequential vector
			TVector::operator+=(*dynamic_cast<const TVector*>(&v));

		//	we're done
			return *this;
		}


	///	copy layouts from another parallel vector
		void copy_layouts(const this_type &v)
		{
			m_pSlaveLayout = v.m_pSlaveLayout;
			m_pMasterLayout = v.m_pMasterLayout;
			m_pVerticalSlaveLayout = v.m_pVerticalSlaveLayout;
			m_pVerticalMasterLayout = v.m_pVerticalMasterLayout;

			m_pCommunicator = v.m_pCommunicator;
			m_processCommunicator = v.m_processCommunicator;
		}
	protected:
	private:
	// 	type of storage  (i.e. consistent, additiv, additiv unique)
	//	holds or-combiation of constants enumerated in ug::ParallelStorageType.
		uint m_type;

	// 	index layout for slave dofs (0 is process-wise (finest grained) partition)
		IndexLayout* m_pSlaveLayout;

	// 	index layout for master dofs
		IndexLayout* m_pMasterLayout;

	// 	index layout for vertical slave dofs
		IndexLayout* m_pVerticalSlaveLayout;

	// 	index layout for vertical master dofs
		IndexLayout* m_pVerticalMasterLayout;

	// 	communicator for direct neighbor communication
		pcl::InterfaceCommunicator<IndexLayout>* m_pCommunicator;

	// 	process communicator (world by default)
		pcl::ProcessCommunicator m_processCommunicator;
};

} // end namespace ug

#include "parallel_vector_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__ */
