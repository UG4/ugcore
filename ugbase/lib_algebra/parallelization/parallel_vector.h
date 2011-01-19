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
 * overwritten to adapted the functionality to parallel (e.g. two_norm, set).
 *
 * \tparam 	TVector		Sequential Vector type
 */
template <typename TVector>
class ParallelVector : public TVector
{
	private:
	// 	disallow copy constructor
		ParallelVector(const ParallelVector&);

	public:
	///	own type
		typedef ParallelVector<TVector> this_type;

	public:
	///	Default constructor
		ParallelVector()
		: TVector(), m_type(PST_UNDEFINED), m_ddlev(0),
			m_pVerticalSlaveLayout(NULL), m_pVerticalMasterLayout(NULL)
		{
			m_vpSlaveLayout.clear();
			m_vpMasterLayout.clear();
			m_vCommunicator.clear();
		}

	///	Constructor setting the Layouts
		ParallelVector(	IndexLayout& slaveLayout, IndexLayout masterLayout,
						IndexLayout& verticalSlaveLayout,
						IndexLayout& verticalMasterLayout)
		: TVector(), m_type(PST_UNDEFINED), m_ddlev(0),
			m_pVerticalSlaveLayout(&verticalSlaveLayout),
			m_pVerticalMasterLayout(&verticalMasterLayout)
		{
			set_slave_layout(slaveLayout);
			set_master_layout(masterLayout);
			m_vCommunicator.clear();
		}

		/////////////////////////////////
		// Layouts and communicator
		/////////////////////////////////

	///	returns the number of domain decomposition level
		size_t num_domain_decomposition_level() const
		{
		//	check that level exists
			if(m_vpSlaveLayout.size() != m_vpMasterLayout.size())
			{
				UG_LOG("ERROR in 'set_domain_decomposition_level':"
						" Different numbers of Master and Slave Layouts.\n");
				throw(UGFatalError("Cannot determine number of Domain Decomp"));
			}
			return m_vpSlaveLayout.size();
		}

	///	sets the domain decomposition level to be used
		bool set_domain_decomposition_level(size_t ddlev)
		{
		//	check that level exists
			if(!(ddlev < num_domain_decomposition_level()))
			{
				UG_LOG("ERROR in 'set_domain_decomposition_level':"
						" accessing level, that does not exist.\n");
				return false;
			}

		//	check that layouts are really set
			if(m_vpSlaveLayout.at(ddlev) == NULL || m_vpMasterLayout.at(ddlev) == NULL)
			{
				UG_LOG("ERROR in 'set_domain_decomposition_level':"
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

	///	sets the vertical slave layout
		void set_vertical_slave_layout(IndexLayout& layout)
			{m_pVerticalSlaveLayout = &layout;}

	///	sets the vertical master layout
		void set_vertical_master_layout(IndexLayout& layout)
			{m_pVerticalMasterLayout = &layout;}

	///	returns the slave layout
		IndexLayout& get_slave_layout(size_t ddlev = 0)	{return *(m_vpSlaveLayout.at(ddlev));}

	///	returns the master layout
		IndexLayout& get_master_layout(size_t ddlev = 0) {return *(m_vpMasterLayout.at(ddlev));}

	///	returns the vertical slave layout
		IndexLayout& get_vertical_slave_layout()	{return *m_pVerticalSlaveLayout;}

	///	returns the vertical slave layout
		IndexLayout& get_vertical_master_layout() {return *m_pVerticalMasterLayout;}

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

		//////////////////////////////////////////////////
		// overwritten functions of sequential vector
		//////////////////////////////////////////////////

	/// two norm (overwrites TVector::two_norm())
	/**
	 * Returns the two norm of the Vector. First, the two norm of each process
	 * is computed using the two_norm() method of the sequential vector. Then,
	 * the norms are summarized over all processes.
	 */
		// todo: should be const
		number two_norm();

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
		bool set_random(number from, number to, ParallelStorageType type=PST_CONSISTENT);

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
			ParallelStorageType mask = get_storage_mask() & v.get_storage_mask();

		//	check mask
			UG_ASSERT(mask != 0, "cannot substract vector v");
			if(mask == 0)
				throw(UG_ERROR_IncompatibleParallelStorageType(
						get_storage_mask(), v.get_storage_mask()));

		//	set this vector to mask
			set_storage_type(mask);

		//	forward
			TVector::operator-=(*dynamic_cast<const TVector*>(&v));

		//	we're done
			return *this;
		}

	///	add a vector
		this_type &operator +=(const this_type &v)
		{
		//	compute parallel storage mask
			ParallelStorageType mask = get_storage_mask() & v.get_storage_mask();

		//	check mask
			UG_ASSERT(mask != 0, "cannot add vector v");
			if(mask == 0)
				throw(UG_ERROR_IncompatibleParallelStorageType(
						get_storage_mask(), v.get_storage_mask()));

		//	set to new mask
			set_storage_type(mask);

		// 	forward to sequential vector
			TVector::operator+=(*dynamic_cast<const TVector*>(&v));

		//	we're done
			return *this;
		}

	protected:
	///	copy layouts from another parallel vector
		void copy_layouts(const this_type &v)
		{
			m_vpSlaveLayout = v.m_vpSlaveLayout;
			m_vpMasterLayout = v.m_vpMasterLayout;
			m_pVerticalSlaveLayout = v.m_pVerticalSlaveLayout;
			m_pVerticalMasterLayout = v.m_pVerticalMasterLayout;

			m_vCommunicator = v.m_vCommunicator;
			m_vProcessCommunicator = v.m_vProcessCommunicator;
		}

	private:
	// 	type of storage  (i.e. consistent, additiv, additiv unique)
		int m_type;

	//	current domain decomposition level
		int m_ddlev;

	// 	index layout for slave dofs (0 is process-wise (finest grained) partition)
		std::vector<IndexLayout*> m_vpSlaveLayout;

	// 	index layout for master dofs
		std::vector<IndexLayout*> m_vpMasterLayout;

	// 	index layout for vertical slave dofs
		IndexLayout* m_pVerticalSlaveLayout;

	// 	index layout for vertical master dofs
		IndexLayout* m_pVerticalMasterLayout;

	// 	communicator for direct neighbor communication
		std::vector<pcl::ParallelCommunicator<IndexLayout>*> m_vCommunicator;

	// 	process communicator (world by default)
		std::vector<pcl::ProcessCommunicator> m_vProcessCommunicator;
};

} // end namespace ug

#include "parallel_vector_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__ */
