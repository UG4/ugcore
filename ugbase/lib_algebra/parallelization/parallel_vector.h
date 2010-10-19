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
#include "../martin_algebra/template_operations/template_expressions.h"
#include "parallel_storage_type.h"

namespace ug
{

/**
 * A ParallelVector is a wrapper around a sequential vector to make it usable in parallel.
 * It has all the function a sequential vector supports, since it is publically derived from it.
 * Furthermore the ParallelStorageType is remembered and can be switched.
 * In addition some functions of the sequential vector are overwritten to adapted the functionality
 * to parallel (e.g. two_norm, set)
 */
template <typename TVector>
class ParallelVector : public TVector, public TE_VEC<ParallelVector<TVector> >
{
	private:
	// 	disallow copy constructor
		ParallelVector(const ParallelVector&);

	public:
		typedef ParallelVector<TVector> this_type;

	public:
		ParallelVector()
		: TVector(), m_type(PST_UNDEFINED),
			m_pSlaveLayout(NULL), m_pMasterLayout(NULL),
			m_pVerticalSlaveLayout(NULL), m_pVerticalMasterLayout(NULL), m_pCommunicator(NULL)
		{}

		ParallelVector(	IndexLayout& slaveLayout, IndexLayout masterLayout,
						IndexLayout& verticalSlaveLayout, IndexLayout& verticalMasterLayout)
		: TVector(), m_type(PST_UNDEFINED),
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


		/////////////////////////
		// overwritten functions of sequential vector
		/////////////////////////

		// two norm (overwrites TVector::two_norm())
		// TODO: should be const
		inline number two_norm();

		// dotprod (overwrites TVector::dotprod())
		inline number dotprod(const this_type& v);

		// set all dofs to value 'w' (overwrites TVector::set(number w))
		bool set(number w, ParallelStorageType type = PST_CONSISTENT);

		/////////////////////////
		// forward hidden (due to overloading in base class) functions of sequential vector
		/////////////////////////
//		bool set(const local_vector_type &u, const local_index_type &ind)
//			{return TVector::set(u, ind);}
		// TODO: Ask Martin why return is void, compile error in other operator=
/*		void operator= (const this_type& v)
			{return TVector::operator=(v);}
*/		// TODO: Ask Martin why double return
		//double operator = (double d)
		//	{return TVector::operator=(d);}

/*		this_type &operator =(const this_type &v)
		{
			copy_storage_type(v);

			(*(TVector*)this) = (TVector)v;
		}
*/
		this_type &operator -=(const this_type &v)
		{
			ParallelStorageType mask = get_storage_mask() & v.get_storage_mask();
			UG_ASSERT(mask != 0, "cannot substract vector v");
			if(mask == 0) throw(UG_ERROR_IncompatibleParallelStorageType(get_storage_mask(), v.get_storage_mask()));
			set_storage_type(mask);

			TVector::operator-=(*dynamic_cast<const TVector*>(&v));
			return *this;
		}
		this_type &operator +=(const this_type &v)
		{
			ParallelStorageType mask = get_storage_mask() & v.get_storage_mask();
			UG_ASSERT(mask != 0, "cannot add vector v");
			if(mask == 0) throw(UG_ERROR_IncompatibleParallelStorageType(get_storage_mask(), v.get_storage_mask()));
			set_storage_type(mask);

			TVector::operator+=(*dynamic_cast<const TVector*>(&v));
			return *this;
		}

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

#include "parallel_vector_impl.h"

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR__ */
