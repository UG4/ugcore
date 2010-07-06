/*
 * parallel_vector_impl.h
 *
 *  Created on: 3.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__

#include "parallel_vector.h"

namespace ug
{

template <typename TVector>
bool
ParallelVector<TVector>::
change_storage_type(ParallelStorageType type)
{
	// can only change if current state is defined
	if(has_storage_type(PST_UNDEFINED)) return false;

	// if already in that type
	if(has_storage_type(type)) return true;

	// else switch to that type
	switch(type)
	{
	case PST_CONSISTENT:
			 if(has_storage_type(PST_UNIQUE)){
				if(m_pMasterLayout == NULL || m_pSlaveLayout == NULL) return false;
				 UniqueToConsistent(this, *m_pMasterLayout, *m_pSlaveLayout, &m_Communicator);
				 set_storage_type(PST_CONSISTENT);
				 break;
			 }
			 else if(has_storage_type(PST_ADDITIVE)){
				if(m_pMasterLayout == NULL || m_pSlaveLayout == NULL) return false;
				AdditiveToConsistent(this, *m_pMasterLayout, *m_pSlaveLayout, &m_Communicator);
				set_storage_type(PST_CONSISTENT);
				break;
			}
			else return false;
	case PST_ADDITIVE:
			if(has_storage_type(PST_UNIQUE)){
				add_storage_type(PST_ADDITIVE);
				break;
			}
			else if(has_storage_type(PST_CONSISTENT)){
				if(m_pMasterLayout == NULL) return false;
				ConsistentToUnique(this, *m_pSlaveLayout);
				set_storage_type(PST_ADDITIVE);
				add_storage_type(PST_UNIQUE);
				break;
			}
			else return false;
	case PST_UNIQUE:
			if(has_storage_type(PST_ADDITIVE)){
				if(m_pMasterLayout == NULL || m_pSlaveLayout == NULL) return false;
				AdditiveToUnique(this, *m_pMasterLayout, *m_pSlaveLayout, &m_Communicator);
				add_storage_type(PST_UNIQUE);
				break;
			}
			else if(has_storage_type(PST_CONSISTENT)){
				if(m_pSlaveLayout == NULL) return false;
				ConsistentToUnique(this, *m_pSlaveLayout);
				set_storage_type(PST_ADDITIVE);
				add_storage_type(PST_UNIQUE);
				break;
			}
			else return false;
	default: return false;
	}
	return true;
}

template <typename TVector>
bool
ParallelVector<TVector>::
set(number w, ParallelStorageType type)
{
	if(type == PST_UNDEFINED) return false;

	// set all local vector to value. Therefore parallel vector is consistent
	if(!TVector::set(w)) return false;
	set_storage_type(PST_CONSISTENT);

	// if w == 0.0, we have all types
	if(w == 0.0){
		add_storage_type(PST_ADDITIVE);
		add_storage_type(PST_UNIQUE);
		return true;
	}

	// consistent required
	if(type & PST_CONSISTENT) return true;

	// additive or additive unique
	return change_storage_type(PST_UNIQUE);
}

template <typename TVector>
inline
number
ParallelVector<TVector>::
two_norm()
{
	// step 1: make vector d additive unique
	if(!change_storage_type(PST_UNIQUE)) {
		UG_ASSERT(0, "Cannot change to unique representation.");
		return -1;
	}

	// step 2: compute new defect norm
	double tNormLocal = (double)TVector::two_norm();
	tNormLocal *= tNormLocal;

	// step 3: sum local norms
	double tNormGlobal;
	// TODO: Replace this by allreduce only between layout
	pcl::AllReduce(&tNormLocal, &tNormGlobal, 1, PCL_DT_DOUBLE, PCL_RO_SUM);

	// step 4: return global norm
	return sqrt((number)tNormGlobal);
}

template <typename TVector>
inline
number
ParallelVector<TVector>::
dotprod(const this_type& v)
{
	// step 0: check that storage type is given
	if(this->has_storage_type(PST_UNDEFINED) || v.has_storage_type(PST_UNDEFINED))
	{
		// TODO: rethink this arrangement, should throw an error
		UG_ASSERT(0, "Parallel storage type not given.\n");
	}

	bool check = false;
	// step 1: Check if good storage type are given (no communication needed)
	//         - additive (unique) <-> consistent is ok
	//         - unique <-> unique is ok
	if(this->has_storage_type(PST_ADDITIVE) && v.has_storage_type(PST_CONSISTENT)) check = true;
	if(this->has_storage_type(PST_CONSISTENT) && v.has_storage_type(PST_ADDITIVE)) check = true;
	if(this->has_storage_type(PST_UNIQUE)   && v.has_storage_type(PST_UNIQUE))     check = true;

	// step 2: fall back
	//         if storage type not as in the upper cases, communicate to correct solution
	//         a user of this function should ideally avoid such a change and do it outside of this function
	if(!check)
	{
		// unique <-> additive => consistent <-> additive
		if(this->has_storage_type(PST_UNIQUE) && v.has_storage_type(PST_ADDITIVE))
			{this->change_storage_type(PST_CONSISTENT);}
		// additive <-> unique => unique <-> unique
		else if(this->has_storage_type(PST_ADDITIVE) && v.has_storage_type(PST_UNIQUE))
			{this->change_storage_type(PST_UNIQUE);}
		// consistent <-> consistent => unique <-> consistent
		else {this->change_storage_type(PST_UNIQUE);}
	}

	// step 3: compute local dot product
	double tSumLocal = (double)TVector::dotprod(v);
	double tSumGlobal;

	// step 4: sum global contributions
	// TODO: Replace this by allreduce only between layout
	pcl::AllReduce(&tSumLocal, &tSumGlobal, 1, PCL_DT_DOUBLE, PCL_RO_SUM);

	// step 5: return result
	return tSumGlobal;
}


} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__ */
