/*
 * parallel_vector_impl.h
 *
 *  Created on: 3.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__

#include <cmath>
#include "parallel_vector.h"
#include "common/common.h" // for UGError
// additions for profiling (14042011ih)
#include "common/profiler/profiler.h"
#define PROFILE_PARALLEL_VECTOR
#ifdef PROFILE_PARALLEL_VECTOR
	#define PARVEC_PROFILE_FUNC()		PROFILE_FUNC()
	#define PARVEC_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define PARVEC_PROFILE_END()		PROFILE_END()
#else
	#define PARVEC_PROFILE_FUNC()
	#define PARVEC_PROFILE_BEGIN(name)
	#define PARVEC_PROFILE_END()
#endif
// additions for profiling - end

namespace ug
{

template <typename TVector>
bool
ParallelVector<TVector>::
change_storage_type(ParallelStorageType type)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	choose communicator to use
	pcl::InterfaceCommunicator<IndexLayout>* pCommunicator;
	if(m_pCommunicator == NULL)
		pCommunicator = const_cast<pcl::InterfaceCommunicator<IndexLayout>*>(m_pCommunicator);
	else
		pCommunicator = NULL;

	// can only change if current state is defined
	if(has_storage_type(PST_UNDEFINED))
		UG_THROW("ParallelVector::change_storage_type: Trying to change"
				" storage type of a vector that has type PST_UNDEFINED.");

	// if already in that type
	if(has_storage_type(type)) return true;

	// else switch to that type
	switch(type)
	{
	case PST_CONSISTENT:
			 if(has_storage_type(PST_UNIQUE)){
			 	PARVEC_PROFILE_BEGIN(ParVec_CSTUnique2Consistent); // added 14042011ih
				if(m_pMasterLayout == NULL || m_pSlaveLayout == NULL)
					UG_THROW("ParallelVector::change_storage_type: No "
							"layouts given but trying to change type.");
				 UniqueToConsistent(this, *const_cast<IndexLayout*>(m_pMasterLayout),
				                    *const_cast<IndexLayout*>(m_pSlaveLayout),
				                    pCommunicator);
				 set_storage_type(PST_CONSISTENT);
				 PARVEC_PROFILE_END(); //ParVec_CSTUnique2Consistent // added 14042011ih
				 break;
			 }
			 else if(has_storage_type(PST_ADDITIVE)){
			 	PARVEC_PROFILE_BEGIN(ParVec_CSTAdditive2Consistent); // added 14042011ih
				if(m_pMasterLayout == NULL || m_pSlaveLayout == NULL)
					UG_THROW("ParallelVector::change_storage_type: No "
							"layouts given but trying to change type.");
				AdditiveToConsistent(this, *const_cast<IndexLayout*>(m_pMasterLayout),
				                     *const_cast<IndexLayout*>(m_pSlaveLayout),
				                     pCommunicator);
				set_storage_type(PST_CONSISTENT);
				PARVEC_PROFILE_END(); //ParVec_CSTAdditive2Consistent // added 14042011ih
				break;
			}
			else return false;
	case PST_ADDITIVE:
			if(has_storage_type(PST_UNIQUE)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTUnique2Additive); // added 14042011ih
				add_storage_type(PST_ADDITIVE);
				PARVEC_PROFILE_END(); //ParVec_CSTUnique2Additive // added 14042011ih
				break;
			}
			else if(has_storage_type(PST_CONSISTENT)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTConsistent2Additive); // added 14042011ih
				if(m_pMasterLayout == NULL)
					UG_THROW("ParallelVector::change_storage_type: No "
							"layouts given but trying to change type.");
				ConsistentToUnique(this, *const_cast<IndexLayout*>(m_pSlaveLayout));
				set_storage_type(PST_ADDITIVE);
				add_storage_type(PST_UNIQUE);
				PARVEC_PROFILE_END(); //ParVec_CSTConsistent2Additive // added 14042011ih
				break;
			}
			else return false;
	case PST_UNIQUE:
			if(has_storage_type(PST_ADDITIVE)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTAdditive2Unique); // added 14042011ih
				if(m_pMasterLayout == NULL || m_pSlaveLayout == NULL)
					UG_THROW("ParallelVector::change_storage_type: No "
							"layouts given but trying to change type.");
				AdditiveToUnique(this, *const_cast<IndexLayout*>(m_pMasterLayout),
				                 *const_cast<IndexLayout*>(m_pSlaveLayout),
				                 pCommunicator);
				add_storage_type(PST_UNIQUE);
				PARVEC_PROFILE_END(); //ParVec_CSTAdditive2Unique // added 14042011ih
				break;
			}
			else if(has_storage_type(PST_CONSISTENT)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTConsistent2Unique); // added 14042011ih
				if(m_pSlaveLayout == NULL)
					UG_THROW("ParallelVector::change_storage_type: No "
							"layouts given but trying to change type.");
				ConsistentToUnique(this, *const_cast<IndexLayout*>(m_pSlaveLayout));
				set_storage_type(PST_ADDITIVE);
				add_storage_type(PST_UNIQUE);
				PARVEC_PROFILE_END(); //ParVec_CSTConsistent2Unique // added 14042011ih
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
	PROFILE_FUNC_GROUP("algebra");
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
bool
ParallelVector<TVector>::
set_random(number from, number to, ParallelStorageType type)
{
	PROFILE_FUNC_GROUP("algebra");
	// set all local vector to value. Therefore parallel vector is consistent
	if(!TVector::set_random(from, to)) return false;
	set_storage_type(PST_ADDITIVE);

	// additive or additive unique
	return change_storage_type(type);
}

template <typename TVector>
inline
number ParallelVector<TVector>::norm() const
{
	PROFILE_FUNC_GROUP("algebra parallelization");
// 	step 1: make vector d additive unique
	if(!const_cast<ParallelVector<TVector>*>(this)->change_storage_type(PST_UNIQUE))
		UG_THROW("ParallelVector::norm(): Cannot change"
							" ParallelStorageType to unique.");

// 	step 2: compute process-local defect norm, square them
	double tNormLocal = (double)TVector::norm();
	tNormLocal *= tNormLocal;

// 	step 3: sum squared local norms
	PARVEC_PROFILE_BEGIN(ParVec_norm_allreduce);
	double tNormGlobal = m_processCommunicator.allreduce(tNormLocal, PCL_RO_SUM);
	PARVEC_PROFILE_END();

// 	step 4: return global norm ( = square root of summed squared local norms)
	return sqrt((number)tNormGlobal);
}

template <typename TVector>
inline
number ParallelVector<TVector>::dotprod(const this_type& v)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
// 	step 0: check that storage type is given
	if(this->has_storage_type(PST_UNDEFINED) || v.has_storage_type(PST_UNDEFINED))
	{
		UG_LOG("ERROR in 'ParallelVector::dotprod': "
				"Parallel storage type of vector not given.\n");
		UG_THROW("ERROR in ParallelVector::dotprod(): No parallel "
				"Storage type given.");
	}

//	step 1: Check if good storage type are given (no communication needed)
//         - additive (unique) <-> consistent is ok
//         - unique <-> unique is ok
	bool check = false;
	if(this->has_storage_type(PST_ADDITIVE)
			&& v.has_storage_type(PST_CONSISTENT)) check = true;
	if(this->has_storage_type(PST_CONSISTENT)
			&& v.has_storage_type(PST_ADDITIVE)) check = true;
	if(this->has_storage_type(PST_UNIQUE)
			&& v.has_storage_type(PST_UNIQUE))     check = true;

// 	step 2: fall back
//         	if storage type not as in the upper cases, communicate to
//			correct solution a user of this function should ideally avoid
//			such a change and do it outside of this function
	if(!check)
	{
		// unique <-> additive => consistent <-> additive
		if(this->has_storage_type(PST_UNIQUE)
				&& v.has_storage_type(PST_ADDITIVE))
			{this->change_storage_type(PST_CONSISTENT);}
		// additive <-> unique => unique <-> unique
		else if(this->has_storage_type(PST_ADDITIVE)
				&& v.has_storage_type(PST_UNIQUE))
			{this->change_storage_type(PST_UNIQUE);}
		// consistent <-> consistent => unique <-> consistent
		else {this->change_storage_type(PST_UNIQUE);}
	}

// 	step 3: compute local dot product
	double tSumLocal = (double)TVector::dotprod(v);
	double tSumGlobal;

// 	step 4: sum global contributions
	m_processCommunicator.allreduce(&tSumLocal, &tSumGlobal, 1,
									PCL_DT_DOUBLE, PCL_RO_SUM);

// 	step 5: return result
	return tSumGlobal;
}




// for template expressions
// d = alpha*v1
template<typename T>
inline void VecScaleAssign(ParallelVector<T> &dest,
		double alpha1, const ParallelVector<T> &v1)
{
	PROFILE_FUNC_GROUP("algebra");
	dest.copy_storage_type(v1);
	VecScaleAssign(*dynamic_cast<T*>(&dest), alpha1, *dynamic_cast<const T*>(&v1));
}

// for template expressions
// d = v1
template<typename T>
inline void VecAssign(ParallelVector<T> &dest,
		const ParallelVector<T> &v1)
{
	PROFILE_FUNC_GROUP("algebra");
	dest.copy_storage_type(v1);
	VecAssign(*dynamic_cast<T*>(&dest), *dynamic_cast<const T*>(&v1));
}


// dest = alpha1*v1 + alpha2*v2
template<typename T>
inline void VecScaleAdd(ParallelVector<T>  &dest,
		double alpha1, const ParallelVector<T> &v1,
		double alpha2, const ParallelVector<T> &v2)
{
	PROFILE_FUNC_GROUP("algebra");
	uint mask = v1.get_storage_mask() & v2.get_storage_mask();
	UG_ASSERT(mask != 0, "VecScaleAdd: cannot add vectors v1 and v2");
	dest.set_storage_type(mask);

	VecScaleAdd(*dynamic_cast<T*>(&dest), 	alpha1, *dynamic_cast<const T*>(&v1),
											alpha2, *dynamic_cast<const T*>(&v2));
}

// dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename T>
inline void VecScaleAdd(ParallelVector<T> &dest,
		double alpha1, const ParallelVector<T> &v1,
		double alpha2, const ParallelVector<T> &v2,
		double alpha3, const ParallelVector<T> &v3)
{
	PROFILE_FUNC_GROUP("algebra");
	uint mask = 	v1.get_storage_mask() &
								v2.get_storage_mask() &
								v3.get_storage_mask();
	UG_ASSERT(mask != 0, "VecScaleAdd: cannot add vectors v1 and v2");
	dest.set_storage_type(mask);

	VecScaleAdd((T&)dest, alpha1, (const T&)v1, alpha2, (const T&)v2, alpha3, (const T&)v3);
}


////////////////////////////////////////////////////////////////////////////////////////

template<typename TVector>
bool CloneVector(ParallelVector<TVector> &dest, const ParallelVector<TVector>& src)
{
	CloneVector((TVector&)dest, (TVector&)src);
	dest.copy_storage_type(src);
	dest.copy_layouts(src);
	return true;
}

template<typename TVector>
ParallelVector<TVector>* ParallelVector<TVector>::virtual_clone() const
{
	ParallelVector<TVector>* pVec = new ParallelVector<TVector>();
	*pVec = *this;
	return pVec;
}

template<typename TVector>
SmartPtr<ParallelVector<TVector> > ParallelVector<TVector>::clone() const
{
	return SmartPtr<ParallelVector<TVector> >(this->virtual_clone());
}

template<typename TVector>
ParallelVector<TVector>* ParallelVector<TVector>::virtual_clone_without_values() const
{
	ParallelVector<TVector>* pVec = new ParallelVector<TVector>(this->size());
	pVec->copy_storage_type(*this); pVec->copy_layouts(*this);
	return pVec;
}

template<typename TVector>
SmartPtr<ParallelVector<TVector> > ParallelVector<TVector>::clone_without_values() const
{
	return SmartPtr<ParallelVector<TVector> >(this->virtual_clone_without_values());
}



} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__ */
