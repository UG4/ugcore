/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_VECTOR_IMPL__

#include "parallel_vector.h"

#include <cmath>

#include "common/common.h" // for UGError
#include "common/profiler/profiler.h" // additions for profiling (14042011ih)
#include "lib_algebra/parallelization/communication_policies.h"
#include "lib_algebra/parallelization/parallelization_util.h"

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

namespace ug {

template <typename TVector>
typename ParallelVector<TVector>::this_type&
ParallelVector<TVector>::operator = (const this_type &v)
{
	//	forward to sequential vectors
	TVector::operator = (*static_cast<const TVector*>(&v));

	//	copy storage type and layouts
	this->set_storage_type(v.get_storage_mask());
	this->set_layouts(v.layouts());

	//	we're done
	return *this;
}

template <typename TVector>
typename ParallelVector<TVector>::this_type&
ParallelVector<TVector>::operator -= (const this_type &v)
{
	//	compute storage mask
	uint mask = get_storage_mask() & v.get_storage_mask();

	//	check mask
	if(mask == 0)
		UG_THROW("ParallelVector::operator -= : Incompatible storage type: "<<
		         get_storage_mask()<<" and "<<v.get_storage_mask());

	//	set this vector to mask
	m_type = mask;

	//	forward
	TVector::operator -= (*static_cast<const TVector*>(&v));

	//	we're done
	return *this;
}

template <typename TVector>
typename ParallelVector<TVector>::this_type&
ParallelVector<TVector>::operator += (const this_type &v)
{
	//	compute parallel storage mask
	uint mask = get_storage_mask() & v.get_storage_mask();

	//	check mask
	if(mask == 0)
		UG_THROW("ParallelVector::operator += : Incompatible storage type: "<<
		         get_storage_mask()<<" and "<<v.get_storage_mask());

	//	set to new mask
	m_type = mask;

	// 	forward to sequential vector
	TVector::operator += (*static_cast<const TVector*>(&v));

	//	we're done
	return *this;
}

template <typename TVector>
bool
ParallelVector<TVector>::
change_storage_type(ParallelStorageType type)
{
	PROFILE_FUNC_GROUP("algebra parallelization");

	// can only change if current state is defined
	if(has_storage_type(PST_UNDEFINED)) [[unlikely]]
		UG_THROW("ParallelVector::change_storage_type: Trying to change"
				" storage type of a vector that has type PST_UNDEFINED.");

	// if already in that type
	if(has_storage_type(type)) return true;

	// check layouts
	if(layouts().invalid())
		UG_THROW("ParallelVector::change_storage_type: No "
					"layouts given but trying to change type.")

	// else switch to that type
	switch(type)
	{
		case PST_CONSISTENT:
			if(has_storage_type(PST_UNIQUE)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTUnique2Consistent);
				UniqueToConsistent(this, layouts()->master(), layouts()->slave(),
				                   &layouts()->comm());
				set_storage_type(PST_CONSISTENT);
				PARVEC_PROFILE_END(); //ParVec_CSTUnique2Consistent
			}
			else if(has_storage_type(PST_ADDITIVE)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTAdditive2Consistent);
				AdditiveToConsistent(this, layouts()->master(), layouts()->slave(),
				                     &layouts()->comm());
				set_storage_type(PST_CONSISTENT);
				PARVEC_PROFILE_END(); //ParVec_CSTAdditive2Consistent
			}
			else return false;

			if(layouts()->overlap_enabled()){
				PARVEC_PROFILE_BEGIN(ParVec_CSTAdditive2Consistent_CopyOverlap);
				CopyValues(this, layouts()->slave_overlap(),
				           layouts()->master_overlap(), &layouts()->comm());
			}

			break;

		case PST_ADDITIVE:
			if(has_storage_type(PST_UNIQUE)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTUnique2Additive);
				add_storage_type(PST_ADDITIVE);
				PARVEC_PROFILE_END(); //ParVec_CSTUnique2Additive
			}
			else if(has_storage_type(PST_CONSISTENT)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTConsistent2Additive);
				ConsistentToUnique(this, layouts()->slave());
				if(layouts()->overlap_enabled())
					SetLayoutValues(this, layouts()->master_overlap(),  typename TVector::value_type(0));
				set_storage_type(PST_ADDITIVE);
				add_storage_type(PST_UNIQUE);
				PARVEC_PROFILE_END(); //ParVec_CSTConsistent2Additive
			}
			else return false;
			break;

		case PST_UNIQUE:
			if(has_storage_type(PST_ADDITIVE)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTAdditive2Unique);
				if(layouts()->overlap_enabled()){
					AdditiveToConsistent(this, layouts()->master(), layouts()->slave(),
				                     	 &layouts()->comm());
					CopyValues(this, layouts()->slave_overlap(),
				           	   layouts()->master_overlap(), &layouts()->comm());
					ConsistentToUnique(this, layouts()->slave());
				}
				else{
					AdditiveToUnique(this, layouts()->master(), layouts()->slave(),
					                 &layouts()->comm());
				}
				add_storage_type(PST_UNIQUE);
				PARVEC_PROFILE_END(); //ParVec_CSTAdditive2Unique
			}
			else if(has_storage_type(PST_CONSISTENT)){
				PARVEC_PROFILE_BEGIN(ParVec_CSTConsistent2Unique);
				if(layouts()->overlap_enabled()){
					CopyValues(this, layouts()->slave_overlap(),
				           	   layouts()->master_overlap(), &layouts()->comm());
				}
				ConsistentToUnique(this, layouts()->slave());
				set_storage_type(PST_ADDITIVE);
				add_storage_type(PST_UNIQUE);
				PARVEC_PROFILE_END(); //ParVec_CSTConsistent2Unique
			}
			else return false;

			if(layouts()->overlap_enabled()){
				SetLayoutValues(this, layouts()->slave_overlap(),   typename TVector::value_type(0));
			}
			break;
		default: return false;
	}
	return true;
}

template <typename TVector>
void
ParallelVector<TVector>::
set(number w, ParallelStorageType type)
{
	PROFILE_FUNC_GROUP("algebra");
	if(type == PST_UNDEFINED) return;

	// set all local vector to value. Therefore parallel vector is consistent
	TVector::set(w);
	set_storage_type(PST_CONSISTENT);

	// if w == 0.0, we have all types
	if(w == 0.0){
		add_storage_type(PST_ADDITIVE);
		add_storage_type(PST_UNIQUE);
		return;
	}

	// consistent required
	if(type & PST_CONSISTENT) return;

	// additive or additive unique
	change_storage_type(PST_UNIQUE);
}

template <typename TVector>
number
ParallelVector<TVector>::
operator = (number d)
{
	this->set(d);
	return d;
	// ø todo check if this value is actually needed to be returnd
	//  ø todo  standard c++ operator= returns this& instead of number
}

template <typename TVector>
void
ParallelVector<TVector>::
set_random(number from, number to, ParallelStorageType type)
{
	PROFILE_FUNC_GROUP("algebra");
	// set all local vector to value. Therefore parallel vector is consistent
	TVector::set_random(from, to);
	set_storage_type(PST_ADDITIVE);

	// additive or additive unique
	change_storage_type(type);
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
	double tNormGlobal;
	if(layouts()->proc_comm().empty())
		tNormGlobal = tNormLocal;
	else
		tNormGlobal = layouts()->proc_comm().allreduce(tNormLocal, PCL_RO_SUM);
	PARVEC_PROFILE_END();

	// 	step 4: return global norm ( = square root of summed squared local norms)
	return sqrt((number)tNormGlobal);
}

template <typename TVector>
inline
number ParallelVector<TVector>::maxnorm() const
{
	PROFILE_FUNC_GROUP("algebra parallelization");

	if(this->has_storage_type(PST_ADDITIVE))
		if(!const_cast<ParallelVector*>(this)->change_storage_type(PST_UNIQUE))
			UG_THROW("ParallelVector::norm(): Cannot change"
					" ParallelStorageType to unique.");

	// 	step 2: compute process-local defect norm, square them
	double tNormLocal = (double)TVector::maxnorm();

	// 	step 3: sum squared local norms
	PARVEC_PROFILE_BEGIN(ParVec_norm_allreduce);
	double tNormGlobal;
	if(layouts()->proc_comm().empty())
		tNormGlobal = tNormLocal;
	else
		tNormGlobal = layouts()->proc_comm().allreduce(tNormLocal, PCL_RO_MAX);
	PARVEC_PROFILE_END();

	// 	step 4: return global norm ( = square root of summed squared local norms)
	return (number)tNormGlobal;
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
	if(layouts()->proc_comm().empty())
		tSumGlobal = tSumLocal;
	else
		layouts()->proc_comm().allreduce(&tSumLocal, &tSumGlobal, 1,
		                                PCL_DT_DOUBLE, PCL_RO_SUM);

	// 	step 5: return result
	return tSumGlobal;
}

template <typename TVector>
void ParallelVector<TVector>::check_storage_type() const
{

//	a) for unique: all slaves must be zero
	if(this->has_storage_type(PST_UNIQUE)){
		//	get slave layout
		const IndexLayout& layout = layouts()->slave();

		//	interface iterators

		auto end = layout.end();

	//	iterate over interfaces
		for(auto layout_iter = layout.begin(); layout_iter != end; ++layout_iter)
		{
		//	get interface
			const IndexLayout::Interface& interface = layout.interface(layout_iter);

		//	loop over indices
			for(auto iter = interface.begin(); iter != interface.end(); ++iter)
			{
			//  get index
				const size_t index = interface.get_element(iter);

			//	set value of vector to zero
				if(BlockNorm2((*this)[index]) != 0.0){
					UG_LOG_ALL_PROCS("Unique vector at slave index "<<index<<
					                 " has block norm: "<<BlockNorm2((*this)[index])<<"\n");
				}
			}
		}
	}

//	c) for consistent: slave values must be equal to master values
	if(this->has_storage_type(PST_CONSISTENT)){
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> com;

	//	step 1: copy master values to slaves
	//	create the required communication policies
		ComPol_CheckConsistency<TVector> cpVecCheckConsistency(this);

	//	perform communication
		com.send_data(layouts()->slave(), cpVecCheckConsistency);
		com.receive_data(layouts()->master(), cpVecCheckConsistency);
		com.communicate();
	}
}

template <typename TVector>
void ParallelVector<TVector>::enforce_consistent_type()
{
	// try to change to unique
	if(this->has_storage_type(PST_ADDITIVE)){
		this->change_storage_type(PST_UNIQUE);
	}

	// if not unique -> make it unique
	// at this point one may change the vector
	if(!this->has_storage_type(PST_UNIQUE)){
		SetLayoutValues(this, layouts()->slave(), typename TVector::value_type(0.0));
		this->set_storage_type(PST_UNIQUE);
	}

	// change to consistent
	this->change_storage_type(PST_CONSISTENT);
}

// for template expressions
// d = alpha*v1
template<typename T>
inline void VecScaleAssign(ParallelVector<T> &dest,
                           double alpha1,
                           const ParallelVector<T> &v1)
{
	PROFILE_FUNC_GROUP("algebra");
	dest.set_storage_type(v1.get_storage_mask());
	VecScaleAssign(*static_cast<T*>(&dest),
						alpha1, *static_cast<const T*>(&v1));
}

// for template expressions
// d = v1
template<typename T>
inline void VecAssign(ParallelVector<T> &dest,
                      const ParallelVector<T> &v1)
{
	PROFILE_FUNC_GROUP("algebra");
	dest.set_storage_type(v1.get_storage_mask());
	VecAssign(*static_cast<T*>(&dest),
					*static_cast<const T*>(&v1));
}


// dest = alpha1*v1 + alpha2*v2
template<typename T>
inline void VecScaleAdd(ParallelVector<T>  &dest,
                        double alpha1, const ParallelVector<T> &v1,
                        double alpha2, const ParallelVector<T> &v2)
{
	PROFILE_FUNC_GROUP("algebra");
	uint mask = v1.get_storage_mask() & v2.get_storage_mask();
	UG_COND_THROW(mask == 0, "VecScaleAdd: cannot add vectors v1 and v2 because their storage masks are incompatible");
	dest.set_storage_type(mask);

	VecScaleAdd(*static_cast<T*>(&dest),
				alpha1, *static_cast<const T*>(&v1),
	            alpha2, *static_cast<const T*>(&v2));
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
	UG_COND_THROW(mask == 0, "VecScaleAdd: cannot add vectors v1 and v2 because their storage masks are incompatible");
	dest.set_storage_type(mask);

	VecScaleAdd((T&)dest,
				alpha1, static_cast<const T &>(v1),
				alpha2, static_cast<const T &>(v2),
				alpha3, static_cast<const T &>(v3));
}

// returns scal<a, b>
template<typename T>
inline double VecProd(const ParallelVector<T> &a, const ParallelVector<T> &b)
{
	return const_cast<ParallelVector<T>* >(&a)->dotprod(b);
}

// Elementwise (Hadamard) product of two vectors
template<typename T>
inline void VecHadamardProd(ParallelVector<T> &dest, const ParallelVector<T> &v1, const ParallelVector<T> &v2)
{
	uint mask = v1.get_storage_mask() & v2.get_storage_mask();
	UG_COND_THROW(mask == 0, "VecHadamardProd: cannot multiply vectors v1 and v2 because their storage masks are incompatible");
	dest.set_storage_type(mask);

	VecHadamardProd(*static_cast<T*>(&dest),
						*static_cast<const T*>(&v1),
						*static_cast<const T*>(&v2));
}

// Elementwise exp of a vector
template<typename T>
inline void VecExp(ParallelVector<T> &dest, const ParallelVector<T> &v)
{
	VecExp(*static_cast<T*>(&dest),
				 *static_cast<const T*>(&v));
}

// Elementwise log (natural logarithm) of a vector
template<typename T>
inline void VecLog(ParallelVector<T> &dest, const ParallelVector<T> &v)
{
	VecLog(*static_cast<T*>(&dest),
				*static_cast<const T*>(&v));
}

////////////////////////////////////////////////////////////////////////////////////////

template<typename TVector>
bool CloneVector(ParallelVector<TVector> &dest, const ParallelVector<TVector>& src)
{
	CloneVector((TVector&)dest, (TVector&)src);
	dest.set_storage_type(src.get_storage_mask());
	dest.set_layouts(src.layouts());
	return true;
}

template<typename TVector>
ParallelVector<TVector>* ParallelVector<TVector>::virtual_clone() const
{
	auto* pVec = new ParallelVector();
	*pVec = *this;
	return pVec;
}

template<typename TVector>
SmartPtr<ParallelVector<TVector> > ParallelVector<TVector>::clone() const
{
	return SmartPtr<ParallelVector >(this->virtual_clone());
}

template<typename TVector>
ParallelVector<TVector>* ParallelVector<TVector>::virtual_clone_without_values() const
{
	auto* pVec = new ParallelVector(this->size());
//	pVec->set_storage_type(this->get_storage_mask()); --> undefined is correct
	pVec->set_layouts(this->layouts());
	return pVec;
}

template<typename TVector>
SmartPtr<ParallelVector<TVector> > ParallelVector<TVector>::clone_without_values() const
{
	return SmartPtr<ParallelVector >(this->virtual_clone_without_values());
}



} // end namespace ug

#endif