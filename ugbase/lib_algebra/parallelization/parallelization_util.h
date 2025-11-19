/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__

#include <utility>
#include <vector>
#include <map>
#include "common/assert.h"
#include "algebra_id.h"
#include "communication_policies.h"
#include "algebra_layouts.h"

// additions for profiling (18042011ih)
#include "common/profiler/profiler.h"
#define PROFILE_PARALLELIZATION_UTIL
#ifdef PROFILE_PARALLELIZATION_UTIL
	#define PU_PROFILE_FUNC()	PROFILE_FUNC()
	#define PU_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define PU_PROFILE_END(name)	PROFILE_END_(name)
#else
	#define PU_PROFILE_FUNC()
	#define PU_PROFILE_BEGIN(name)
	#define PU_PROFILE_END(name)
#endif
// additions for profiling - end

namespace ug{

/**
 * \brief Util Functions for parallel Algebra
 *
 * In order to simplify the use of Communication Policies for Algebras some
 * util functions are provided.
 *
 * \defgroup lib_algebra_parallelization_util Parallel Algebra Util
 * \ingroup lib_algebra_parallelization
 */

/// \ingroup lib_algebra_parallelization_util
/// @{


///	Generates a set of unique global algebra ids.
/**	Horizontal slaves have the same algebra-id as their associated
 * horizontal masters.
 * Make sure that masterLayout and slaveLayout do not reference
 * indices >= numIDs.
 */
template <class TLayout>
void GenerateGlobalAlgebraIDs(pcl::InterfaceCommunicator<TLayout>& communicator,
		std::vector<AlgebraID>& idsOut,
		size_t numIDs,
		const TLayout& masterLayout,
		const TLayout& slaveLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	generate an id for each entry.
	idsOut.resize(numIDs);
	int localProc = pcl::ProcRank();
	for(size_t i = 0; i < numIDs; ++i)
		idsOut[i] = AlgebraID(localProc, i);

//	copy all ids from master to slave interfaces
	ComPol_VecCopy<std::vector<AlgebraID> >	copyPol(&idsOut);

	communicator.send_data(masterLayout, copyPol);
	communicator.receive_data(slaveLayout, copyPol);
	communicator.communicate();

//	a set of global ids has now been generated.
}


template <typename TMatrix>
void MatAddSlaveRowsToMasterRowOverlap0(TMatrix& mat)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	using namespace std;
	vector<AlgebraID> globalIDs;
	const IndexLayout& masters = mat.layouts()->master();
	const IndexLayout& slaves = mat.layouts()->slave();
	pcl::InterfaceCommunicator<IndexLayout>& comm = mat.layouts()->comm();

	GenerateGlobalAlgebraIDs(comm, globalIDs, mat.num_rows(), masters, slaves);

//	global ids are applied, now communicate...
	ComPol_MatAddRowsOverlap0<TMatrix> comPolMatAdd(mat, globalIDs);
	comm.send_data(slaves, comPolMatAdd);
	comm.receive_data(masters, comPolMatAdd);
	comm.communicate();
}


template <typename TMatrix>
void MatMakeConsistentOverlap0(TMatrix& mat)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	using namespace std;
	vector<AlgebraID> globalIDs;
	const IndexLayout& masters = mat.layouts()->master();
	const IndexLayout& slaves = mat.layouts()->slave();
	pcl::InterfaceCommunicator<IndexLayout>& comm = mat.layouts()->comm();

	GenerateGlobalAlgebraIDs(comm, globalIDs, mat.num_rows(), masters, slaves);

//	global ids are applied, now communicate...
	ComPol_MatAddRowsOverlap0<TMatrix> comPolMatAdd(mat, globalIDs);
	comm.send_data(slaves, comPolMatAdd);
	comm.receive_data(masters, comPolMatAdd);
	comm.communicate();

	ComPol_MatCopyRowsOverlap0<TMatrix> comPolMatCopy(mat, globalIDs);
	comm.send_data(masters, comPolMatCopy);
	comm.receive_data(slaves, comPolMatCopy);
	comm.communicate();
}

/// changes parallel storage type from additive to consistent
/**
 * This function changes the storage type of a parallel vector from additive
 * to consistent. A InterfaceCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void AdditiveToConsistent(	TVector* pVec,
                          	const IndexLayout& masterLayout, const IndexLayout& slaveLayout,
                          	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	step 1: add slave values to master
	//	create the required communication policies
		ComPol_VecAdd<TVector> cpVecAdd(pVec);

		PU_PROFILE_BEGIN(AdditiveToConsistent_step1);
	//	perform communication
		com.send_data(slaveLayout, cpVecAdd);
		com.receive_data(masterLayout, cpVecAdd);
		com.communicate();
		PU_PROFILE_END(AdditiveToConsistent_step1);
	//	step 2: copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVec);

		PU_PROFILE_BEGIN(AdditiveToConsistent_step2);
	//	perform communication
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
		PU_PROFILE_END(AdditiveToConsistent_step2);
}

/// changes parallel storage type from unique to consistent
/**
 * This function changes the storage type of a parallel vector from unique
 * to consistent. A InterfaceCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void UniqueToConsistent(	TVector* pVec,
							const IndexLayout& masterLayout, const IndexLayout& slaveLayout,
							pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	create a new communicator if required.
	pcl::InterfaceCommunicator<IndexLayout> tCom;
	if(!pCom)
		pCom = &tCom;
	pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

//	step 1: copy master values to slaves
//	create the required communication policies
	ComPol_VecCopy<TVector> cpVecCopy(pVec);

//	perform communication
	com.send_data(masterLayout, cpVecCopy);
	com.receive_data(slaveLayout, cpVecCopy);
	com.communicate();
}


///	Copies values from the source to the target layout
template <typename TVector>
void CopyValues(	TVector* pVec,
					const IndexLayout& sourceLayout, const IndexLayout& targetLayout,
					pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	create a new communicator if required.
	pcl::InterfaceCommunicator<IndexLayout> tCom;
	if(!pCom)
		pCom = &tCom;
	pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

//	step 1: copy master values to slaves
//	create the required communication policies
	ComPol_VecCopy<TVector> cpVecCopy(pVec);

//	perform communication
	com.send_data(sourceLayout, cpVecCopy);
	com.receive_data(targetLayout, cpVecCopy);
	com.communicate();
}


/// changes parallel storage type from additive to unique
/**
 * This function changes the storage type of a parallel vector from additive
 * to unique. A InterfaceCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void AdditiveToUnique(	TVector* pVec,
						const IndexLayout& masterLayout, const IndexLayout& slaveLayout,
						pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	step 1: add slave values to master and set slave values to zero
	//	create the required communication policies
		ComPol_VecAddSetZero<TVector> cpVecAddSetZero(pVec);

	//	perform communication
		com.send_data(slaveLayout, cpVecAddSetZero);
		com.receive_data(masterLayout, cpVecAddSetZero);
		com.communicate();
}

/// sets the values of a vector to a given number only on the interface indices
/**
 * \param[in,out]		pVec			Vector
 * \param[in]			layout			Layout
 * \param[in]			val				Value to set on layout indices
 */
template <typename TVector>
void SetInterfaceValues(TVector* pVec,
                     	const IndexLayout::Interface& interface,
                     	typename TVector::value_type val)
{
	PROFILE_FUNC_GROUP("algebra parallelization");

//	loop over indices
	for(typename IndexLayout::Interface::const_iterator iter = interface.begin();
			iter != interface.end(); ++iter)
	{
	//  get index
		const size_t index = interface.get_element(iter);

	//	set value of vector to zero
		(*pVec)[index] = val;
	}
}


/// sets the values of a vector to a given number only on the layout indices
/**
 * \param[in,out]		pVec			Vector
 * \param[in]			layout			Layout
 * \param[in]			val				Value to set on layout indices
 */
template <typename TVector>
void SetLayoutValues(	TVector* pVec,
                     	const IndexLayout& layout,
                     	typename TVector::value_type val)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	interface iterators
	typename IndexLayout::const_iterator iter = layout.begin();
	typename IndexLayout::const_iterator end = layout.end();

//	iterate over interfaces
	for(; iter != end; ++iter)
	{
	//	get interface
		const typename IndexLayout::Interface& interface = layout.interface(iter);

	//	loop over indices
		for(typename IndexLayout::Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		//  get index
			const size_t index = interface.get_element(iter);

		//	set value of vector to zero
			(*pVec)[index] = val;
		}
	}
}

/// scales the values of a vector by a given number only on the layout indices
/**
 * \param[in,out]		pVec			Vector
 * \param[in]			layout			Layout
 * \param[in]			scale			Scaling factor
 */
template <typename TVector>
void ScaleLayoutValues(	TVector* pVec,
                     	const IndexLayout& layout,
                     	number scale)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	interface iterators
	typename IndexLayout::const_iterator iter = layout.begin();
	typename IndexLayout::const_iterator end = layout.end();

//	iterate over interfaces
	for(; iter != end; ++iter)
	{
	//	get interface
		const typename IndexLayout::Interface& interface = layout.interface(iter);

	//	loop over indices
		for(typename IndexLayout::Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		//  get index
			const size_t index = interface.get_element(iter);

		//	set value of vector to zero
			(*pVec)[index] *= scale;
		}
	}
}


/// changes parallel storage type from consistent to unique
/**
 * This function changes the storage type of a parallel vector from consistent
 * to unique. Note, that no communication is needed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			slaveLayout		Slave Layout
 */
template <typename TVector>
void ConsistentToUnique(	TVector* pVec,
                        	const IndexLayout& slaveLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	SetLayoutValues(pVec, slaveLayout, 0.0);
}

/// subtracts values of slave layout from master layout and sets slave layouts to negative of difference
/**
 * This function subtracts all slave values from the master value. Then, the
 * slave values are set to the negative of the computed difference.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void VecSubtractOnLayout(	TVector* pVec,
                         	const IndexLayout& masterLayout, const IndexLayout& slaveLayout,
                         	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	step 1: subtract slave values from master
	//	create the required communication policies
		ComPol_VecSubtract<TVector> cpVecSubtract(pVec);

	//	Subtract slave values from master values
		com.send_data(slaveLayout, cpVecSubtract);
		com.receive_data(masterLayout, cpVecSubtract);
		com.communicate();

	//	step 2: Copy values to slaves
		ComPol_VecScaleCopy<TVector> cpVecScaleCopy(pVec, -1.0);

		com.send_data(masterLayout, cpVecScaleCopy);
		com.receive_data(slaveLayout, cpVecScaleCopy);
		com.communicate();
}

/// subtracts values of only one slave dof per master on layout
/**
 * This function subtracts one slave value dof per master dof. Communication is performed.
 * Information flow is from slaves to masters but not vice versa.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void VecSubtractOneSlaveFromMaster(	TVector* pVec,
                                   	const IndexLayout& masterLayout,
                                   	const IndexLayout& slaveLayout,
                                   	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	create the required communication policies
		ComPol_VecSubtractOnlyOneSlave<TVector> cpVecSubtractOOS(pVec);

	//	sending: slaves, receiving: masters; masters subtract the value of only
	//	one slave on reception (according to the policy used)
		com.send_data(slaveLayout, cpVecSubtractOOS);
		com.receive_data(masterLayout, cpVecSubtractOOS);
		com.communicate();

}


/// copy a vector only at a layout
/**
 * This function copies a vector only at a layout. Communication is performed.
 * Information flow is from masters to slaves but not vice versa.
 * Note: This function is basically 'UniqueToConsistent()', so "sender" and
 * "receiver" in this function are the same as in 'UniqueToConsistent()').
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void VecCopy(	TVector* pVec,
             	const IndexLayout& masterLayout, const IndexLayout& slaveLayout,
             	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVec);

	//	perform communication
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
}

//	returns the highest referenced index of the elements in the layout.
int GetHighestReferencedIndex(IndexLayout& layout);

////////////////////////////////////////////////////////////////////////
///	fills a connection list, which gives the connected processes to each interface.
/**	The i-th entry in connectionsOut is a vector which contains the process-ids
 * of all connected processes of the layout-entry with value i.
 * This includes connections to all slaves of the associated master.
 * The first entry for each connection is the process on which the
 * master-element lies, followed by the processes where associated slaves lie.
 *
 * \param	connectionsOut will have the size highestReferencedIndex + 1 when the
 * 			method is done. By indexing connectionsOut with an entry of one of the
 * 			layouts interfaces, you will obtain a vector of all processes on
 * 			which copies of that entry reside.
 *
 * \param 	highestReferencedIndex has to hold the highest index which is referenced
 * 			either by entries in the interfaces of masterLayout or by entries in the
 * 			interfaces of slaveLayout. Note that it can be obtained by a call to
 * 			std::max(GetHighestReferencedIndex(masterLayout), GetHighestReferencedIndex(slaveLayout)).
 *
 * \todo	The method should probably automatically find the highestReferencedIndex,
 * 			to avoid misuse.
 */
void CommunicateConnections(std::vector<std::vector<int> >& connectionsToProcsOut,
							std::vector<std::vector<int> >& connectionsToSubDomsOut,
							IndexLayout& masterLayout,
							IndexLayout& slaveLayout,
							int highestReferencedIndex, pcl::IDomainDecompositionInfo& ddinfo);

/**	given a layout which defines relations between neighbours, this method
 * creates a layout which connect the elements in the given layouts with
 * newly created elements on a root process and establishes a OneToMany connection.
 *
 * All newly created elements will be contained in the masterInterfaceOut on rootProc.
 * The associated vector-indices are assigned so that the work with a mxm matrix,
 * where m is the returned number of newly created elements.
 *
 * \param	pNewMasterIDsOut (out)an optional vector to which the
 * 			associated new element ids for each
 * 			element in the given masterLayout and slaveLayout will be written.
 * 			If not large enough it will be resized to the highest
 * 			referenced index and will contain -1 in each entry which
 * 			is not referenced by masterLayout or slaveLayout.
 * \return	the number of newly created elements (!= 0 only on rootProc).
 */
int BuildOneToManyLayout(IndexLayout& masterLayoutOut,
						  IndexLayout& slaveLayoutOut,
						  int rootProcID,
						  IndexLayout& masterLayout,
						  IndexLayout& slaveLayout,
						  pcl::ProcessCommunicator procComm,
						  std::vector<int>* pNewMasterIDsOut = nullptr);



/**	Given standard master and slave index layouts (as created e.g. by
 * ug::CreateIndexLayout), this method constructs layouts as required
 * by domain decomposition methods such as the Feti solver.
 *
 * Note that the created layouts are horizontal layouts which allow
 * to communicate between neighbours of the same subdomain (process-layouts),
 * between neighbours in different subdomains (subdom-layouts) and
 * between neighbours in the same subdomain, which lie in a subdomain layout
 * (deltaNbr-layouts).
 */
void BuildDomainDecompositionLayouts(
		IndexLayout& subdomMastersOut, IndexLayout& subdomSlavesOut,
		IndexLayout& processMastersOut, IndexLayout& processSlavesOut,
		IndexLayout& deltaNbrMastersOut, IndexLayout& deltaNbrSlavesOut,
		IndexLayout& crossPointMastersOut, IndexLayout& crossPointSlavesOut,
		const IndexLayout& standardMasters, const IndexLayout& standardSlaves,
		int highestReferencedIndex, pcl::IDomainDecompositionInfo& ddinfo);
/// @}

/// extracts diagonal of a matrix for interface indices
/**
 * This function extracts diagonal values of a matrix for all interface indices.
 * No communication is performed.
 *
 * \param[out]		pDiagVector	Vector with diagonal entries
 * \param[in]		pMatrix		Matrix
 * \param[in]		Layout		Index Layout
 */
template <typename TMatrix, typename TVector>
void MatExtractDiagOnLayout(	TVector* pDiagVector,
								const TMatrix* pMatrix,
								const IndexLayout& Layout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	interface iterator
	typename IndexLayout::const_iterator iter = Layout.begin();
	typename IndexLayout::const_iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
	//	get interface
		const typename IndexLayout::Interface& interface = Layout.interface(iter);

		for(typename IndexLayout::Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		// 	get index
			const size_t index = interface.get_element(iter);

		//	copy values
			const typename TMatrix::value_type& block = (*pMatrix)(index, index);
			for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
			{
				BlockRef((*pDiagVector)[index], beta) = BlockRef(block, beta, beta);
			}
		}
	}
}

/// writes diagonal of a matrix for interface indices
/**
 * This function writes diagonal values of a matrix for all interface indices.
 * No communication is performed.
 *
 * \param[out]		pMatrix		Matrix
 * \param[in]		pDiagVector	Vector with diagonal entries
 * \param[in]		Layout		Index Layout
 */
template <typename TMatrix, typename TVector>
void MatWriteDiagOnLayout(	TMatrix* pMatrix,
                          	const TVector* pDiagVector,
                          	const IndexLayout& Layout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	interface iterator
	typename IndexLayout::const_iterator iter = Layout.begin();
	typename IndexLayout::const_iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
	//	get interface
		const typename IndexLayout::Interface& interface = Layout.interface(iter);

		for(typename IndexLayout::Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		// 	get index
			const size_t index = interface.get_element(iter);

		//	copy values
			typename TMatrix::value_type& block = (*pMatrix)(index, index);
			for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
			{
				BlockRef(block, beta, beta) = BlockRef((*pDiagVector)[index], beta);
			}
		}
	}
}


/// changes parallel storage type from additive to consistent on diagonal of a matrix
/**
 * This function changes the storage type of a matrix from additive
 * to consistent on the diagonal. A InterfaceCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TAlgebra>
void MatAdditiveToConsistentOnDiag(	typename TAlgebra::matrix_type* pMat,
                                   	const IndexLayout& masterLayout, const IndexLayout& slaveLayout,
                                   	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	\todo: We could work on the matrix directly here, without temporary vector

//	create a vector of length of the diagonal
	typename TAlgebra::vector_type vecDiag;

//	resize the vector to correct size
	vecDiag.resize(pMat->num_rows());

//	copy diag values
	MatExtractDiagOnLayout(&vecDiag, pMat, masterLayout);
	MatExtractDiagOnLayout(&vecDiag, pMat, slaveLayout);

//	change vector to consistent
	AdditiveToConsistent(&vecDiag, masterLayout, slaveLayout, pCom);

//	write consistent values back
	MatWriteDiagOnLayout(pMat, &vecDiag, masterLayout);
	MatWriteDiagOnLayout(pMat, &vecDiag, slaveLayout);
}

/// gathers all values in master indices of a second vector
/**
 * This function gathers the slave values from a source vector in the master
 * values of a second vector. Slave values are added to the current master values.
 * A InterfaceCommunicator is created iff no communicator passed.
 *
 * \param[out]			pVecDest		destination Parallel Vector
 * \param[in]			pVecSrc			source Parallel Vector
 * \param[in]			masterLayouDest Master Layout
 * \param[in]			slaveLayoutSrc	Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void VecGather(	TVector* pVecDest, const TVector* pVecSrc,
               	const IndexLayout& masterLayoutDest, const IndexLayout& slaveLayoutSrc,
               	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	step 1: add slave values to master
	//	create the required communication policies
		ComPol_VecAdd<TVector> cpVecAdd(pVecDest, pVecSrc);

	//	perform communication
		com.send_data(slaveLayoutSrc, cpVecAdd);
		com.receive_data(masterLayoutDest, cpVecAdd);
		com.communicate();
}

/// broadcasts all values from master indices to slave values in a second vector
/**
 * This function broadcasts the master values from a source vector to the slave
 * values of a second vector. Slave values are overwritten.
 * A InterfaceCommunicator is created iff no communicator passed.
 *
 * \param[out]			pVecDest		destination Parallel Vector
 * \param[in]			pVecSrc			source Parallel Vector
 * \param[in]			masterLayouDest Master Layout
 * \param[in]			slaveLayoutSrc	Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void VecBroadcast(	TVector* pVecDest, const TVector* pVecSrc,
                  	const IndexLayout& slaveLayoutDest, const IndexLayout& masterLayoutSrc,
                  	pcl::InterfaceCommunicator<IndexLayout>* pCom = nullptr)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	step 1: copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVecDest, pVecSrc);

	//	perform communication
		com.send_data(masterLayoutSrc, cpVecCopy);
		com.receive_data(slaveLayoutDest, cpVecCopy);
		com.communicate();
}

inline bool PrintLayouts(const HorizontalAlgebraLayouts &layout)
{
	return TestLayout(layout.proc_comm(), layout.comm(), layout.master(), layout.slave(), true);
}

/**
 * @param layout layout to search
 * @param pid proc_id to find
 * @return layout.end() if not found, else it so that layout.proc_id(it) == pid
 */
template<typename TLayout>
typename TLayout::iterator find_pid(TLayout &layout, int pid)
{
	for(typename TLayout::iterator it = layout.begin(); it != layout.end(); ++it)
		if(layout.proc_id(it) == pid) return it;
	return layout.end();
}

/**
 * @return an empty AlgebraLayout with proc_comm pcl::PCD_LOCAL
 */
SmartPtr<AlgebraLayouts> CreateLocalAlgebraLayouts();


///	Generates a set of global consecutive indices
/**	TIndVec has to be an std::vector compatible type, e.g. std::vector<size_t>.*/
template <class TIndVec>
void GenerateGlobalConsecutiveIndices(TIndVec& indsOut, size_t numLocalInds,
									  const AlgebraLayouts& layouts);

///	Tests layouts by matching master and slave interfaces and by comparing global id's
/**	Checks normal H-Master and H-Slave interfaces as well as H-Master-Overlap and
 * H-Slave-Overlap interfaces.
 *
 * \param mat			A parallel matrix.
 * \param algebraIDs	may be specified optionally. Make sure that those are global
 *						consistent algebraID's-*/
template <class TMatrix>
void TestHorizontalAlgebraLayouts(
			const TMatrix& mat,
			std::vector<AlgebraID>* algebraIDs = nullptr,
			bool verbose = false);
}//	end of namespace


///	include implementation
#include "parallelization_util_impl.h"

#endif