/*
 * parallelization_util.h
 *
 *  Created on: 14.6.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__

#include <utility>
#include <vector>
#include <map>
#include "common/assert.h"
#include "communication_policies.h"

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

///	this type is used to identify distributed objects.
typedef std::pair<int, size_t> AlgebraID;

///	Generates a set of global algebra ids.
/**	Make sure that masterLayout and slaveLayout do not reference
 * indices >= numIDs.
 */
template <class TLayout>
void GenerateGlobalAlgebraIDs(std::vector<AlgebraID>& idsOut,
							  size_t numIDs,
							  TLayout& masterLayout,
							  TLayout& slaveLayout)
{
//	generate an id for each entry.
	idsOut.resize(numIDs);
	int localProc = pcl::GetProcRank();
	for(size_t i = 0; i < numIDs; ++i)
		idsOut[i] = std::make_pair(localProc, i);

//	copy all ids from master to slave interfaces
	ComPol_VecCopy<std::vector<AlgebraID> >	copyPol(&idsOut);

	pcl::ParallelCommunicator<TLayout> communicator;
	communicator.send_data(masterLayout, copyPol);
	communicator.receive_data(slaveLayout, copyPol);
	communicator.communicate();

//	a set of global ids has now been generated.
}

/// changes parallel storage type from additive to consistent
/**
 * This function changes the storage type of a parallel vector from additive
 * to consistent. A ParallelCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void AdditiveToConsistent(	TVector* pVec,
							IndexLayout& masterLayout, IndexLayout& slaveLayout,
							pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::ParallelCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

	//	step 1: add slave values to master
	//	create the required communication policies
		ComPol_VecAdd<TVector> cpVecAdd(pVec);

	//	perform communication
		com.send_data(slaveLayout, cpVecAdd);
		com.receive_data(masterLayout, cpVecAdd);
		com.communicate();

	//	step 2: copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVec);

	//	perform communication
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
}

/// changes parallel storage type from unique to consistent
/**
 * This function changes the storage type of a parallel vector from unique
 * to consistent. A ParallelCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void UniqueToConsistent(	TVector* pVec,
							IndexLayout& masterLayout, IndexLayout& slaveLayout,
							pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::ParallelCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

	//	step 1: copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVec);

	//	perform communication
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
}


/// changes parallel storage type from additive to unique
/**
 * This function changes the storage type of a parallel vector from additive
 * to unique. A ParallelCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TVector>
void AdditiveToUnique(	TVector* pVec,
						IndexLayout& masterLayout, IndexLayout& slaveLayout,
						pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::ParallelCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

	//	step 1: add slave values to master and set slave values to zero
	//	create the required communication policies
		ComPol_VecAddSetZero<TVector> cpVecAddSetZero(pVec);

	//	perform communication
		com.send_data(slaveLayout, cpVecAddSetZero);
		com.receive_data(masterLayout, cpVecAddSetZero);
		com.communicate();
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
							IndexLayout& slaveLayout)
{
//	interface iterators
	typename IndexLayout::iterator iter = slaveLayout.begin();
	typename IndexLayout::iterator end = slaveLayout.end();

//	iterate over interfaces
	for(; iter != end; ++iter)
	{
	//	get interface
		typename IndexLayout::Interface& interface = slaveLayout.interface(iter);

	//	loop over indices
		for(typename IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		//  get index
			const size_t index = interface.get_element(iter);

		//	set value of vector to zero
			(*pVec)[index] = 0.0;
		}
	}
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
							IndexLayout& masterLayout, IndexLayout& slaveLayout,
							pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::ParallelCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

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


/// adds one vector to another only at a layout
/**
 * This function adds to vectors only on a layout. No communication is performed.
 *
 * \param[in,out]	pVecDest		Destination vector of sum
 * \param[in]		pVecSrc			vector to be added
 * \param[in]		scale			Scaling factor
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecScaleAddOnLayout(	TVector* pVecDest, const TVector* pVecSrc,
							number scale, IndexLayout& Layout)
{
//	interface iterators
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
	//	get interface
		typename IndexLayout::Interface& interface = Layout.interface(iter);

		for(typename IndexLayout::Interface::iterator iter = interface.begin();
			iter != interface.end(); ++iter)
		{
		//	get index
			const size_t index = interface.get_element(iter);

		// 	get entry
			typename TVector::value_type entry = (*pVecSrc)[index];

		//	scale entry
			entry *= scale;

		//	add value
			(*pVecDest)[index] += entry;
		}
	}
}

/// sets entries of one vector to a given value only at a layout
/**
 * This function sets the entries of a vector to a given value only at a layout. No communication is performed.
 *
 * \param[in,out]	pVec			vector to set
 * \param[in]		value			value for entries to be set
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecSetOnLayout(	TVector* pVec, number value, IndexLayout& Layout)
{
//	interface iterators
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
	//	get interface
		typename IndexLayout::Interface& interface = Layout.interface(iter);

	//	loop over indices
		for(typename IndexLayout::Interface::iterator iter = interface.begin();
			iter != interface.end(); ++iter)
		{
		//	get index
			const size_t index = interface.get_element(iter);

		//	set value
			(*pVec)[index] = value;
		}
	}
}

/// sets entries of one vector to a given value *except* at the interface
/**
 * This function sets the entries of a vector to a given value *except* on a layout. No communication is performed.
 *
 * \param[in,out]	pVec			vector to set
 * \param[in]		value			value for entries to be set
 * \param[in]		Layout			Index Layout to exclude
 */
template <typename TVector>
void VecSetExcludingLayout(	TVector* pVec, number value,
							IndexLayout& Layout)
{
//	interface iterators
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();


	for(size_t i=0; i < pVec->size(); i++)
	{
		if (iter == end)
		{
		//	set value
			(*pVec)[i] = value;
		} else
		{
			for(; iter != end; ++iter)
			{
			//	get interface
				typename IndexLayout::Interface& interface = Layout.interface(iter);
	
			//	loop over indices
				for(typename IndexLayout::Interface::iterator iter = interface.begin();
					iter != interface.end(); ++iter)
				{
				//	get index
					const size_t index = interface.get_element(iter);
	
					if (index == i) continue;
	
				//	set value
					(*pVec)[i] = value;
				}
			}
		}
	}
}

/// scale entries of one vector by a given value only at a layout
/**
 * This function scales the entries of a vector to a given value only at a layout. No communication is performed.
 *
 * \param[in,out]	pVec			vector to set
 * \param[in]		scale			Scaling factor
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecScaleOnLayout(	TVector* pVec, number scale,
						IndexLayout& Layout)
{
//	interface iterators
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
	//	get interface
		typename IndexLayout::Interface& interface = Layout.interface(iter);

	//	loop over indices
		for(typename IndexLayout::Interface::iterator iter = interface.begin();
			iter != interface.end(); ++iter)
		{
		//	get index
			const size_t index = interface.get_element(iter);

		//	set value
			(*pVec)[index] *= scale;
		}
	}
}

/// scale entries of one vector by a given value *except* at the interface
/**
 * This function scales the entries of a vector to a given value *except* on a layout. No communication is performed.
 *
 * \param[in,out]	pVec			vector to set
 * \param[in]		scale			Scaling factor
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecScaleExcludingLayout(	TVector* pVec, number scale,
								IndexLayout& Layout)
{
//	interface iterators
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();


	for(size_t i=0; i < pVec->size(); i++)
	{
		if (iter == end)
		{
		//	set value
			(*pVec)[i] *= scale;
		} else
		{
			for(; iter != end; ++iter)
			{
			//	get interface
				typename IndexLayout::Interface& interface = Layout.interface(iter);
	
			//	loop over indices
				for(typename IndexLayout::Interface::iterator iter = interface.begin();
					iter != interface.end(); ++iter)
				{
				//	get index
					const size_t index = interface.get_element(iter);
	
					if (index == i) continue;
	
				//	set value
					(*pVec)[i] *= scale;
				}
			}
		}
	}
}

/// sets dirichlet rows for interface indices
/**
 * This function sets matrix rows to identity for all interface indices.
 * No communication is performed.
 *
 * \param[in,out]	pMatrix		Matrix
 * \param[in]		Layout		Index Layout
 */
template <typename TMatrix>
void MatSetDirichletOnLayout(	TMatrix* pMatrix,
											IndexLayout& Layout)
{
//	interface iterator
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
	//	get interface
		typename IndexLayout::Interface& interface = Layout.interface(iter);

		for(typename IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		// 	get index
			const size_t index = interface.get_element(iter);

		//	set identity row
			SetDirichletRow(*pMatrix, index);
		}
	}
}


/**	given a layout which defines relations between neighbours, this method
 * creates a layout which connect the elements in the given layouts with
 * newly created elements on a root process and establishes a OneToMany connection.
 *
 * All newly created elements will be contained in the masterInterfaceOut on rootProc.
 * The associated vector-indices are assigned so that the work with a mxm matrix,
 * where m is the returned number of newly created elements.
 *
 * \param 	highestReferencedIndex The highest index which is referenced by
 * 			masterLayout or by slaveLayout. -1 if there is none.
 * \param	pNewMasterIDsOut (out)an optional vector to whichthe
 * 			associated new element ids for each
 * 			element in the given masterLayout and slaveLayout will be written.
 * 			It will be resized to highestReferencedIndex+1 and will contain
 * 			-1 in each entry which is not referenced by masterLayout or slaveLayout.
 * \return	the number of newly created elements (!= 0 only on rootProc).
 */
int BuildOneToManyLayout(IndexLayout& masterLayoutOut,
						  IndexLayout& slaveLayoutOut,
						  int rootProcID,
						  IndexLayout& masterLayout,
						  IndexLayout& slaveLayout,
						  int highestReferencedIndex,
						  pcl::ProcessCommunicator procComm,
						  std::vector<int>* pNewMasterIDsOut = NULL);



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
		IndexLayout& standardMasters, IndexLayout& standardSlaves,
		int highestReferencedIndex, pcl::IDomainDecompositionInfo& ddinfo);
/// @}

}//	end of namespace


#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__ */
