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
void VecSubtractOneSlaveFromMasterOnLayout(	TVector* pVec,
											IndexLayout& masterLayout,
											IndexLayout& slaveLayout,
											pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::ParallelCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

	//	create the required communication policies
		ComPol_VecSubtractOnlyOneSlave<TVector> cpVecSubtractOOS(pVec);

	//	sending: slaves, receiving: masters; masters subtract the value of only
	//	one slave on reception (according to the policy used)
		com.send_data(slaveLayout, cpVecSubtractOOS);
		com.receive_data(masterLayout, cpVecSubtractOOS);
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
void VecScaleAppendOnLayout(	TVector* pVecDest, const TVector* pVecSrc,
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

/// computes linear combination of two vectors.
/**
 * This function computes linear combination of two vectors (cf. 'VecScaleAdd()'
 * of parallel vector): dest = alpha1*v1 + alpha2*v2
 * (overloaded version of 'VecScaleAddOnLayout()' above)
 *
 * \param[in,out]	pVecDest		Destination vector of sum
 * \param[in]		pVecSrc1		first vector to be combined
 * \param[in]		pVecSrc2		second vector to be combined
 * \param[in]		alpha1			Scaling factor of src vec 1
 * \param[in]		alpha2			Scaling factor of src vec 2
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecScaleAddOnLayout(	TVector* pVecDest,
							number alpha1, const TVector* pVecSrc1,
							number alpha2, const TVector* pVecSrc2,
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

		// 	get entries
			typename TVector::value_type entry1 = (*pVecSrc1)[index];
			typename TVector::value_type entry2 = (*pVecSrc2)[index];

		//	scale entries
			entry1 *= alpha1;
			entry2 *= alpha2;

		//	add value
			entry1 += entry2;
			(*pVecDest)[index] = entry1; //alpha1*entry1 + alpha2*entry2;
		}
	}
}

/// computes the scalar product of two vectors only at a layout.
/**
 * This function the scalar product of two vectors only at a layout 
 * (cf. 'VecProd()' of parallel vector).
 * No communication is performed.
 *
 * \param[out]		result			scalar product computed
 * \param[in]		pVecSrc1		first vector
 * \param[in]		pVecSrc2		second vector
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecProdOnLayout(	number& result, const TVector* pVecSrc1,
						const TVector* pVecSrc2,
						IndexLayout& Layout)
{
//	interface iterators
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	result = 0.0;
	for(; iter != end; ++iter)
	{
	//	get interface
		typename IndexLayout::Interface& interface = Layout.interface(iter);

		for(typename IndexLayout::Interface::iterator iter = interface.begin();
			iter != interface.end(); ++iter)
		{
		//	get index
			const size_t index = interface.get_element(iter);

		// 	get entries
			const typename TVector::value_type& entry1 = (*pVecSrc1)[index];
			const typename TVector::value_type& entry2 = (*pVecSrc2)[index];

		//	add value
			result += VecProd(entry1,entry2);
		}
	}
}

/// sets entries of one vector to a given value only at a layout
/**
 * This function sets the entries of a vector to a given value only at a layout.
 * No communication is performed.
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
 * This function sets the entries of a vector to a given value *except* on a layout.
 * No communication is performed.
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
			for(iter = Layout.begin(); iter != end; ++iter)
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
 * This function scales the entries of a vector to a given value only at a layout.
 * No communication is performed.
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
 * This function scales the entries of a vector to a given value *except* on a layout.
 * No communication is performed.
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
void VecCopyOnLayout(	TVector* pVec,
						IndexLayout& masterLayout, IndexLayout& slaveLayout,
						pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::ParallelCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

	//	copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVec);

	//	perform communication
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
}

/// scale a vector and copy the result to another vector only at a layout
/**
 * This function scales the values of one vector and copies it to another
 * only at a layout. No communication is performed.
 *
 * \param[in,out]	pVecDest		Destination vector of scaled copy
 * \param[in]		pVecSrc			source vector
 * \param[in]		scale			Scaling factor
 * \param[in]		Layout			Index Layout
 */
template <typename TVector>
void VecScaledCopyOnLayout(	TVector* pVecDest,
							const TVector* pVecSrc, number scale,
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

		// 	get entry
			typename TVector::value_type entry = (*pVecSrc)[index];

		//	scale entry
			entry *= scale;

		//	set/copy value
			(*pVecDest)[index] = entry;
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


//	returns the highest referenced index of the elements in the layout.
int GetHighestReferencedIndex(IndexLayout& layout);

////////////////////////////////////////////////////////////////////////
///	fills a connection list, which gives the connected processes to each entry.
/**	the first entry for each connection is the process on which the
 * master-element lies, followed by the processes where associated slaves lie.
 */
void CommunicateConnections(std::vector<std::vector<int> >& connectionsOut,
							IndexLayout& masterLayout,
							IndexLayout& slaveLayout,
							int highestReferencedIndexs);

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
		IndexLayout& crossPointMastersOut, IndexLayout& crossPointSlavesOut,
		IndexLayout& standardMasters, IndexLayout& standardSlaves,
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
 * to consistent on the diagonal. A ParallelCommunicator is created iff no communicator passed.
 *
 * \param[in,out]		pVec			Parallel Vector
 * \param[in]			masterLayout	Master Layout
 * \param[in]			slaveLayout		Slave Layout
 * \param[in]			pCom			Parallel Communicator
 */
template <typename TAlgebra>
void MatAdditiveToConsistentOnDiag(	typename TAlgebra::matrix_type* pMat,
                                   	IndexLayout& masterLayout, IndexLayout& slaveLayout,
                                   	pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
{
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


}//	end of namespace


#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__ */
