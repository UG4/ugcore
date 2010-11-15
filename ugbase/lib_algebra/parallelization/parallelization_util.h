/*
 * parallelization_util.h
 *
 *  Created on: 14.6.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__


namespace ug{

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

	//	perform communication on the level
		com.send_data(slaveLayout, cpVecAdd);
		com.receive_data(masterLayout, cpVecAdd);
		com.communicate();

	//	step 2: copy master values to slaves
	//	create the required communication policies
		ComPol_VecCopy<TVector> cpVecCopy(pVec);

	//	perform communication on the level
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
}

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

	//	perform communication on the level
		com.send_data(masterLayout, cpVecCopy);
		com.receive_data(slaveLayout, cpVecCopy);
		com.communicate();
}


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

	//	perform communication on the level
		com.send_data(slaveLayout, cpVecAddSetZero);
		com.receive_data(masterLayout, cpVecAddSetZero);
		com.communicate();
}


template <typename TVector>
void ConsistentToUnique(	TVector* pVec,
							IndexLayout& slaveLayout)
{
	typename IndexLayout::iterator iter = slaveLayout.begin();
	typename IndexLayout::iterator end = slaveLayout.end();

	for(; iter != end; ++iter)
	{
		typename IndexLayout::Interface& interface = slaveLayout.interface(iter);
		{
			for(typename IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				(*pVec)[index] = 0.0;
			}
		}
	}
}

/** subtracts values from one layout from the other (bidirectional)
 *
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


/** adds one vector to another only on the interface
 *
 */
template <typename TVector>
void VecAddOnLayoutWithoutCommunication(	TVector* pVecDest, const TVector* pVecSrc,
											IndexLayout& Layout)
{
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
		typename IndexLayout::Interface& interface = Layout.interface(iter);
		{
			for(typename IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				(*pVecDest)[index] += (*pVecSrc)[index];
			}
		}
	}
}


/** sets dirichlet rows for interface indices
 *
 */
template <typename TMatrix>
void MatSetDirichletWithoutCommunication(	TMatrix* pMatrix,
											IndexLayout& Layout)
{
	typename IndexLayout::iterator iter = Layout.begin();
	typename IndexLayout::iterator end = Layout.end();

	for(; iter != end; ++iter)
	{
		typename IndexLayout::Interface& interface = Layout.interface(iter);
		{
			for(typename IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				SetDirichletRow(*pMatrix, index);
			}
		}
	}
}

}//	end of namespace


#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL__ */
