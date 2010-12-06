/*
 * parallelization.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLELIZATION_UTIL__

#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "lib_discretization/lib_discretization.h"
#include "lib_grid/parallelization/parallelization.h"
#include "lib_discretization/dof_manager/dof_distribution.h"

namespace ug
{

///	Adds dof-indices of elements in elemLayout to the specified IndexLayout.
/**	Make sure that TLayout holds elements of type VertexBase*, EdgeBase*,
 *  Face* or Volume*.
 *
 *  \todo: replace IndexLayout with TDoFManager::IndexLayout.
 */
template <class TDoFDistr, class TLayout>
bool AddEntriesToIndexLayout(IndexLayout& indexLayoutOut,
							TDoFDistr& dofDistr,
							TLayout& elemLayout)
{
	typedef typename TLayout::iterator InterfaceIterator;
	typedef typename TLayout::Interface ElemInterface;
	typedef typename ElemInterface::iterator ElemIterator;

	typedef IndexLayout::Interface IndexInterface;

//	iterate over all interfaces
	for(InterfaceIterator iIter = elemLayout.begin();
		iIter != elemLayout.end(); ++iIter)
	{
		ElemInterface& elemInterface = elemLayout.interface(iIter);
		IndexInterface& indexInterface = indexLayoutOut.interface(
											elemLayout.proc_id(iIter));

	//	iterate over entries in the elemInterface and add associated
	//	dofs to the indexInterface
		for(ElemIterator eIter = elemInterface.begin();
			eIter != elemInterface.end(); ++eIter)
		{
			typename ElemInterface::Element elem = elemInterface.get_element(eIter);
			typename TDoFDistr::algebra_index_vector_type indices;
			dofDistr.get_inner_algebra_indices(elem, indices);
			for(size_t i = 0; i < indices.size(); ++i)
			{
				indexInterface.push_back(indices[i]);
			}
		}
	}
	return true;
}


template <class TDoFDistribution>
bool CreateIndexLayout(	IndexLayout& layoutOut,
						TDoFDistribution& dofDistr,
						GridLayoutMap& layoutMap,
						int keyType, int level)
{
//TODO: clear the layout!
	bool bRetVal = true;
	if(layoutMap.has_layout<VertexBase>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofDistr,
								layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level));
	}
/*
	if(layoutMap.has_layout<EdgeBase>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofManager,
								layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level));
	}
	if(layoutMap.has_layout<Face>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofManager,
								layoutMap.get_layout<Face>(keyType).layout_on_level(level));
	}
	if(layoutMap.has_layout<Volume>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofManager,
								layoutMap.get_layout<Volume>(keyType).layout_on_level(level));
	}
*/
	return bRetVal;
}

template <typename TMatrix, typename TDoFDistr>
void CopyLayoutsAndCommunicatorIntoMatrix(TMatrix& mat, IDoFDistribution<TDoFDistr>& dofDistr)
{
	mat.set_slave_layout(dofDistr.get_slave_layout());
	mat.set_master_layout(dofDistr.get_master_layout());
	mat.set_communicator(dofDistr.get_communicator());
	mat.set_process_communicator(dofDistr.get_process_communicator());
}

}//	end of namespace

#endif
