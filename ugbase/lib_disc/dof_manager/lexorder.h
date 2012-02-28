/*
 * cuthill_mckee.h
 *
 *  Created on: 21.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__
#define __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__

#include <vector>
#include <utility> // for pair

#include "lib_disc/function_spaces/approximation_space.h"

namespace ug{

/*
template<int dim>
struct ComparePosition {
///	constructor, passing field with connections for each index
	ComparePosition(const std::vector<std::pair<MathVector<dim>, size_t> >& vPos) : m_vPos(vPos) {}

///	comparison operator
	bool operator() (size_t i,size_t j)
	{
		UG_ASSERT(i < m_vPos.size(), "Invalid index.");
		UG_ASSERT(j < m_vPos.size(), "Invalid index.");
		return ComparePosDim(m_vPos[i].first, m_vPos[j].first);
	}

private:
	///	storage field for teh position of each index
	const std::vector<MathVector<dim> >& m_vPos;


};*/

// Order for 1D
template<class MVT>
bool ComparePosDim(const std::pair<MVT, size_t> &p1,
		const std::pair<MVT, size_t> &p2)
	{return false;}

template<>
inline bool ComparePosDim(const std::pair<MathVector<1>, size_t> &p1,
                          const std::pair<MathVector<1>, size_t> &p2)
	{return (p1.first[0]<p2.first[0]);};

	// Order for 2D
template<>
inline bool ComparePosDim(const std::pair<MathVector<2>, size_t> &p1,
      	                   const std::pair<MathVector<2>, size_t> &p2)
	{ if (p1.first[1]==p2.first[1]) return p1.first[0] < p2.first[0];
	else return (p1.first[1] < p2.first[1]); };

	// Order for 3D
template<>
inline bool ComparePosDim(const std::pair<MathVector<3>, size_t> &p1,
                          const std::pair<MathVector<3>, size_t> &p2)
	{return (p1.first[0]==p2.first[0] && p1.first[1]==p2.first[1] && p1.first[2]<p2.first[2]);};

	// computes ordering using Cuthill-McKee algorithm

template<class MVT>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
					std::vector<std::pair<MVT, size_t> >& vPos,
                    int order, int sign)
{
//	list of sorted indices
	const size_t numVec = vPos.size();
	vNewIndex.resize(numVec);

// sort indices based on their position
//	ComparePosition<dim> myCompPos(vPos);
	std::sort(vPos.begin(), vPos.end(), ComparePosDim<MVT>);


	for (size_t i=0; i<numVec; ++i) {
		vNewIndex[vPos[i].second] = i;
		//UG_LOG("Mapping: (" << vPos[i].first << ", " << vPos[i].second << ") ->" << i << "\n");

	}
}




/// returns an array describing the needed index mapping for Cuthill-McKee ordering
/**
 * This function computes a index mapping, that transforms a index-graph into
 * Cuthill-McKee ordering. For each index an vector of all adjacent indices
 * must be passed. (If no adjacent index is passed for an index, this index is
 * skipped and not sorted). On exit the index field vNewIndex is filled with
 * the index mapping: newInd = vNewIndex[oldInd]
 *
 * \param[out]	vNewIndex		vector returning new index for old index
 * \param[in]	vvNeighbour		vector of adjacent indices for each index
 * \param[in]	bReverse		flag if "reverse Cuthill-McKee" is used
 * \returns		flag if ordering was successful
 */
//template<int dim>
//bool ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
 //            std::vector<std::pair<MathVector<dim>, size_t> >& vPos,
  //           int order, int sign);

///	writes positions of vertex dofs into a std::vector
template <typename TDD, typename TDomain>
void ExtractPos(	const TDD& dofDistr,
						typename TDomain::position_accessor_type& aaPos,
						std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPositions)
{
//	number of total dofs
	int nr = dofDistr.num_indices();

//	resize positions
	vPositions.resize(nr);

	typedef typename TDD::template traits<VertexBase>::const_iterator const_iterator;

//	loop all subsets
	for(int si = 0; si < dofDistr.num_subsets(); ++si)
	{
	//	loop all vertices
		const_iterator iter = dofDistr.template begin<VertexBase>(si);
		const_iterator iterEnd = dofDistr.template begin<VertexBase>(si);

		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			VertexBase* v = *iter;

		//	algebra indices vector
			std::vector<size_t> ind;

		//	load indices associated with vertex
			dofDistr.inner_algebra_indices(v, ind);

		//	write position
			for(size_t i = 0; i < ind.size(); ++i)
			{
				const size_t index = ind[i];
				vPositions[index].first = aaPos[v];
				vPositions[index].second = index;
			}
		}
	}
}

/// orders the dof distribution using Cuthill-McKee
template <typename TDD, typename TDomain>
void OrderLexForDofDist(TDD& dofDistr,
						typename TDomain::position_accessor_type& aaPos,
						int orderflag)
{
	//	get position attachment
	typedef MathVector<TDomain::dim> vec_type;
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> pos_type;

	std::vector<pos_type> vPositions;
	ExtractPos<TDD, TDomain>(dofDistr, aaPos, vPositions);

	//	get mapping: old -> new index
	std::vector<size_t> vNewIndex(dofDistr.num_indices());
	ComputeLexicographicOrder<vec_type>(vNewIndex, vPositions, orderflag, 0);

	//	reorder indices
	dofDistr.permute_indices(vNewIndex);
}


/// orders the all DofDistributions of the ApproximationSpace using lexicographic order
template <typename TDomain>
void OrderLex(ApproximationSpace<TDomain>& approxSpace,
              const char *order)
{
	// TODO: decode input
	int flag= 0;

	//	get position attachment
	typedef TDomain domain_type;
	typename domain_type::position_accessor_type& aaPos
			= approxSpace.domain()->position_accessor();

	//	order levels
	for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
		OrderLexForDofDist<LevelDoFDistribution, TDomain>(*approxSpace.level_dof_distribution(lev), aaPos, flag);

	//	order surface
	OrderLexForDofDist<SurfaceDoFDistribution, TDomain>(*approxSpace.surface_dof_distribution(), aaPos, flag);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__ */
