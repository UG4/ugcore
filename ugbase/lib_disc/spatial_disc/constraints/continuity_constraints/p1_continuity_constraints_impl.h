/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__

#include "lib_disc/spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.h"
#include "lib_disc/domain.h"

namespace ug {

template <typename TMatrix>
static void setConstrainedRow
(
	TMatrix& J,
	size_t constrdInd,
	const std::vector<size_t>& vConstrgInd
)
{
	// set a unity row
	SetRow(J, constrdInd, 0.0);
	J(constrdInd, constrdInd) = 1.0;
}


template <typename TVector>
void InterpolateValues(TVector& u,
                       std::vector<size_t>& constrainedIndex,
                       std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	typedef typename TVector::value_type block_type;

	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	scaling factor
	const number frac = 1./(vConstrainingIndex.size());

//	loop constrained indices
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	get value to be interpolated
		block_type& val = u[constrainedIndex[i]];

	//	reset value
		val = 0.0;

	// 	add equally from all constraining indices
		for(size_t j=0; j < vConstrainingIndex.size(); ++j)
			VecScaleAdd(val, 1.0, val, frac, u[vConstrainingIndex[j][i]]);
			//val += frac * u[vConstrainingIndex[j][i]];
	}
}



template <typename TVector>
void SplitAddRhs_Symmetric(TVector& rhs,
                         std::vector<size_t> & constrainedIndex,
                         std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	typedef typename TVector::value_type block_type;

	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	scaling factor
	const number frac = 1./(vConstrainingIndex.size());

//	loop constrained indices
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	get constrained rhs
	//	modify block directly since set to zero afterwards
		block_type& val = rhs[constrainedIndex[i]];
		val *= frac;

	// 	split equally on all constraining indices
		for(size_t j=0; j < vConstrainingIndex.size(); ++j)
			rhs[vConstrainingIndex[j][i]] += val;

	//	set rhs to zero for constrained index
		val = 0.0;
	}
}

template <typename TVector>
void SplitAddRhs_OneSide(TVector& rhs,
                       std::vector<size_t> & constrainedIndex,
                       std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	typedef typename TVector::value_type block_type;

	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		// choose randomly the first dof to add whole
		// (must be the same as for row)
		const size_t addTo = vConstrainingIndex[0][i];

		block_type& val = rhs[constrainedIndex[i]];
		rhs[addTo] += val;
		val = 0.0;
	}
}


/**
 * @brief Extract algebra indices for constrained and constraining indices from DoF distribution
 *
 * @param dd                DoF distribution
 * @param hgVrt             the hanging vertex
 * @param vConstrainingVrt  vector of constraining vertices
 * @param constrainedInd    vector of algebra indices on hanging vertex
 * @param vConstrainingInd  vector of (vector of constraining algebra indices) for constraining vertices
 */
inline void get_algebra_indices(ConstSmartPtr<DoFDistribution> dd,
						 ConstrainedVertex* hgVrt,
						 std::vector<Vertex*>& vConstrainingVrt,
						 std::vector<size_t>& constrainedInd,
						 std::vector<std::vector<size_t> >& vConstrainingInd)
{
// get subset index
	const int si = dd->subset_handler()->get_subset_index(hgVrt);

//	get constraining vertices
	CollectConstraining(vConstrainingVrt, *dd->multi_grid(), hgVrt);

// clear constrainedInd
	constrainedInd.clear();

//	resize constraining indices
	vConstrainingInd.clear();
	vConstrainingInd.resize(vConstrainingVrt.size());

	if (dd->grouped())  // block algebra (assuming constrainers have exactly the same functions as constrained)
	{
	//	get indices for constrained vertex
		dd->inner_algebra_indices(hgVrt, constrainedInd, false);

	//	get indices for constraining vertices
		for (size_t i = 0; i < vConstrainingVrt.size(); ++i)
			dd->inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i], false);
	}
	else  // scalar algebra
	{
	// 	get algebra indices for constrained and constraining vertices
		for (size_t fct = 0; fct < dd->num_fct(); ++fct)
		{
		//	check that function is defined on subset
			if (!dd->is_def_in_subset(fct, si)) continue;

		//	get indices for constrained vertex
			dd->inner_algebra_indices_for_fct(hgVrt, constrainedInd, false, fct);

		//	get indices for constraining vertices
			for (size_t i = 0; i < vConstrainingVrt.size(); ++i)
			{
				const int siC = dd->subset_handler()->get_subset_index(vConstrainingVrt[i]);

			//	check that function is defined on subset
				UG_COND_THROW(!dd->is_def_in_subset(fct, siC),
					"Function " << fct << " is defined for a constrained vertex, "
					"but not for one of its constraining vertices!");

				dd->inner_algebra_indices_for_fct(vConstrainingVrt[i], vConstrainingInd[i], false, fct);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	Sym P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, int type, number time,
              ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			  const std::vector<number>* vScaleMass,
			  const std::vector<number>* vScaleStiff)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t> constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

	//	adapt rhs
		SplitAddRhs_Symmetric(d, constrainedInd, vConstrainingInd);
	}
}


template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t> constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

	//	adapt rhs
		SplitAddRhs_Symmetric(rhs, constrainedInd, vConstrainingInd);

		const size_t nCstrd = constrainedInd.size();
		for (size_t i = 0; i < nCstrd; ++i)
			rhs[constrainedInd[i]] = u[constrainedInd[i]];
	}
}

static void fillConstraintMapSymmetric
(
	std::map<size_t, std::vector<size_t> >& constraintMap,
	ConstSmartPtr<DoFDistribution> dd
)
{
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t> constrainedInd;
	std::vector<size_t> constrainers;
	std::vector<Vertex*> vConstrainingVrt;

	// get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

	// loop constrained vertices
	for (; iter != iterEnd; ++iter)
	{
		// get hanging vert
		ConstrainedVertex* hgVrt = *iter;

		// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

		// save in constraint map
		const size_t nInd = constrainedInd.size();
		const size_t nConstrainers = vConstrainingInd.size();
		constrainers.resize(nConstrainers);
		for (size_t i = 0; i < nInd; ++i)
		{
			for (size_t j = 0; j < nConstrainers; ++j)
				constrainers[j] = vConstrainingInd[j][i];
			constraintMap[constrainedInd[i]] = constrainers;
		}
	}
}

template <typename TMatrix>
static void splitConstrainedRowSymmetric
(
	TMatrix& J,
	size_t constrdInd,
	const std::vector<size_t>& vConstrgInd
)
{
	typedef typename TMatrix::value_type block_type;
	typedef typename TMatrix::row_iterator row_iterator;

	const size_t nConstrg = vConstrgInd.size();
	const number frac = 1.0 / nConstrg;

	// distribute row equally to constrainers
	for (row_iterator conn = J.begin_row(constrdInd); conn != J.end_row(constrdInd); ++conn)
	{
		const size_t connInd = conn.index();
		block_type block = conn.value();
		block *= frac;
		for (size_t k = 0; k < nConstrg; ++k)
			J(vConstrgInd[k], connInd) += block;
	}
}

template <typename TVector>
static void splitConstrainedRhsEntrySymmetric
(
	TVector& rhs,
	size_t constrdInd,
	const std::vector<size_t>& vConstrgInd
)
{
	typedef typename TVector::value_type block_type;

	const size_t nConstrg = vConstrgInd.size();
	const number frac = 1.0 / nConstrg;

	// get constrained rhs
	// modify block directly since set to zero afterwards
	block_type& val = rhs[constrdInd];
	val *= frac;

	// split equally on all constraining indices
	for (size_t i = 0; i < nConstrg; ++i)
		rhs[vConstrgInd[i]] += val;

	// set rhs to zero for constrained index
	val = 0.0;
}

template <typename TMatrix>
static void splitConstrainedCols
(
	TMatrix& J,
	const std::map<size_t, std::vector<size_t> >& constraintMap
)
{
	typedef typename TMatrix::value_type block_type;
	typedef typename TMatrix::row_iterator row_iterator;
	typedef typename std::map<size_t, std::vector<size_t> >::const_iterator constraint_iter;

	const size_t nRows = J.num_rows();
	std::vector<constraint_iter> tmpConstraints;
	for (size_t r = 0; r < nRows; ++r)
	{
		// leave out constrained rows (already final)
		if (constraintMap.find(r) != constraintMap.end())
			continue;

		// loop row
		tmpConstraints.clear();
		for (row_iterator conn = J.begin_row(r); conn != J.end_row(r); ++conn)
		{
			const size_t connInd = conn.index();

			// remember constrained columns
			constraint_iter it = constraintMap.find(connInd);
			if (it != constraintMap.end())
				tmpConstraints.push_back(it);
		}

		// split constrained columns to constrainers (separate iteration to not throw off iterators)
		const size_t nConstraints = tmpConstraints.size();
		for (size_t i = 0; i < nConstraints; ++i)
		{
			const size_t cdi = tmpConstraints[i]->first;
			const std::vector<size_t>& vConstrgInd = tmpConstraints[i]->second;
			const size_t nConstrg = vConstrgInd.size();
			const number frac = 1.0 / nConstrg;

			// add to constraining column entries
			block_type block = J(r, cdi);
			block *= frac;
			for (size_t k = 0; k < nConstrg; ++k)
				J(r, vConstrgInd[k]) += block;

			// set constrained column entry zero
			J(r, cdi) = 0.0;
		}
	}
}



template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd, int type, number time,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const number s_a0)
{
	if (this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

	std::map<size_t, std::vector<size_t> > constraintMap;
	typedef typename std::map<size_t, std::vector<size_t> >::const_iterator constraint_iter;
	fillConstraintMapSymmetric(constraintMap, dd);

	// split constrained rows to constrainers and set interpolation row (or identity)
	constraint_iter it = constraintMap.begin();
	constraint_iter itEnd = constraintMap.end();
	for (; it != itEnd; ++it)
	{
		size_t constrdInd = it->first;
		const std::vector<size_t>& vConstrgInd = it->second;

		// split original row
		splitConstrainedRowSymmetric(J, constrdInd, vConstrgInd);

		// set identity row (or interpolating row in linear case)
		setConstrainedRow(J, constrdInd, vConstrgInd);
	}

	// apply chain rule, i.e., also split constrained columns
	// (unfortunately, we have to iterate the whole matrix for this)
	splitConstrainedCols(J, constraintMap);
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	if (this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

	std::map<size_t, std::vector<size_t> > constraintMap;
	typedef typename std::map<size_t, std::vector<size_t> >::const_iterator constraint_iter;
	fillConstraintMapSymmetric(constraintMap, dd);

	// split constrained rows to constrainers and set interpolation row (or identity)
	constraint_iter it = constraintMap.begin();
	constraint_iter itEnd = constraintMap.end();
	for (; it != itEnd; ++it)
	{
		size_t constrdInd = it->first;
		const std::vector<size_t>& vConstrgInd = it->second;

		// split original row
		splitConstrainedRowSymmetric(mat, constrdInd, vConstrgInd);

		// split original rhs entry
		splitConstrainedRhsEntrySymmetric(rhs, constrdInd, vConstrgInd);
		rhs[constrdInd] = u[constrdInd];

		// set identity row (or interpolating row in linear case)
		setConstrainedRow(mat, constrdInd, vConstrgInd);
	}

	// apply chain rule, i.e., also split constrained columns
	// (unfortunately, we have to iterate the whole matrix for this)
	splitConstrainedCols(mat, constraintMap);
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
				int type, number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t> constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constraining edges
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

	// 	Interpolate values
		InterpolateValues(u, constrainedInd, vConstrainingInd);
	}
}


template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_prolongation
(
	matrix_type& P,
	ConstSmartPtr<DoFDistribution> ddFine,
	ConstSmartPtr<DoFDistribution> ddCoarse,
	int type,
	number time
)
{
	if (this->m_spAssTuner->single_index_assembling_enabled())
			UG_THROW("index-wise assemble routine is not "
					"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t> constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = ddFine->begin<ConstrainedVertex>();
	iterEnd = ddFine->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
		get_algebra_indices(ddFine, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

	//	set zero row
		size_t sz = constrainedInd.size();
		for (size_t i = 0; i < sz; ++i)
			SetRow(P, constrainedInd[i], 0.0);
	}
}


template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_restriction
(
	matrix_type& R,
	ConstSmartPtr<DoFDistribution> ddCoarse,
	ConstSmartPtr<DoFDistribution> ddFine,
	int type,
	number time
)
{

}



template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_correction
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	//typedef typename vector_type::value_type block_type;

	if (this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints \n");

	// storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

	// get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

	// loop constrained vertices
	for (; iter != iterEnd; ++iter)
	{
		// get hanging vert
		ConstrainedVertex* hgVrt = *iter;

		// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

		// set all entries corresponding to constrained dofs to zero
		for (size_t i = 0; i < constrainedInd.size(); ++i)
			u[constrainedInd[i]] = 0.0;
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	OneSide P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<int dim>
struct SortVertexPos {

		SortVertexPos(ConstSmartPtr<Domain<dim, MultiGrid, MGSubsetHandler> > spDomain)
			: m_aaPos(spDomain->position_accessor())
		{}

		inline bool operator() (Vertex* vrt1, Vertex* vrt2)
			{UG_THROW(dim <<" not implemented.");}

	protected:
		const typename Domain<dim, MultiGrid, MGSubsetHandler>::position_accessor_type& m_aaPos;
};

template<>
inline bool SortVertexPos<1>::operator() (Vertex* vrt1, Vertex* vrt2)
{
	if(m_aaPos[vrt1][0] < m_aaPos[vrt2][0]) {
		return true;
	}
	return false;
}

template<>
inline bool SortVertexPos<2>::operator() (Vertex* vrt1, Vertex* vrt2)
{
	if(m_aaPos[vrt1][0] < m_aaPos[vrt2][0]) {
		return true;
	}
	else if(m_aaPos[vrt1][0] == m_aaPos[vrt2][0]) {
		if(m_aaPos[vrt1][1] < m_aaPos[vrt2][1])
			return true;
	}
	return false;
}

template<>
inline bool SortVertexPos<3>::operator() (Vertex* vrt1, Vertex* vrt2)
{
	if(m_aaPos[vrt1][0] < m_aaPos[vrt2][0]) {
		return true;
	}
	else if(m_aaPos[vrt1][0] == m_aaPos[vrt2][0]) {
		if(m_aaPos[vrt1][1] < m_aaPos[vrt2][1]){
			return true;
		}
		else if(m_aaPos[vrt1][1] == m_aaPos[vrt2][1]){
			if(m_aaPos[vrt1][2] < m_aaPos[vrt2][2])
				return true;
		}
		return false;
	}
	return false;
}


/**
 * @brief Extract DoF indices for constrained and constraining indices from DoF distribution
 *
 * One cannot simply use algebra indices as constrainers and constrained vertices
 * might not have the same number of functions defined on them. Mapping correct
 * indices is only possible through DoF indices.
 *
 * @param dd                DoF distribution
 * @param hgVrt             the hanging vertex
 * @param vConstrainingVrt  vector of constraining vertices
 * @param constrainedInd    vector of DoFs indices on hanging vertex
 * @param vConstrainingInd  vector of (vector of constraining DoF indices) for constraining vertices
 * @param sortVertexPos     sorting functional for constrainers
 */
template <typename TDomain>
inline void get_algebra_indices(ConstSmartPtr<DoFDistribution> dd,
						 ConstrainedVertex* hgVrt,
						 std::vector<Vertex*>& vConstrainingVrt,
						 std::vector<size_t>& constrainedInd,
						 std::vector<std::vector<size_t> >& vConstrainingInd,
						 const SortVertexPos<TDomain::dim>& sortVertexPos)
{
// get subset index
	const int si = dd->subset_handler()->get_subset_index(hgVrt);

//	get constraining vertices
	CollectConstraining(vConstrainingVrt, *dd->multi_grid(), hgVrt);

#ifdef UG_PARALLEL
	std::sort(vConstrainingVrt.begin(), vConstrainingVrt.end(), sortVertexPos);
#endif

// clear constrainedInd
	constrainedInd.clear();

//	resize constraining indices
	vConstrainingInd.clear();
	vConstrainingInd.resize(vConstrainingVrt.size());

// 	get algebra indices for constrained and constraining vertices
	if (dd->grouped())  // block algebra (assuming constrainers have exactly the same functions as constrained)
	{
	//	get indices for constrained vertex
		dd->inner_algebra_indices(hgVrt, constrainedInd, false);

	//	get indices for constraining vertices
		for (size_t i = 0; i < vConstrainingVrt.size(); ++i)
			dd->inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i], false);
	}
	else  // scalar algebra
	{
		for (size_t fct = 0; fct < dd->num_fct(); fct++)
		{
		//	check that function is defined on subset
			if (!dd->is_def_in_subset(fct, si)) continue;

		//	get indices for constrained vertex
			dd->inner_algebra_indices_for_fct(hgVrt, constrainedInd, false, fct);

		//	get indices for constraining vertices
			for (size_t i = 0; i < vConstrainingVrt.size(); ++i)
			{
				const int siC = dd->subset_handler()->get_subset_index(vConstrainingVrt[i]);

			//	check that function is defined on subset
				UG_COND_THROW(!dd->is_def_in_subset(fct, siC),
					"Function " << fct << " is defined for a constrained vertex, "
					"but not for one of its constraining vertices!");

				dd->inner_algebra_indices_for_fct(vConstrainingVrt[i], vConstrainingInd[i], false, fct);
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, int type, number time,
              ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
     		  const std::vector<number>* vScaleMass,
              const std::vector<number>* vScaleStiff)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

#ifdef UG_PARALLEL
	SortVertexPos<TDomain::dim> sortVertexPos(this->approximation_space()->domain());
#endif

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
#ifdef UG_PARALLEL
		get_algebra_indices<TDomain>(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd, sortVertexPos);
#else
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);
#endif

	//	adapt rhs
		SplitAddRhs_OneSide(d, constrainedInd, vConstrainingInd);
	}
}


template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

#ifdef UG_PARALLEL
	SortVertexPos<TDomain::dim> sortVertexPos(this->approximation_space()->domain());
#endif

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
#ifdef UG_PARALLEL
		get_algebra_indices<TDomain>(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd, sortVertexPos);
#else
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);
#endif

	//	adapt rhs
		SplitAddRhs_OneSide(rhs, constrainedInd, vConstrainingInd);

		const size_t nCstrd = constrainedInd.size();
		for (size_t i = 0; i < nCstrd; ++i)
			rhs[constrainedInd[i]] = u[constrainedInd[i]];
	}
}


template <typename TDomain>
static void fillConstraintMapOneSide
(
	std::map<size_t, std::vector<size_t> >& constraintMap,
	ConstSmartPtr<DoFDistribution> dd,
	ConstSmartPtr<TDomain> dom
)
{
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t> constrainedInd;
	std::vector<size_t> constrainers;
	std::vector<Vertex*> vConstrainingVrt;

#ifdef UG_PARALLEL
	SortVertexPos<TDomain::dim> sortVertexPos(dom);
#endif

	// get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

	// loop constrained vertices
	for (; iter != iterEnd; ++iter)
	{
		// get hanging vert
		ConstrainedVertex* hgVrt = *iter;

		// get algebra indices for constrained and constraining vertices
#ifdef UG_PARALLEL
		get_algebra_indices<TDomain>(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd, sortVertexPos);
#else
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);
#endif

		// save in constraint map
		const size_t nInd = constrainedInd.size();
		const size_t nConstrainers = vConstrainingInd.size();
		constrainers.resize(nConstrainers);
		for (size_t i = 0; i < nInd; ++i)
		{
			for (size_t j = 0; j < nConstrainers; ++j)
				constrainers[j] = vConstrainingInd[j][i];
			constraintMap[constrainedInd[i]] = constrainers;
		}
	}
}

template <typename TMatrix>
static void splitConstrainedRowOneSide
(
	TMatrix& J,
	size_t constrdInd,
	const std::vector<size_t>& vConstrgInd
)
{
	typedef typename TMatrix::row_iterator row_iterator;

	// add complete row to first constrainer
	for (row_iterator conn = J.begin_row(constrdInd); conn != J.end_row(constrdInd); ++conn)
	{
		const size_t connInd = conn.index();
		J(vConstrgInd[0], connInd) += conn.value();
	}
}

template <typename TVector>
static void splitConstrainedRhsEntryOneSide
(
	TVector& rhs,
	size_t constrdInd,
	const std::vector<size_t>& vConstrgInd
)
{
	// add complete entry to first constrainer index
	typename TVector::value_type& val = rhs[constrdInd];
	rhs[vConstrgInd[0]] += val;

	// set rhs to zero for constrained index
	val = 0.0;
}


template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd, int type, number time,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const number s_a0)
{
	if (this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints.");

	std::map<size_t, std::vector<size_t> > constraintMap;
	typedef typename std::map<size_t, std::vector<size_t> >::const_iterator constraint_iter;
	fillConstraintMapOneSide<TDomain>(constraintMap, dd, this->approximation_space()->domain());

	// split constrained rows to constrainers and set interpolation row (or identity)
	constraint_iter it = constraintMap.begin();
	constraint_iter itEnd = constraintMap.end();
	for (; it != itEnd; ++it)
	{
		size_t constrdInd = it->first;
		const std::vector<size_t>& vConstrgInd = it->second;

		// split original row
		splitConstrainedRowOneSide(J, constrdInd, vConstrgInd);

		// set identity row (or interpolating row in linear case)
		setConstrainedRow(J, constrdInd, vConstrgInd);
	}

	// apply chain rule, i.e., also split constrained columns
	// (unfortunately, we have to iterate the whole matrix for this)
	splitConstrainedCols(J, constraintMap);
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, int type, number time)
{
	if (this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints");

	std::map<size_t, std::vector<size_t> > constraintMap;
	typedef typename std::map<size_t, std::vector<size_t> >::const_iterator constraint_iter;
	fillConstraintMapOneSide<TDomain>(constraintMap, dd, this->approximation_space()->domain());

	// split constrained rows to constrainers and set interpolation row (or identity)
	constraint_iter it = constraintMap.begin();
	constraint_iter itEnd = constraintMap.end();
	for (; it != itEnd; ++it)
	{
		size_t constrdInd = it->first;
		const std::vector<size_t>& vConstrgInd = it->second;

		// split original row
		splitConstrainedRowOneSide(mat, constrdInd, vConstrgInd);

		// split original rhs entry
		splitConstrainedRhsEntryOneSide(rhs, constrdInd, vConstrgInd);
		rhs[constrdInd] = u[constrdInd];

		// set identity row (or interpolating row in linear case)
		setConstrainedRow(mat, constrdInd, vConstrgInd);
	}

	// apply chain rule, i.e., also split constrained columns
	// (unfortunately, we have to iterate the whole matrix for this)
	splitConstrainedCols(mat, constraintMap);
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
				int type, number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constraining edges
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

	// 	Interpolate values
		InterpolateValues(u, constrainedInd, vConstrainingInd);
	}
}



template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_prolongation
(
	matrix_type& P,
	ConstSmartPtr<DoFDistribution> ddFine,
	ConstSmartPtr<DoFDistribution> ddCoarse,
	int type,
	number time
)
{
	if (this->m_spAssTuner->single_index_assembling_enabled())
			UG_THROW("index-wise assemble routine is not "
					"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = ddFine->begin<ConstrainedVertex>();
	iterEnd = ddFine->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

	// get algebra indices for constrained and constraining vertices
		get_algebra_indices(ddFine, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

	//	set zero row
		size_t sz = constrainedInd.size();
		for (size_t i = 0; i < sz; ++i)
			SetRow(P, constrainedInd[i], 0.0);
	}
}


template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_restriction
(
	matrix_type& R,
	ConstSmartPtr<DoFDistribution> ddCoarse,
	ConstSmartPtr<DoFDistribution> ddFine,
	int type,
	number time
)
{}




template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_correction
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	//typedef typename vector_type::value_type block_type;

	if (this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for OneSideP1Constraints \n");

	// storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<Vertex*> vConstrainingVrt;

	// get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

	// loop constrained vertices
	for (; iter != iterEnd; ++iter)
	{
		// get hanging vert
		ConstrainedVertex* hgVrt = *iter;

		// get algebra indices for constrained and constraining vertices
		get_algebra_indices(dd, hgVrt, vConstrainingVrt, constrainedInd, vConstrainingInd);

		// set all entries corresponding to constrained dofs to zero
		for (size_t i = 0; i < constrainedInd.size(); ++i)
			u[constrainedInd[i]] = 0.0;
	}
}


}; // namespace ug



#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__ */
