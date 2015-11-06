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

///	sets a matrix row corresponding to averaging the constrained index
template <typename TMatrix>
void SetInterpolation(TMatrix& A,
                      std::vector<size_t> & constrainedIndex,
                      std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	typedef typename TMatrix::row_iterator row_iterator;

	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	loop all constrained dofs
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	remove all couplings
		{
		const row_iterator iterEnd = A.end_row(constrainedIndex[i]);
		for(row_iterator conn = A.begin_row(constrainedIndex[i]); conn != iterEnd; ++conn)
			conn.value() = 0.0;
		}

	//	set diag of row to identity
		A(constrainedIndex[i], constrainedIndex[i]) = 1.0;

	//	set coupling to all contraining dofs the inverse of the
	//	number of contraining dofs
		const number frac = -1.0/(vConstrainingIndex.size());
		for(size_t j=0; j < vConstrainingIndex.size();++j)
			A(constrainedIndex[i], vConstrainingIndex[j][i]) = frac;
	}
}

template <typename TVector>
void InterpolateValues(TVector& u,
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
	//	get value to be interpolated
		block_type& val = u[constrainedIndex[i]];

	//	reset value
		val = 0.0;

	// 	add equally from all constraining indices
		for(size_t j=0; j < vConstrainingIndex.size(); ++j)
			VecScaleAdd(val, 1.0, val, frac, u[vConstrainingIndex[j][i]]);
	}
}



template <typename TMatrix>
void SplitAddRow_Symmetric(TMatrix& A,
                           std::vector<size_t> & constrainedIndex,
                           std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	typedef typename TMatrix::value_type block_type;
	typedef typename TMatrix::row_iterator row_iterator;

	UG_ASSERT(!vConstrainingIndex.empty(), "There have to be constraining indices!");

	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	scaling factor
	const number frac = 1./(vConstrainingIndex.size());

//	handle each contrained index
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	add coupling constrained dof -> constrained dof
	//	we don't have to adjust the block itself, since the row of
	//	constraints will be set to interpolation afterwards
		block_type block = A(constrainedIndex[i], constrainedIndex[i]);

	//	scale by weight
		block *= frac*frac;

	//	add coupling
		for(size_t k = 0; k < vConstrainingIndex.size(); ++k)
			for(size_t m = 0; m < vConstrainingIndex.size(); ++m)
		{
			A(vConstrainingIndex[k][i],
			  vConstrainingIndex[m][i]) += block;
		}

	//	loop coupling between constrained dof -> constraining dof
		for(row_iterator conn = A.begin_row(constrainedIndex[i]);
							conn != A.end_row(constrainedIndex[i]); ++conn)
		{
		//	skip self-coupling (already handled)
			const size_t j = conn.index();
			if(j == constrainedIndex[i]) continue;

		//	get coupling entry
			block_type block = conn.value();
			block_type blockT = A(j, constrainedIndex[i]);

		//	multiply the cpl value by the inverse number of constraining
			block *= frac;
			blockT *= frac;

		//	add the coupling to the constraining indices rows
			for(size_t k = 0; k < vConstrainingIndex.size(); ++k)
			{
				UG_ASSERT(vConstrainingIndex[k][i]!=constrainedIndex[i],
						"Modifying 'this' (=conn referenced) matrix row is invalid!" << constrainedIndex[i]);
				A(vConstrainingIndex[k][i], j) += block;
				A(j, vConstrainingIndex[k][i]) += blockT;
			}

		//	set the splitted coupling to zero
		//	this must only be done in columns, since the row associated to
		//	the contrained index will be set to an interpolation.
			A(j, constrainedIndex[i]) = 0.0;
		}
	}
}

template <typename TMatrix>
void SplitAddRow_OneSide(TMatrix& A,
                         std::vector<size_t> & constrainedIndex,
                         std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	typedef typename TMatrix::row_iterator row_iterator;
	typedef typename TMatrix::value_type block_type;

	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	scaling factor
	const number frac = 1./(vConstrainingIndex.size());

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	// choose randomly the first dof to add whole row
		const size_t addTo = vConstrainingIndex[0][i];

		block_type block = A(constrainedIndex[i], constrainedIndex[i]);

	//	scale by weight
		block *= frac;

	//	add coupling
		for(size_t k = 0; k < vConstrainingIndex.size(); ++k)
			A(addTo, vConstrainingIndex[k][i]) += block;

		for(row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
		//	skip self-coupling (already handled)
			const size_t j = conn.index();
			if(j == constrainedIndex[i]) continue;

		//	get transposed coupling entry
			block_type blockT = A(j, constrainedIndex[i]);
			blockT *= frac;

		//	add the coupling to the constraining indices rows
			for(size_t k = 0; k < vConstrainingIndex.size(); ++k)
				A(j, vConstrainingIndex[k][i]) += blockT;

		//	coupling due to one side adding
			block_type block = conn.value();
			A(addTo, j) += block;

		//	set the splitted coupling to zero
		//	this must only be done in columns, since the row associated to
		//	the contrained index will be set to an interpolation.
			A(j, constrainedIndex[i]) = 0.0;
		}
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

// 	get algebra indices for constrained and constraining vertices
	for (size_t fct = 0; fct < dd->num_fct(si); fct++)
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
			if (!dd->is_def_in_subset(fct, siC))
			{
				UG_THROW("Function " << fct << " is defined for a constrained vertex, "
						 "but not for one of its constraining vertices!");
			}

			dd->inner_algebra_indices_for_fct(vConstrainingVrt[i], vConstrainingInd[i], false, fct);
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
              ConstSmartPtr<DoFDistribution> dd, number time,
              ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			  const std::vector<number>* vScaleMass,
			  const std::vector<number>* vScaleStiff)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
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
           ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
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
	}
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd, number time,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const number s_a0)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
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

	// 	Split using indices
		SplitAddRow_Symmetric(J, constrainedInd, vConstrainingInd);

	//	set interpolation
		SetInterpolation(J, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
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

	// 	Split using indices
		SplitAddRow_Symmetric(mat, constrainedInd, vConstrainingInd);

	//	set interpolation
		SetInterpolation(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		SplitAddRhs_Symmetric(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
                number time)
{
	if(this->m_spAssTuner->single_index_assembling_enabled())
		UG_THROW("index-wise assemble routine is not "
				"implemented for SymP1Constraints \n");

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





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	OneSide P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<int dim>
struct SortVertexPos {

		SortVertexPos(SmartPtr<Domain<dim, MultiGrid, MGSubsetHandler> > spDomain)
			: m_aaPos(spDomain->position_accessor())
		{}

		inline bool operator() (Vertex* vrt1, Vertex* vrt2)
			{UG_THROW(dim <<" not implemented.");}

	protected:
  	  typename Domain<dim, MultiGrid, MGSubsetHandler>::position_accessor_type& m_aaPos;
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
			if (!dd->is_def_in_subset(fct, siC))
			{
				UG_THROW("Function " << fct << " is defined for a constrained vertex, "
						 "but not for one of its constraining vertices!");
			}

			dd->inner_algebra_indices_for_fct(vConstrainingVrt[i], vConstrainingInd[i], false, fct);
		}
	}
}


template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, number time,
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
           ConstSmartPtr<DoFDistribution> dd, number time)
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
	}
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd, number time,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const number s_a0)
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

	// 	Split using indices
		SplitAddRow_OneSide(J, constrainedInd, vConstrainingInd);

	//	set interpolation
		SetInterpolation(J, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<DoFDistribution> dd, number time)
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

//	loop constraining edges
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

	// 	Split using indices
		SplitAddRow_OneSide(mat, constrainedInd, vConstrainingInd);

	//	Set interpolation
		SetInterpolation(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		SplitAddRhs_OneSide(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
                number time)
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


}; // namespace ug



#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__ */
