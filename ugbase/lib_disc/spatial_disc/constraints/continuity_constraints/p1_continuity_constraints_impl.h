/*
 * p1_continuity_constraints_impl.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__

#include "lib_disc/spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.h"

namespace ug {

template <typename TMatrix>
void SplitAddRow_Symmetric(TMatrix& A,
                           std::vector<size_t> & constrainedIndex,
                           std::vector<std::vector<size_t> >& vConstrainingIndices)
{
//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.");

//	handle each contrained index
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	add coupling constrained dof -> constrained dof
	//	get entry
		typename TMatrix::value_type& block
				= A(constrainedIndex[i], constrainedIndex[i]);

	//	scale by weight
		block *= (1./(vConstrainingIndices.size()))*
					(1./(vConstrainingIndices.size()));

	//	add coupling
		for(size_t k = 0; k < vConstrainingIndices.size(); ++k)
			for(size_t m = 0; m < vConstrainingIndices.size(); ++m)
		{
			A(vConstrainingIndices[k][i],
			  vConstrainingIndices[m][i]) += block;
		}

	//	reset block
		 A(constrainedIndex[i], constrainedIndex[i]) = 0.0;

	//	loop coupling between constrained dof -> constraining dof
		for(typename TMatrix::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
		//	skip self-coupling (already handled)
			const size_t j = conn.index();
			if(j == constrainedIndex[i]) continue;

		//	get coupling entry
			typename TMatrix::value_type block = conn.value();

		//	get transposed coupling entry
			typename TMatrix::value_type blockT = A(j, constrainedIndex[i]);

		//	multiply the cpl value by the inverse number of constraining
		//	indices
			block *= 1./(vConstrainingIndices.size());
			blockT *= 1./(vConstrainingIndices.size());

		//	add the coupling to the constraining indices rows
			for(size_t k = 0; k < vConstrainingIndices.size(); ++k)
			{
				A(vConstrainingIndices[k][i], j) += block;

				A(j, vConstrainingIndices[k][i]) += blockT;
			}

		//	set the splitted coupling to zero
			conn.value() = 0.0;
			A(j, constrainedIndex[i]) = 0.0;
		}
	}
}

///	sets a matrix row corresponding to averaging the constrained index
template <typename TMatrix>
void SetInterpolation(TMatrix& A,
                      std::vector<size_t> & constrainedIndex,
                      std::vector<std::vector<size_t> >& vConstrainingIndex)
{
//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		if(vConstrainingIndex[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.\n");

//	loop all constrained dofs
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	remove all couplings
		for(typename TMatrix::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
			conn.value() = 0.0;
		}

	//	set diag of row to identity
		A(constrainedIndex[i], constrainedIndex[i]) = 1.0;

	//	set coupling to all contraining dofs the inverse of the
	//	number of contraining dofs
		number frac = -1.0/(vConstrainingIndex.size());
		for(size_t j=0; j < vConstrainingIndex.size();++j)
			A(constrainedIndex[i], vConstrainingIndex[j][i]) = frac;
	}
}

template <typename TVector>
void HandleRhs_Symmetric(TVector& rhs,
                         std::vector<size_t> & constrainedIndex,
                         std::vector<std::vector<size_t> >& vConstrainingIndices)
{
//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.");

//	loop constrained indices
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	get constrained rhs
		typename TVector::value_type& val = rhs[constrainedIndex[i]];
		val *= 1./(vConstrainingIndices.size());

	// 	split equally on all constraining indices
		for(size_t j=0; j < vConstrainingIndices.size(); ++j)
			rhs[vConstrainingIndices[j][i]] += val;

	//	set rhs to zero for contrained index
		val = 0.0;
	}
}

template <typename TVector>
void InterpolateValues_Symmetric(TVector& u,
                                 std::vector<size_t> & constrainedIndex,
                                 std::vector<std::vector<size_t> >& vConstrainingIndices)
{
//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.\n");

//	loop constrained indices
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	get constrained rhs
		typename TVector::value_type& val = u[constrainedIndex[i]];
		const number scale = 1./(vConstrainingIndices.size());

		val = 0.0;

	// 	split equally on all constraining indices
		for(size_t j=0; j < vConstrainingIndices.size(); ++j)
		{
			typename TVector::value_type entry = u[vConstrainingIndices[j][i]];
			entry *= scale;
			val += entry;
		}
	}
}


template <typename TMatrix>
void SplitAddRow_OneSide(TMatrix& A,
                         std::vector<size_t> & constrainedIndex,
                         std::vector<std::vector<size_t> >& vConstrainingIndices)
{
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wring number of indices. Cannot split row.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		for(typename TMatrix::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
			typename TMatrix::value_type block = conn.value();
			const size_t j = conn.index();

			// choose randomly the first dof to add whole row
			A(vConstrainingIndices[0][i], j) += block;
			A(constrainedIndex[i], j) = 0.0;
		}
	}
}

template <typename TVector>
void HandleRhs_OneSide(TVector& rhs,
                       std::vector<size_t> & constrainedIndex,
                       std::vector<std::vector<size_t> >& vConstrainingIndices)
{
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		typename TVector::value_type& val = rhs[constrainedIndex[i]];

		// choose randomly the first dof to add whole rhs (must be the same as for row)
		rhs[vConstrainingIndices[0][i]] += val;
		val = 0.0;
	}
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	Sym P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<TDD> dd, number time)
{
//	algebra indices of constraining vertex
	std::vector<std::vector<size_t> > vConstrainingInd;

//	algebra indices of constrained vertex
	std::vector<size_t>  constrainedInd;

//	vector of constraining vertices
	std::vector<VertexBase*> vConstrainingVrt;

//	iterators for hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;

//	get begin end of hanging vertices
	iter = dd->template begin<HangingVertex>();
	iterEnd = dd->template end<HangingVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		HangingVertex* hgVrt = *iter;

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

	// 	get algebra indices for constraining vertices
		for(size_t i=0; i < vConstrainingVrt.size(); ++i)
			dd->inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

	// 	get algebra indices constrained vertices
		dd->inner_algebra_indices(hgVrt, constrainedInd);

	// 	Split using indices
		SplitAddRow_Symmetric(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		HandleRhs_Symmetric(rhs, constrainedInd, vConstrainingInd);
	}

//	get begin end of hanging vertices
	iter = dd->template begin<HangingVertex>();
	iterEnd = dd->template end<HangingVertex>();

//	second loop to set the constraints
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		HangingVertex* hgVrt = *iter;

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

	// 	get algebra indices for constraining vertices
		for(size_t i=0; i < vConstrainingVrt.size(); ++i)
			dd->inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

	// 	get algebra indices constrained vertices
		dd->inner_algebra_indices(hgVrt, constrainedInd);

	//	set interpolation
		SetInterpolation(mat, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd,
                number time)
{
//	algebra indices of constraining vertex
	std::vector<std::vector<size_t> > vConstrainingInd;

//	algebra indices of constrained vertex
	std::vector<size_t>  constrainedInd;

//	vector of constraining vertices
	std::vector<VertexBase*> vConstrainingVrt;

//	iterators for hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;

//	get begin end of hanging vertices
	iter = dd->template begin<HangingVertex>();
	iterEnd = dd->template end<HangingVertex>();

//	loop constraining edges
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		HangingVertex* hgVrt = *iter;

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

	// 	get algebra indices for constraining vertices
		for(size_t i=0; i < vConstrainingVrt.size(); ++i)
			dd->inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

	// 	get algebra indices constrained vertices
		dd->inner_algebra_indices(hgVrt, constrainedInd);

	// 	Interpolate values
		InterpolateValues_Symmetric(u, constrainedInd, vConstrainingInd);
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	OneSide P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<TDD> dd, number time)
{
//	algebra indices of constraining vertex
	std::vector<std::vector<size_t> > vConstrainingInd;

//	algebra indices of constrained vertex
	std::vector<size_t>  constrainedInd;

//	vector of constraining vertices
	std::vector<VertexBase*> vConstrainingVrt;

//	iterators for hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;

//	get begin end of hanging vertices
	iter = dd->template begin<HangingVertex>();
	iterEnd = dd->template end<HangingVertex>();

//	loop constraining edges
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		HangingVertex* hgVrt = *iter;

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

	// 	get algebra indices for constraining vertices
		for(size_t i=0; i < vConstrainingVrt.size(); ++i)
			dd->inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

	// 	get algebra indices constrained vertices
		dd->inner_algebra_indices(hgVrt, constrainedInd);

	// 	Split using indices
		SplitAddRow_OneSide(mat, constrainedInd, vConstrainingInd);

	//	Set interpolation
		SetInterpolation(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		HandleRhs_OneSide(rhs, constrainedInd, vConstrainingInd);
	}
}

}; // namespace ug



#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__ */
