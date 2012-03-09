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

///	sets a matrix row corresponding to averaging the constrained index
template <typename TMatrix>
void SetInterpolation(TMatrix& A,
                      std::vector<size_t> & constrainedIndex,
                      std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	loop all constrained dofs
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	remove all couplings
		for(typename TMatrix::row_iterator conn = A.begin_row(constrainedIndex[i]);
									conn != A.end_row(constrainedIndex[i]); ++conn)
			conn.value() = 0.0;

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
	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	loop constrained indices
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	get value to be interpolated
		typename TVector::value_type& val = u[constrainedIndex[i]];

	//	scaling factor
		const number frac = 1./(vConstrainingIndex.size());

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
	//	get entry
		typename TMatrix::value_type& block
				= A(constrainedIndex[i], constrainedIndex[i]);

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
			block *= frac;
			blockT *= frac;

		//	add the coupling to the constraining indices rows
			for(size_t k = 0; k < vConstrainingIndex.size(); ++k)
			{
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
	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		for(typename TMatrix::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
			typename TMatrix::value_type block = conn.value();
			const size_t j = conn.index();

			// choose randomly the first dof to add whole row
			A(vConstrainingIndex[0][i], j) += block;
			A(constrainedIndex[i], j) = 0.0;
		}
	}
}

template <typename TVector>
void SplitAddRhs_Symmetric(TVector& rhs,
                         std::vector<size_t> & constrainedIndex,
                         std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

//	loop constrained indices
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	get constrained rhs
		typename TVector::value_type& val = rhs[constrainedIndex[i]];
		val *= 1./(vConstrainingIndex.size());

	// 	split equally on all constraining indices
		for(size_t j=0; j < vConstrainingIndex.size(); ++j)
			rhs[vConstrainingIndex[j][i]] += val;

	//	set rhs to zero for contrained index
		val = 0.0;
	}
}

template <typename TVector>
void SplitAddRhs_OneSide(TVector& rhs,
                       std::vector<size_t> & constrainedIndex,
                       std::vector<std::vector<size_t> >& vConstrainingIndex)
{
	//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndex.size(); ++i)
		UG_ASSERT(vConstrainingIndex[i].size() == constrainedIndex.size(),
				  "Wrong number of indices.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		typename TVector::value_type& val = rhs[constrainedIndex[i]];

		// choose randomly the first dof to add whole rhs (must be the same as for row)
		rhs[vConstrainingIndex[0][i]] += val;
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
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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

	//	adapt rhs
		SplitAddRhs_Symmetric(d, constrainedInd, vConstrainingInd);
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u,
           ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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

	//	adapt rhs
		SplitAddRhs_Symmetric(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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
		SplitAddRow_Symmetric(J, constrainedInd, vConstrainingInd);

	//	set interpolation
		SetInterpolation(J, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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

	//	set interpolation
		SetInterpolation(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		SplitAddRhs_Symmetric(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd,
                number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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
		InterpolateValues(u, constrainedInd, vConstrainingInd);
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	OneSide P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	Sym P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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

	//	adapt rhs
		SplitAddRhs_OneSide(d, constrainedInd, vConstrainingInd);
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u,
           ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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

	//	adapt rhs
		SplitAddRhs_OneSide(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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
		SplitAddRow_OneSide(J, constrainedInd, vConstrainingInd);

	//	set interpolation
		SetInterpolation(J, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<TDD> dd, number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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
		SplitAddRhs_OneSide(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd,
                number time)
{
//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	typename TDD::template traits<HangingVertex>::const_iterator iter, iterEnd;
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
		InterpolateValues(u, constrainedInd, vConstrainingInd);
	}
}


}; // namespace ug



#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__ */
