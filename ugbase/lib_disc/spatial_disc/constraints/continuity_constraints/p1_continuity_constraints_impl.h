/*
 * p1_continuity_constraints_impl.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
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
		const row_iterator iterEnd = A.end_row(constrainedIndex[i]);
		for(row_iterator conn = A.begin_row(constrainedIndex[i]); conn != iterEnd; ++conn)
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
	//	we can work directly on the entry (modifying it) since row of
	//	constraints will be set to interpolation afterwards
		block_type& block = A(constrainedIndex[i], constrainedIndex[i]);

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
			const block_type& block = conn.value();
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

	//	set rhs to zero for contrained index
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	Sym P1 Constraints
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

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
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

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
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

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
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constrained vertices
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

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
void
SymP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
                number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for SymP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constraining edges
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

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

template<int dim>
struct SortVertexPos {

		SortVertexPos(SmartPtr<Domain<dim, MultiGrid, MGSubsetHandler> > spDomain)
			: m_aaPos(spDomain->position_accessor())
		{}

		inline bool operator() (VertexBase* vrt1, VertexBase* vrt2)
			{UG_THROW(dim <<" not implemented.");}

	protected:
  	  typename Domain<dim, MultiGrid, MGSubsetHandler>::position_accessor_type& m_aaPos;
};

template<>
inline bool SortVertexPos<1>::operator() (VertexBase* vrt1, VertexBase* vrt2)
{
	if(m_aaPos[vrt1][0] < m_aaPos[vrt2][0]) {
		return true;
	}
	return false;
}

template<>
inline bool SortVertexPos<2>::operator() (VertexBase* vrt1, VertexBase* vrt2)
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
inline bool SortVertexPos<3>::operator() (VertexBase* vrt1, VertexBase* vrt2)
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



template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

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

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

#ifdef UG_PARALLEL
		std::sort(vConstrainingVrt.begin(), vConstrainingVrt.end(), sortVertexPos);
#endif

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
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

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

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

#ifdef UG_PARALLEL
		std::sort(vConstrainingVrt.begin(), vConstrainingVrt.end(), sortVertexPos);
#endif

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
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

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

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

#ifdef UG_PARALLEL
		std::sort(vConstrainingVrt.begin(), vConstrainingVrt.end(), sortVertexPos);
#endif

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
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

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

	//	get constraining vertices
		CollectConstraining(vConstrainingVrt, hgVrt);

	//	resize constraining indices
		vConstrainingInd.clear();
		vConstrainingInd.resize(vConstrainingVrt.size());

#ifdef UG_PARALLEL
		std::sort(vConstrainingVrt.begin(), vConstrainingVrt.end(), sortVertexPos);
#endif

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
void
OneSideP1Constraints<TDomain,TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
                number time)
{
	if(this->m_AssIndex.index_set)
		UG_THROW("index-wise assemble routine is not "
				"yet implemented for OneSideP1Constraints \n");

//	storage for indices and vertices
	std::vector<std::vector<size_t> > vConstrainingInd;
	std::vector<size_t>  constrainedInd;
	std::vector<VertexBase*> vConstrainingVrt;

//	get begin end of hanging vertices
	DoFDistribution::traits<ConstrainedVertex>::const_iterator iter, iterEnd;
	iter = dd->begin<ConstrainedVertex>();
	iterEnd = dd->end<ConstrainedVertex>();

//	loop constraining edges
	for(; iter != iterEnd; ++iter)
	{
	//	get hanging vert
		ConstrainedVertex* hgVrt = *iter;

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
