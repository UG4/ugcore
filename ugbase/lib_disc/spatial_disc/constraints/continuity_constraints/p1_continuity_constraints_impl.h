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

/// returns the vertices of the object constraining a hanging vertex
void CollectConstraining(std::vector<VertexBase*>& vConstrainingVrt,
                         HangingVertex* hgVrt,
                         bool bClearContainer = true)
{
//	clear container
	if(bClearContainer) vConstrainingVrt.clear();

//	switch constraining parent
	switch(hgVrt->get_parent_base_object_type_id())
	{
	case EDGE:
	{
	//	cast to constraining edge
		ConstrainingEdge* constrainingEdge =
				dynamic_cast<ConstrainingEdge*>(hgVrt->get_parent());

	//	check that edge is correct
		if(constrainingEdge == NULL)
			UG_THROW_FATAL("Parent element should be "
						"constraining edge, but is not.");

	//	get constraining vertices
		for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
		{
		//	get constrained edge
			ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
												constrainingEdge->constrained_edge(i_cde));

		//	check
			if(constrainedEdge == NULL)
				UG_THROW_FATAL("Child element should be "
							"constrained edge, but is not.");

		//	get non-hanging vertex
			VertexBase* vrt = GetConnectedVertex(constrainedEdge, hgVrt);

		//	push back in list of interpolation vertices
			vConstrainingVrt.push_back(vrt);
		}
	}
		break;
	case FACE:
	{
	//	cast to constraining quadrilateral
		ConstrainingQuadrilateral* bigQuad =
				dynamic_cast<ConstrainingQuadrilateral*>(hgVrt->get_parent());

	//	check that quad is correct
		if(bigQuad == NULL)
			UG_THROW_FATAL("Parent element should be "
							"constraining quad, but is not.");

	//	get constraining vertices
	//	\todo: This is only valid for a surface grid!!!
		for(size_t i_cf=0; i_cf < bigQuad->num_constrained_faces(); ++i_cf)
		{
			Face* face = bigQuad->constrained_face(i_cf);

			VertexBase* vrt = NULL;
			size_t i_vrt = 0;
			for(i_vrt = 0; i_vrt < face->num_vertices(); ++i_vrt)
			{
				vrt = face->vertex(i_vrt);
				if(hgVrt != vrt && dynamic_cast<HangingVertex*>(vrt) == NULL)
					break;
			}
			if(i_vrt == face->num_vertices())
				UG_THROW_FATAL("ERROR: Vertex not detected.\n");

			vConstrainingVrt.push_back(vrt);
		}
	}
		break;
	default: UG_THROW_FATAL("Parent element of hang. vertex wrong.");
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
		SplitAddRow(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		HandleRhs(rhs, constrainedInd, vConstrainingInd);
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
		InterpolateValues(u, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
SplitAddRow(matrix_type& A	,
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
		typename matrix_type::value_type& block
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
		for(typename matrix_type::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
		//	skip self-coupling (already handled)
			const size_t j = conn.index();
			if(j == constrainedIndex[i]) continue;

		//	get coupling entry
			typename matrix_type::value_type block = conn.value();

		//	get transposed coupling entry
			typename matrix_type::value_type blockT = A(j, constrainedIndex[i]);

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

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
SetInterpolation(matrix_type& A,
                 std::vector<size_t> & constrainedIndex,
                 std::vector<std::vector<size_t> >& vConstrainingIndices)
{
//	check number of indices passed
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.\n");

//	loop all constrained dofs
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
	//	remove all couplings
		for(typename matrix_type::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
			conn.value() = 0.0;
		}

	//	set diag of row to identity
		A(constrainedIndex[i], constrainedIndex[i]) = 1.0;

	//	set coupling to all contraining dofs the inverse of the
	//	number of contraining dofs
		number frac = -1.0/(vConstrainingIndices.size());
		for(size_t j=0; j < vConstrainingIndices.size();++j)
			A(constrainedIndex[i], vConstrainingIndices[j][i]) = frac;
	}
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
HandleRhs(vector_type& rhs,
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
		typename vector_type::value_type& val = rhs[constrainedIndex[i]];
		val *= 1./(vConstrainingIndices.size());

	// 	split equally on all constraining indices
		for(size_t j=0; j < vConstrainingIndices.size(); ++j)
			rhs[vConstrainingIndices[j][i]] += val;

	//	set rhs to zero for contrained index
		val = 0.0;
	}
}

template <typename TDomain, typename TAlgebra>
void
SymP1Constraints<TDomain,TAlgebra>::
InterpolateValues(vector_type& u,
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
		typename vector_type::value_type& val = u[constrainedIndex[i]];
		const number scale = 1./(vConstrainingIndices.size());

		val = 0.0;

	// 	split equally on all constraining indices
		for(size_t j=0; j < vConstrainingIndices.size(); ++j)
		{
			typename vector_type::value_type entry = u[vConstrainingIndices[j][i]];
			entry *= scale;
			val += entry;
		}
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
		SplitAddRow(mat, constrainedInd, vConstrainingInd);

	//	Set interpolation
		SetInterpolation(mat, constrainedInd, vConstrainingInd);

	//	adapt rhs
		HandleRhs(rhs, constrainedInd, vConstrainingInd);
	}
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
SplitAddRow(matrix_type& A,
            std::vector<size_t> & constrainedIndex,
            std::vector<std::vector<size_t> >& vConstrainingIndices)
{
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wring number of indices. Cannot split row.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		for(typename matrix_type::row_iterator conn = A.begin_row(constrainedIndex[i]);
				conn != A.end_row(constrainedIndex[i]); ++conn)
		{
			typename matrix_type::value_type block = conn.value();
			const size_t j = conn.index();

			// choose randomly the first dof to add whole row
			A(vConstrainingIndices[0][i], j) += block;
			A(constrainedIndex[i], j) = 0.0;
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
SetInterpolation(matrix_type& A,
                 std::vector<size_t> & constrainedIndex,
                 std::vector<std::vector<size_t> >& vConstrainingIndices)
{
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.");

	const number scale = -1./(vConstrainingIndices.size());
	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
			A(constrainedIndex[i], constrainedIndex[i]) = 1.0;
			for(size_t j = 0; j < vConstrainingIndices.size(); ++j)
			{
				A(constrainedIndex[i], vConstrainingIndices[j][i]) = scale;
			}
	}
}

template <typename TDomain, typename TAlgebra>
void
OneSideP1Constraints<TDomain,TAlgebra>::
HandleRhs(vector_type& rhs,
          std::vector<size_t> & constrainedIndex,
          std::vector<std::vector<size_t> >& vConstrainingIndices)
{
	for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
		if(vConstrainingIndices[i].size() != constrainedIndex.size())
			UG_THROW_FATAL("Wrong number of indices. Cannot split row.");

	for(size_t i = 0; i < constrainedIndex.size(); ++i)
	{
		typename vector_type::value_type& val = rhs[constrainedIndex[i]];

		// choose randomly the first dof to add whole rhs (must be the same as for row)
		rhs[vConstrainingIndices[0][i]] += val;
		val = 0.0;
	}
}

}; // namespace ug



#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS_IMPL__ */
