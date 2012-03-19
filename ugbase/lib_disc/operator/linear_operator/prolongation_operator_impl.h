/*
 * prolongation_operator_impl.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR_IMPL__

#include "prolongation_operator.h"

namespace ug{

template <typename TDD, typename TAlgebra>
void AssembleVertexProlongation(typename TAlgebra::matrix_type& mat,
                                const TDD& coarseDD, const TDD& fineDD,
								std::vector<bool>& vIsRestricted)
{
// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.local_finite_element_id(fct) != LFEID(LFEID::LAGRANGE, 1))
			UG_THROW_FATAL("AssembleVertexProlongation:"
				"Interpolation only implemented for Lagrange P1 functions.");

//  get subsethandler and grid
	const MultiGrid& grid = coarseDD.multi_grid();

//  get number of dofs on different levels
	const size_t numFineDoFs = fineDD.num_indices();
	const size_t numCoarseDoFs = coarseDD.num_indices();

//	check if grid distribution has dofs, otherwise skip creation since father
//	elements may not exist in parallel.
	if(numFineDoFs == 0 || numCoarseDoFs == 0) return;

//  resize matrix
	if(!mat.resize(numFineDoFs, numCoarseDoFs))
		UG_THROW_FATAL("AssembleVertexProlongation:"
				"Cannot resize Interpolation Matrix.");

//	clear restricted vector
	vIsRestricted.clear(); vIsRestricted.resize(numCoarseDoFs, false);

	std::vector<MultiIndex<2> > coarseMultInd, fineMultInd;

//  iterators
	typedef typename TDD::template traits<VertexBase>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<VertexBase>(si);
		iterEnd = fineDD.template end<VertexBase>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//  get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);
			EdgeBase* edge = dynamic_cast<EdgeBase*>(geomObj);
			Quadrilateral* quad = dynamic_cast<Quadrilateral*>(geomObj);
			ConstrainingFace* cFace = dynamic_cast<ConstrainingFace*>(geomObj);
			Face* face = dynamic_cast<Face*>(geomObj);
			Hexahedron* hexaeder = dynamic_cast<Hexahedron*>(geomObj);

			for(size_t fct = 0; fct < fineDD.num_fct(); fct++)
			{
				if(!fineDD.is_def_in_subset(fct, si)) continue;

			//  get global indices
				fineDD.inner_multi_indices(*iter, fct, fineMultInd);

			//  Check if father is Vertex
				if(vert != NULL)
				{
				//  get global indices
					coarseDD.inner_multi_indices(vert, fct, coarseMultInd);

					vIsRestricted[coarseMultInd[0][0]] = true;

					BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
								fineMultInd[0][1], coarseMultInd[0][1]) = 1.0;
					continue;
				}

			//  Check if father is Edge
				if(edge != NULL)
				{
					for(int i = 0; i < 2; ++i)
					{
						VertexBase* v = edge->vertex(i);

					//  get global indices
						coarseDD.inner_multi_indices(v, fct, coarseMultInd);

						vIsRestricted[coarseMultInd[0][0]] = true;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.5;
					}
					continue;
				}

			//  Check if father is Quad
				if(quad != NULL || cFace != NULL)
				{
					for(int i = 0; i < 4; ++i)
					{
						VertexBase* v = face->vertex(i);

					//  get global indices
						coarseDD.inner_multi_indices(v, fct, coarseMultInd);

						vIsRestricted[coarseMultInd[0][0]] = true;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.25;
					}
					continue;
				}

			//  Check if father is Hexaeder
				if(hexaeder != NULL)
				{
					for(int i = 0; i < 8; ++i)
					{
						VertexBase* v = hexaeder->vertex(i);

					//  get global indices
						coarseDD.inner_multi_indices(v, fct, coarseMultInd);

						vIsRestricted[coarseMultInd[0][0]] = true;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.125;
					}
					continue;
				}

				UG_THROW_FATAL("AssembleVertexProlongation: Element Father not detected.");
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 	P1Prolongation
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void P1Prolongation<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	m_spApproxSpace = approxSpace;
}

template <typename TDomain, typename TAlgebra>
void P1Prolongation<TDomain, TAlgebra>::set_levels(GridLevel coarseLevel, GridLevel fineLevel)
{
	m_fineLevel = fineLevel;
	m_coarseLevel = coarseLevel;
	if(m_fineLevel.level() - m_coarseLevel.level() != 1)
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::set_levels:"
				" Can only project between successive level.");
	if(m_fineLevel.type() != GridLevel::LEVEL ||
	   m_coarseLevel.type() != GridLevel::LEVEL)
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::set_levels:"
				" Can only project between level dof distributions, but fine="
				<<m_fineLevel<<", coarse="<<m_coarseLevel);
}

template <typename TDomain, typename TAlgebra>
void P1Prolongation<TDomain, TAlgebra>::init()
{
	if(!m_spApproxSpace.valid())
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::init: "
				"Approximation Space not set. Cannot init Projection.");

	m_matrix.resize(0,0);

	try{
	if(m_coarseLevel.type() == GridLevel::LEVEL)
		AssembleVertexProlongation<LevelDoFDistribution, algebra_type>
		(m_matrix,
		 *m_spApproxSpace->level_dof_distribution(m_coarseLevel.level()),
		 *m_spApproxSpace->level_dof_distribution(m_fineLevel.level()),
		 m_vIsRestricted);
	} UG_CATCH_THROW("P1Prolongation<TDomain, TAlgebra>::init:"
				"Cannot assemble interpolation matrix.");

	#ifdef UG_PARALLEL
		m_matrix.set_storage_type(PST_CONSISTENT);
	#endif

	m_bInit = true;
}

template <typename TDomain, typename TAlgebra>
void P1Prolongation<TDomain, TAlgebra>::apply(vector_type& uFineOut, const vector_type& uCoarseIn)
{
//	Check, that operator is initiallized
	if(!m_bInit)
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::apply:"
				" Operator not initialized.");

//	Some Assertions
	UG_ASSERT(uFineOut.size() == m_matrix.num_rows(),
				  "Vector and Row sizes have to match!");
	UG_ASSERT(uCoarseIn.size() == m_matrix.num_cols(),
				  "Vector and Column sizes have to match!");

//	Apply Matrix
	if(!m_matrix.apply(uFineOut, uCoarseIn))
	{
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::apply: " <<
					"Cannot apply matrix. "
#ifdef UG_PARALLEL
					<< "(Type uCoarse = " <<uCoarseIn.get_storage_mask()
#endif
		);
	}

//	Set dirichlet nodes to zero again
//	todo: We could handle this by eliminating dirichlet rows as well
	try{
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		m_vConstraint[i]->adjust_defect(uFineOut, uFineOut, m_fineLevel);

	}UG_CATCH_THROW("P1Prolongation<TDomain, TAlgebra>::apply: "
					"Error while setting dirichlet defect to zero.");
}

template <typename TDomain, typename TAlgebra>
void P1Prolongation<TDomain, TAlgebra>::apply_transposed(vector_type& uCoarseOut, const vector_type& uFineIn)
{
//	Check, that operator is initiallized
	if(!m_bInit)
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::apply_transposed:"
				"Operator not initialized.");

	vector_type	uTmp; uTmp.resize(uCoarseOut.size());

//	Some Assertions
	UG_ASSERT(uFineIn.size() == m_matrix.num_rows(),
				  "Vector and Row sizes have to match!");
	UG_ASSERT(uCoarseOut.size() == m_matrix.num_cols(),
				  "Vector and Column sizes have to match!");

//	Apply transposed matrix
	if(!m_matrix.apply_transposed(uTmp, uFineIn))
		UG_THROW_FATAL("P1Prolongation<TDomain, TAlgebra>::apply_transposed:"
				" Cannot apply transposed matrix.");

	uTmp *= m_dampRes;

//	Copy only restricted values
//	This is needed in adaptive situations, where the defect that has been
//	projected from the surface should remain and only hidden (i.e.
//	indices with children) should be changed.
	for(size_t i = 0; i < uTmp.size(); ++i)
		if(m_vIsRestricted[i])
			uCoarseOut[i] = uTmp[i];

//	Set dirichlet nodes to zero again
//	todo: We could handle this by eliminating dirichlet columns as well
	try{
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		m_vConstraint[i]->adjust_defect(uCoarseOut, uCoarseOut, m_coarseLevel);
	} UG_CATCH_THROW("ProjectionOperator::apply_transposed: "
					"Error while setting dirichlet defect to zero.");
}

template <typename TDomain, typename TAlgebra>
SmartPtr<IProlongationOperator<TAlgebra> >
P1Prolongation<TDomain, TAlgebra>::clone()
{
	SmartPtr<P1Prolongation> op(new P1Prolongation);
	op->set_approximation_space(m_spApproxSpace);
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		op->add_constraint(m_vConstraint[i]);
	op->set_restriction_damping(m_dampRes);
	return op;
}

template <typename TDomain, typename TAlgebra>
void
P1Prolongation<TDomain, TAlgebra>::add_constraint(SmartPtr<IConstraint<algebra_type> > pp)
{
//	add only once
	if(std::find(m_vConstraint.begin(), m_vConstraint.end(), pp) !=
			m_vConstraint.end()) return;
	m_vConstraint.push_back(pp);
}

template <typename TDomain, typename TAlgebra>
void
P1Prolongation<TDomain, TAlgebra>::remove_constraint(SmartPtr<IConstraint<algebra_type> > pp)
{
	m_vConstraint.erase(m_vConstraint.begin(),
	                     std::remove(m_vConstraint.begin(), m_vConstraint.end(), pp));
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR_IMPL__ */
