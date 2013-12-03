/*
 * projection_operator_impl.h
 *
 *  Created on: 01.08.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR_IMPL__

#include "projection_operator.h"

namespace ug{

/**
 * This functions assembles the interpolation matrix between to
 * grid levels using only the Vertex degrees of freedom.
 *
 * \param[out]	mat 			Assembled interpolation matrix that interpolates u -> v
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	coarseLevel		Coarse Level index
 * \param[in]	fineLevel		Fine Level index
 */
template <typename TAlgebra>
void AssembleInjectionForP1Lagrange(typename TAlgebra::matrix_type& mat,
                                    const DoFDistribution& coarseDD,
                                    const DoFDistribution& fineDD)
{
	PROFILE_FUNC_GROUP("gmg");
//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.local_finite_element_id(fct).type() != LFEID::LAGRANGE ||
			fineDD.local_finite_element_id(fct).order() != 1)
			UG_THROW("AssembleInjectionForP1Lagrange: "
					"Interpolation only implemented for Lagrange P1 functions.");

// 	get MultiGrid
	const MultiGrid& grid = *coarseDD.multi_grid();

// 	get number of dofs on different levels
	const size_t numFineDoFs = fineDD.num_indices();
	const size_t numCoarseDoFs = coarseDD.num_indices();

// 	resize matrix
	mat.resize_and_clear(numCoarseDoFs, numFineDoFs);

	std::vector<size_t> coarseInd, fineInd;

// 	Vertex iterators
	typedef DoFDistribution::traits<Vertex>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

	iterBegin = fineDD.template begin<Vertex>();
	iterEnd = fineDD.template end<Vertex>();

// 	loop nodes of fine subset
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	// 	get father
		GeometricObject* geomObj = grid.get_parent(*iter);
		VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);

	//	Check if father is Vertex
		if(vert != NULL)
		{
			// get global indices
			coarseDD.inner_algebra_indices(vert, coarseInd);
		}
		else continue;

	// 	get global indices
		fineDD.inner_algebra_indices(*iter, fineInd);

		for(size_t i = 0; i < coarseInd.size(); ++i)
			mat(coarseInd[i], fineInd[i]) = 1.0;
	}
}

template <int dim, typename TAlgebra>
void AssembleInjectionByAverageOfChildren(typename TAlgebra::matrix_type& mat,
                                          const DoFDistribution& coarseDD,
                                          const DoFDistribution& fineDD)
{
	PROFILE_FUNC_GROUP("gmg");
// 	get MultiGrid
	const MultiGrid& grid = *coarseDD.multi_grid();

	std::vector<size_t> coarseInd, fineInd;

// 	Vertex iterators
	typedef typename DoFDistribution::dim_traits<dim>::const_iterator const_iterator;
	typedef typename DoFDistribution::dim_traits<dim>::geometric_base_object Element;
	const_iterator iter, iterBegin, iterEnd;

	iterBegin = coarseDD.template begin<Element>();
	iterEnd = coarseDD.template end<Element>();

// 	loop elements of coarse subset
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
	//	get element
		Element* coarseElem = *iter;

	//  get children
		const size_t numChild = grid.num_children<Element, Element>(coarseElem);
		if(numChild == 0) continue;

	// get global indices
		coarseDD.inner_algebra_indices(coarseElem, coarseInd);

		for(size_t c = 0; c < numChild; ++c)
		{
			Element* child = grid.get_child<Element, Element>(coarseElem,  c);

			fineDD.inner_algebra_indices(child, fineInd);

			for(size_t i = 0; i < coarseInd.size(); ++i)
				mat(coarseInd[i], fineInd[i]) = 1.0/(numChild);
		}
	}
}

template <typename TAlgebra>
void AssembleInjectionByAverageOfChildren(typename TAlgebra::matrix_type& mat,
                                          const DoFDistribution& coarseDD,
                                          const DoFDistribution& fineDD)
{
// 	get number of dofs on different levels
	const size_t numFineDoFs = fineDD.num_indices();
	const size_t numCoarseDoFs = coarseDD.num_indices();

// 	resize matrix
	mat.resize_and_clear(numCoarseDoFs, numFineDoFs);

	if(coarseDD.max_dofs(VERTEX)) AssembleInjectionByAverageOfChildren<0, TAlgebra>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(EDGE)) AssembleInjectionByAverageOfChildren<1, TAlgebra>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(FACE)) AssembleInjectionByAverageOfChildren<2, TAlgebra>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(VOLUME)) AssembleInjectionByAverageOfChildren<3, TAlgebra>(mat, coarseDD, fineDD);
}

template <typename TDomain, typename TAlgebra>
template <typename TElem>
void InjectionTransfer<TDomain, TAlgebra>::
set_identity_on_pure_surface(matrix_type& mat,
                             const DoFDistribution& coarseDD, const DoFDistribution& fineDD)
{
	PROFILE_FUNC_GROUP("gmg");

	std::vector<size_t> vCoarseIndex, vFineIndex;
	const MultiGrid& mg = *coarseDD.multi_grid();

//  iterators
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	for(int si = 0; si < coarseDD.num_subsets(); ++si)
	{
		iterBegin = coarseDD.template begin<TElem>(si);
		iterEnd = coarseDD.template end<TElem>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			TElem* coarseElem = *iter;

			const size_t numChild = mg.num_children<TElem, TElem>(coarseElem);
			if(numChild != 0) continue;

		//	get indices
			coarseDD.inner_algebra_indices(coarseElem, vCoarseIndex);
			fineDD.inner_algebra_indices(coarseElem, vFineIndex);
			UG_ASSERT(vCoarseIndex.size() == vFineIndex.size(), "Size mismatch");

		//	set identity
			for(size_t i = 0; i < vCoarseIndex.size(); ++i)
				mat(vCoarseIndex[i], vFineIndex[i]) = 1.0;
		}
	}
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
set_identity_on_pure_surface(matrix_type& mat,
                             const DoFDistribution& coarseDD, const DoFDistribution& fineDD)
{
	if(coarseDD.max_dofs(VERTEX)) set_identity_on_pure_surface<VertexBase>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(EDGE)) set_identity_on_pure_surface<EdgeBase>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(FACE)) set_identity_on_pure_surface<Face>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(VOLUME)) set_identity_on_pure_surface<Volume>(mat, coarseDD, fineDD);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 	InjectionTransfer
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	m_spApproxSpace = approxSpace;
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
set_levels(GridLevel coarseLevel, GridLevel fineLevel)
{
	m_fineLevel = fineLevel;
	m_coarseLevel = coarseLevel;

	if(m_fineLevel.level() - m_coarseLevel.level() != 1)
		UG_THROW("InjectionTransfer::set_levels:"
				" Can only project between successive level.");
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::init()
{
	PROFILE_FUNC_GROUP("gmg");
	if(!m_spApproxSpace.valid())
		UG_THROW("InjectionTransfer::init: "
				"Approximation Space not set. Cannot init Projection.");

// 	check only lagrange P1 functions
	bool P1LagrangeOnly = true;
	for(size_t fct = 0; fct < m_spApproxSpace->num_fct(); ++fct)
		if(m_spApproxSpace->local_finite_element_id(fct).type() != LFEID::LAGRANGE ||
				m_spApproxSpace->local_finite_element_id(fct).order() != 1)
			P1LagrangeOnly = false;

	try{
		if(P1LagrangeOnly)
		{
			AssembleInjectionForP1Lagrange<TAlgebra>
			(m_matrix,
			 *m_spApproxSpace->dof_distribution(m_coarseLevel),
			 *m_spApproxSpace->dof_distribution(m_fineLevel));
		}
		else
		{
			AssembleInjectionByAverageOfChildren<TAlgebra>
			(m_matrix,
			 *m_spApproxSpace->dof_distribution(m_coarseLevel),
			 *m_spApproxSpace->dof_distribution(m_fineLevel));
		}

	} UG_CATCH_THROW("InjectionTransfer::init():"
						" Cannot assemble interpolation matrix.");

	if(m_coarseLevel.is_surface()){
		set_identity_on_pure_surface(m_matrix, *m_spApproxSpace->dof_distribution(m_coarseLevel), *m_spApproxSpace->dof_distribution(m_fineLevel));
	}

	#ifdef UG_PARALLEL
		m_matrix.set_storage_type(PST_CONSISTENT);
	#endif

	m_bInit = true;
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
prolongate(vector_type& uFine, const vector_type& uCoarse)
{
	PROFILE_FUNC_GROUP("gmg");
//	Check, that operator is initiallized
	if(!m_bInit)
		UG_THROW("InjectionTransfer::apply:"
				" Operator not initialized.");

//	Some Assertions
	UG_ASSERT(uFine.size() >= m_matrix.num_rows(),
			  "Vector [size= " << uFine.size() << "] and Rows [size= "
			  << m_matrix.num_rows() <<"] sizes have to match!");
	UG_ASSERT(uCoarse.size() >= m_matrix.num_cols(),	"Vector [size= "
			  << uCoarse.size() << "] and Cols [size= " <<
			  m_matrix.num_cols() <<"] sizes have to match!");

//	Apply matrix
	if(!m_matrix.apply_transposed(uFine, uCoarse))
		UG_THROW("InjectionTransfer::apply: Cannot apply matrix.");
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
do_restrict(vector_type& uCoarse, const vector_type& uFine)
{
	PROFILE_FUNC_GROUP("gmg");
//	Check, that operator is initialized
	if(!m_bInit)
		UG_THROW("InjectionTransfer::apply_transposed:"
				"Operator not initialized.");

//	Some Assertions
	UG_ASSERT(uFine.size() >= m_matrix.num_cols(),
			  "Vector [size= " << uFine.size() << "] and Cols [size= "
			  << m_matrix.num_cols() <<"] sizes have to match!");
	UG_ASSERT(uCoarse.size() >= m_matrix.num_rows(),	"Vector [size= "
			  << uCoarse.size() << "] and Rows [size= " <<
			  m_matrix.num_rows() <<"] sizes have to match!");

//	Apply matrix
	try{
		m_matrix.apply_ignore_zero_rows(uCoarse, 1.0, uFine);
	}
	UG_CATCH_THROW("InjectionTransfer::apply_transposed:"
						" Cannot apply transposed matrix.");
}


template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferOperator<TDomain, TAlgebra> >
InjectionTransfer<TDomain, TAlgebra>::clone()
{
	SmartPtr<InjectionTransfer> op(new InjectionTransfer);
	op->set_approximation_space(m_spApproxSpace);
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		op->add_constraint(m_vConstraint[i]);
	return op;
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
add_constraint(SmartPtr<IConstraint<TAlgebra> > pp)
{
	UG_THROW("Not Implemented.");
//	add only once
	if(std::find(m_vConstraint.begin(), m_vConstraint.end(), pp) !=
			m_vConstraint.end()) return;
	m_vConstraint.push_back(pp);
}

template <typename TDomain, typename TAlgebra>
void InjectionTransfer<TDomain, TAlgebra>::
remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp)
{
	UG_THROW("Not Implemented.");
	m_vConstraint.erase(m_vConstraint.begin(),
	                     std::remove(m_vConstraint.begin(), m_vConstraint.end(), pp));
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR_IMPL__ */
