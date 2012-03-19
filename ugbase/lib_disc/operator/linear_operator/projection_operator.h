/*
 * projection_operator.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "transfer_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

// TODO: This function should be put to an util file
/**
 * This functions assembles the interpolation matrix between to
 * grid levels using only the Vertex degrees of freedom.
 *
 * \param[out]	mat 			Assembled interpolation matrix that interpolates u -> v
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	coarseLevel		Coarse Level index
 * \param[in]	fineLevel		Fine Level index
 */
template <typename TDD, typename TAlgebra>
void AssembleVertexProjection(typename TAlgebra::matrix_type& mat,
                              const TDD& coarseDD, const TDD& fineDD)
{
//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.local_finite_element_id(fct) != LFEID(LFEID::LAGRANGE, 1))
			UG_THROW_FATAL("AssembleVertexProjection: "
					"Interpolation only implemented for Lagrange P1 functions.");

// 	get MultiGrid
	const MultiGrid& grid = coarseDD.multi_grid();

// 	get number of dofs on different levels
	const size_t numFineDoFs = fineDD.num_indices();
	const size_t numCoarseDoFs = coarseDD.num_indices();

// 	resize matrix
	if(!mat.resize(numCoarseDoFs, numFineDoFs))
		UG_THROW_FATAL("AssembleVertexProjection: "
				"Cannot resize Interpolation Matrix.");

	std::vector<size_t> coarseInd, fineInd;

// 	Vertex iterators
	typedef typename TDD::template traits<VertexBase>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<Vertex>(si);
		iterEnd = fineDD.template end<Vertex>(si);

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
}



template <typename TDomain, typename TAlgebra>
class P1Projection :
	virtual public IProjectionOperator<	typename TAlgebra::vector_type,
										typename TAlgebra::vector_type>
{
	public:
	// 	Type of algebra
		typedef TAlgebra algebra_type;

	//	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	//	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	//	Type of Domain
		typedef TDomain domain_type;

	public:
	//	Constructor
		P1Projection() : m_bInit(false) {}

	//	Constructor
		P1Projection(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
			m_spApproxSpace(approxSpace), m_bInit(false)
		{}

	//	Set Approximation Space
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
		{
			m_spApproxSpace = approxSpace;
		}

	//	Set approximation level
		void set_levels(GridLevel coarseLevel, GridLevel fineLevel)
		{
			m_fineLevel = fineLevel;
			m_coarseLevel = coarseLevel;
			if(m_fineLevel.level() - m_coarseLevel.level() != 1)
				UG_THROW_FATAL("P1Projection::set_levels:"
						" Can only project between successive level.");
			if(m_fineLevel.type() != GridLevel::LEVEL ||
			   m_coarseLevel.type() != GridLevel::LEVEL)
				UG_THROW_FATAL("P1Projection::set_levels:"
						" Can only project between level dof distributions.");
		}

	public:
	//	Init operator
		virtual void init()
		{
			if(!m_spApproxSpace.valid())
				UG_THROW_FATAL("P1Projection::init: "
						"Approximation Space not set. Cannot init Projection.");

			m_matrix.resize(0,0);

			try{
			if(m_coarseLevel.type() == GridLevel::LEVEL)
				AssembleVertexProjection<LevelDoFDistribution, algebra_type>
				(m_matrix,
				 *m_spApproxSpace->level_dof_distribution(m_coarseLevel.level()),
				 *m_spApproxSpace->level_dof_distribution(m_fineLevel.level()));
			} UG_CATCH_THROW("TransferOperator::init():"
						" Cannot assemble interpolation matrix.");

			#ifdef UG_PARALLEL
				m_matrix.set_storage_type(PST_CONSISTENT);
			#endif

			m_bInit = true;
		}

		virtual void init(const vector_type& u){init();}

	// 	Project uFine to uCoarse; uCoarse = P(uFine);
		virtual void apply(vector_type& uCoarseOut, const vector_type& uFineIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
				UG_THROW_FATAL("P1Projection::apply:"
						" Operator not initialized.");

		//	Some Assertions
			UG_ASSERT(uCoarseOut.size() == m_matrix.num_rows(),
			          "Vector [size= " << uCoarseOut.size() << "] and Row [size= "
			          << m_matrix.num_rows() <<"] sizes have to match!");
			UG_ASSERT(uFineIn.size() == m_matrix.num_cols(),	"Vector [size= "
			          << uFineIn.size() << "] and Column [size= " <<
			          m_matrix.num_cols() <<"] sizes have to match!");

		//	Apply matrix
			if(!m_matrix.apply(uCoarseOut, uFineIn))
				UG_THROW_FATAL("P1Projection::apply: Cannot apply matrix.");
		}

	// 	Apply Transposed Operator u = L^T*f
		virtual void apply_transposed(vector_type& uFineOut, const vector_type& uCoarseIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
				UG_THROW_FATAL("P1Projection::apply_transposed:"
						"Operator not initialized.");

		//	Some Assertions
			UG_ASSERT(uCoarseIn.size() == m_matrix.num_rows(),
					  "Vector [size= " << uCoarseIn.size() << "] and Cols [size= "
					  << m_matrix.num_rows() <<"] sizes have to match!");
			UG_ASSERT(uFineOut.size() == m_matrix.num_cols(),	"Vector [size= "
					  << uFineOut.size() << "] and Rows [size= " <<
					  m_matrix.num_cols() <<"] sizes have to match!");

		//	Apply matrix
			if(!m_matrix.apply_transposed(uFineOut, uCoarseIn))
				UG_THROW_FATAL("P1Projection::apply_transposed:"
						" Cannot apply transposed matrix.");
		}

	// 	Apply sub not implemented
		virtual void apply_sub(vector_type& u, const vector_type& v)
		{
			UG_THROW_FATAL("Not Implemented.");
		}

	//	Destructor
		~P1Projection()
		{
		}

		virtual SmartPtr<IProjectionOperator<vector_type> > clone()
		{
			SmartPtr<P1Projection> op(new P1Projection);
			op->set_approximation_space(m_spApproxSpace);
			return op;
		}

	protected:
		matrix_type m_matrix;

		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		GridLevel m_fineLevel;
		GridLevel m_coarseLevel;

		bool m_bInit;
};


}

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__ */
