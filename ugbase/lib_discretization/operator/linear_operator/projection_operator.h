/*
 * projection_operator.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_discretization/spacial_discretization/post_process/post_process_interface.h"

#ifdef UG_PARALLEL
#include "lib_discretization/parallelization/parallelization_util.h"
#endif

namespace ug{

// TODO: This function should be put to an util file
/** AssembleVertexProjection
 *
 * This functions assembles the interpolation matrix between to
 * grid levels using only the Vertex degrees of freedom.
 *
 * \param[in] 	uCoarse			Grid function on coarse level
 * \param[in] 	uFine 			Grid function on fine level
 * \param[out]	mat 			Assembled interpolation matrix that interpolates u -> v
 *
 */
template <typename TApproximationSpace, typename TMatrix>
bool AssembleVertexProjection(TMatrix& mat, TApproximationSpace& approxSpace, size_t coarseLevel, size_t fineLevel)
{
//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& coarseDoFDistr = approxSpace.get_level_dof_distribution(coarseLevel);
	const typename TApproximationSpace::dof_distribution_type& fineDoFDistr = approxSpace.get_level_dof_distribution(fineLevel);

//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDoFDistr.num_fct(); ++fct)
		if(fineDoFDistr.local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
			{UG_LOG("Interpolation only implemented for Lagrange P1 functions.\n"); return false;}

	// get MultiGrid
	MultiGrid& grid = approxSpace.get_domain().get_grid();

	// get number of dofs on different levels
	const size_t numFineDoFs = fineDoFDistr.num_dofs();
	const size_t numCoarseDoFs = coarseDoFDistr.num_dofs();

	// create matrix
	if(!mat.destroy())
		{UG_LOG("Cannot destroy Interpolation Matrix.\n"); return false;}
	if(!mat.create(numCoarseDoFs, numFineDoFs))
		{UG_LOG("Cannot create Interpolation Matrix.\n"); return false;}

	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type fineInd;
	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type coarseInd;

	// Vertex iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

	// loop subsets of fine level
	for(int si = 0; si < fineDoFDistr.num_subsets(); ++si)
	{
		iterBegin = fineDoFDistr.template begin<Vertex>(si);
		iterEnd = fineDoFDistr.template end<Vertex>(si);

		// loop nodes of fine subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);

			// Check if father is Vertex
			if(vert != NULL)
			{
				// get global indices
				coarseDoFDistr.get_inner_algebra_indices(vert, coarseInd);
			}
			else continue;

			// get global indices
			fineDoFDistr.get_inner_algebra_indices(*iter, fineInd);

			for(size_t i = 0; i < coarseInd.size(); ++i)
				mat(coarseInd[i], fineInd[i]) = 1.0;
		}
	}
	return true;
}



template <typename TApproximationSpace, typename TAlgebra>
class P1ProjectionOperator :
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

	//	Type of Approximation Space
		typedef TApproximationSpace approximation_space_type;

	//	Type of DoFDistribution
		typedef typename TApproximationSpace::dof_distribution_type dof_distribution_type;

	public:
	//	Constructor
		P1ProjectionOperator() :
			m_pApproximationSpace(NULL), m_fineLevel(0), m_coarseLevel(0), m_bInit(false)
		{}

	//	Set Approximation Space
		void set_approximation_space(approximation_space_type& approxSpace)
		{
			m_pApproximationSpace = &approxSpace;
		}

	//	Set approximation level
		bool set_levels(size_t coarseLevel, size_t fineLevel)
		{
			m_fineLevel = fineLevel;
			m_coarseLevel = coarseLevel;
			if(m_fineLevel - m_coarseLevel != 1)
			{
				UG_LOG("ERROR in ProjectionOperator::set_levels:"
						" Can only project between successiv level.\n");
				return false;
			}
			return true;
		}

	public:
	//	Init operator
		virtual bool init()
		{
			if(m_pApproximationSpace == NULL)
			{
				UG_LOG("ERROR in 'ProjectionOperator::init': Approximation Space not set. Cannot init Projection.\n");
				return false;
			}

			if(m_fineLevel - m_coarseLevel != 1)
			{
				UG_LOG("ERROR in ProjectionOperator::init:"
						" Can only project between successiv level.\n");
				return false;
			}

			if(!AssembleVertexProjection(m_matrix, *m_pApproximationSpace, m_coarseLevel, m_fineLevel))
			{
				UG_LOG("ERROR in 'TransferOperator::prepare(u,v)': Cannot assemble interpolation matrix.\n");
				return false;
			}

			#ifdef UG_PARALLEL
				m_matrix.set_storage_type(PST_CONSISTENT);
			#endif

			m_bInit = true;

			return true;
		}

		virtual bool init(const vector_type& u)
		{
			return init();
		}

	// 	Project uFine to uCoarse; uCoarse = P(uFine);
		virtual bool apply(vector_type& uCoarseOut, const vector_type& uFineIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
			{
				UG_LOG("ERROR in 'ProjectionOperator::apply':Operator not initialized.\n");
				return false;
			}

		//	Some Assertions
			UG_ASSERT(uCoarseOut.size() == m_matrix.num_rows(),	"Vector [size= " << uCoarseOut.size() << "] and Row [size= " << m_matrix.num_rows() <<"] sizes have to match!");
			UG_ASSERT(uFineIn.size() == m_matrix.num_cols(),	"Vector [size= " << uFineIn.size() << "] and Column [size= " << m_matrix.num_cols() <<"] sizes have to match!");

		//	Apply matrix
			if(!m_matrix.apply(uCoarseOut, uFineIn))
			{
				UG_LOG("ERROR in 'P1ProjectionOperator::apply': Cannot apply matrix.\n");
				return false;
			}

		//	we're done
			return true;
		}

	// 	Apply sub not implemented
		virtual bool apply_sub(vector_type& u, const vector_type& v)
		{
			UG_ASSERT(0, "Not Implemented.");
			return false;
		}

	//	Destructor
		~P1ProjectionOperator()
		{
			m_matrix.destroy();
		}

		virtual IProjectionOperator<vector_type, vector_type>* clone()
		{
			P1ProjectionOperator* op = new P1ProjectionOperator;
			op->set_approximation_space(*m_pApproximationSpace);
			return op;
		}

	protected:
		matrix_type m_matrix;

		TApproximationSpace* m_pApproximationSpace;

		size_t m_fineLevel;
		size_t m_coarseLevel;

		bool m_bInit;
};


}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__ */
