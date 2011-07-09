/*
 * prolongation_operator.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_discretization/spatial_discretization/post_process/post_process_interface.h"

#ifdef UG_PARALLEL
#include "lib_discretization/parallelization/parallelization_util.h"
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
template <typename TApproximationSpace, typename TAlgebra>
bool AssembleVertexProlongation(typename TAlgebra::matrix_type& mat,
								TApproximationSpace& approxSpace,
								size_t coarseLevel, size_t fineLevel,
								std::vector<bool>& vIsRestricted)
{
//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& coarseDoFDistr
			= approxSpace.get_level_dof_distribution(coarseLevel);
	const typename TApproximationSpace::dof_distribution_type& fineDoFDistr
			= approxSpace.get_level_dof_distribution(fineLevel);

// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDoFDistr.num_fct(); ++fct)
		if(fineDoFDistr.local_shape_function_set_id(fct)
				!= LSFSID(LSFSID::LAGRANGE, 1))
	{
		UG_LOG("ERROR in 'AssembleVertexProlongation':"
				"Interpolation only implemented for Lagrange P1 functions.\n");
		return false;
	}

//  get subsethandler and grid
	MultiGrid& grid = approxSpace.get_domain().get_grid();

//  get number of dofs on different levels
	const size_t numFineDoFs = fineDoFDistr.num_dofs();
	const size_t numCoarseDoFs = coarseDoFDistr.num_dofs();

//	check if grid distribution has dofs, otherwise skip creation since father
//	elements may not exist in parallel.
	if(numFineDoFs == 0 || numCoarseDoFs == 0)
		return true;

//  resize matrix
	if(!mat.resize(numFineDoFs, numCoarseDoFs))
	{
		UG_LOG("ERROR in 'AssembleVertexProlongation':"
				"Cannot resize Interpolation Matrix.\n");
		return false;
	}

//	clear restricted vector
	vIsRestricted.clear(); vIsRestricted.resize(numCoarseDoFs, false);

	typename TApproximationSpace::dof_distribution_type::multi_index_vector_type fineMultInd;
	typename TApproximationSpace::dof_distribution_type::multi_index_vector_type coarseMultInd;

//  iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	for(int si = 0; si < fineDoFDistr.num_subsets(); ++si)
	{
		iterBegin = fineDoFDistr.template begin<VertexBase>(si);
		iterEnd = fineDoFDistr.template end<VertexBase>(si);

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

			for(size_t fct = 0; fct < fineDoFDistr.num_fct(); fct++)
			{
				if(!fineDoFDistr.is_def_in_subset(fct, si)) continue;

			//  get global indices
				if(fineDoFDistr.inner_multi_indices(*iter, fct, fineMultInd) != 1)
					return false;

			//  Check if father is Vertex
				if(vert != NULL)
				{
				//  get global indices
					if(coarseDoFDistr.inner_multi_indices(vert, fct, coarseMultInd) != 1)
						return false;

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
						if(coarseDoFDistr.inner_multi_indices(v, fct, coarseMultInd) != 1)
							return false;

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
						if(coarseDoFDistr.inner_multi_indices(v, fct, coarseMultInd) != 1)
							return false;

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
						if(coarseDoFDistr.inner_multi_indices(v, fct, coarseMultInd) != 1)
							return false;

						vIsRestricted[coarseMultInd[0][0]] = true;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.125;
					}
					continue;
				}

				UG_LOG("ERROR in 'AssembleVertexProlongation':"
						" Element Father not detected.\n");
				return false;
			}
		}
	}

//	we're done
	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
class P1ProlongationOperator :
	virtual public IProlongationOperator<	typename TAlgebra::vector_type,
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
		// Transfer Operator acts on level -> level + 1
		P1ProlongationOperator() :
			m_pApproxSpace(NULL), m_fineLevel(0), m_coarseLevel(0),
			m_bInit(false), m_dampRes(1.0)
		{
			m_vPostProcess.clear();
		};

	//	set dirichlet values
		// todo: This should be a IPostProcess, indicating dirichlet value, only
		void set_dirichlet_post_process(IPostProcess<typename dof_distribution_type::implementation_type, algebra_type>& pp)
		{
			m_vPostProcess.push_back(&pp);
		}

	//	Set approximation space
		void set_approximation_space(approximation_space_type& approxSpace)
		{
			m_pApproxSpace = &approxSpace;
		}

	//	set interpolation damping
		void set_restriction_damping(number damp)
		{
			m_dampRes = damp;
		}

	//	Set levels
		virtual bool set_levels(size_t coarseLevel, size_t fineLevel)
		{
			m_fineLevel = fineLevel;
			m_coarseLevel = coarseLevel;
			if(m_fineLevel - m_coarseLevel != 1)
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::set_levels':"
						" Can only project between successive level.\n");
				return false;
			}
			return true;
		}

	public:
		virtual bool init(const vector_type& u)
		{
			return init();
		}

		virtual bool init()
		{
			if(m_pApproxSpace == NULL)
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::init':"
						"Approximation Space not set. Cannot init Projection.\n");
				return false;
			}

			if(m_fineLevel - m_coarseLevel != 1)
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::init':"
						" Can only project between successiv level.\n");
				return false;
			}

			m_matrix.resize(0,0);

			if(!AssembleVertexProlongation<approximation_space_type, algebra_type>
				(m_matrix, *m_pApproxSpace, m_coarseLevel, m_fineLevel, m_vIsRestricted))
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::init':"
						"Cannot assemble interpolation matrix.\n");
				return false;
			}

			#ifdef UG_PARALLEL
				m_matrix.set_storage_type(PST_CONSISTENT);
			#endif

			m_bInit = true;

			return true;
		}

		// apply Operator, interpolate function
		virtual bool apply(vector_type& uFineOut, const vector_type& uCoarseIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::apply':"
						" Operator not initialized.\n");
				return false;
			}

		//	Some Assertions
			UG_ASSERT(uFineOut.size() == m_matrix.num_rows(),
			          	  "Vector and Row sizes have to match!");
			UG_ASSERT(uCoarseIn.size() == m_matrix.num_cols(),
			          	  "Vector and Column sizes have to match!");

		//	Apply Matrix
			if(!m_matrix.apply(uFineOut, uCoarseIn))
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::apply': "
						"Cannot apply matrix. "
#ifdef UG_PARALLEL
						"(Type uCoarse = " <<uCoarseIn.get_storage_mask() <<".\n");
#else
				".\n");
#endif
				return false;
			}

		//	Set dirichlet nodes to zero again
		//	todo: We could handle this by eliminating dirichlet rows as well
			for(size_t i = 0; i < m_vPostProcess.size(); ++i)
			{
				const dof_distribution_type& dofDistr =
						m_pApproxSpace->get_level_dof_distribution(m_fineLevel);
				if(m_vPostProcess[i]->post_process_defect(uFineOut, uFineOut, dofDistr) != true)
				{
					UG_LOG("ERROR in 'P1ProlongationOperator::apply': "
							"Error while setting dirichlet defect nr " << i << " to zero.\n");
					return false;
				}
			}

		//	we're done
			return true;
		}

		// apply transposed Operator, restrict function
		bool apply_transposed(vector_type& uCoarseOut, const vector_type& uFineIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::apply_transposed':"
						"Operator not initialized.\n");
				return false;
			}

			vector_type	uTmp; uTmp.resize(uCoarseOut.size());

		//	Some Assertions
			UG_ASSERT(uFineIn.size() == m_matrix.num_rows(),
			          	  "Vector and Row sizes have to match!");
			UG_ASSERT(uCoarseOut.size() == m_matrix.num_cols(),
			          	  "Vector and Column sizes have to match!");

		//	Apply transposed matrix
			if(!m_matrix.apply_transposed(uTmp, uFineIn))
			{
				UG_LOG("ERROR in 'P1ProlongationOperator::apply_transposed':"
						" Cannot apply transposed matrix.\n");
				return false;
			}

			uTmp *= m_dampRes;

		//	Copy only restricted values
			for(size_t i = 0; i < uTmp.size(); ++i)
				if(m_vIsRestricted[i])
					uCoarseOut[i] = uTmp[i];

		//	Set dirichlet nodes to zero again
		//	todo: We could handle this by eliminating dirichlet columns as well
			for(size_t i = 0; i < m_vPostProcess.size(); ++i)
			{
				const dof_distribution_type& dofDistr =
						m_pApproxSpace->get_level_dof_distribution(m_coarseLevel);
				if(m_vPostProcess[i]->post_process_defect(uCoarseOut, uCoarseOut, dofDistr) != true)
				{
					UG_LOG("ERROR in 'ProjectionOperator::apply_transposed': "
							"Error while setting dirichlet defect nr " << i << " to zero.\n");
					return false;
				}
			}

		//	we're done
			return true;
		}

		// apply Operator, i.e. v := v - L(u);
		virtual bool apply_sub(vector_type& u, const vector_type& v)
		{
			UG_ASSERT(0, "Not Implemented.");
			return true;
		}

		virtual IProlongationOperator<vector_type, vector_type>* clone()
		{
			P1ProlongationOperator* op = new P1ProlongationOperator;
			op->set_approximation_space(*m_pApproxSpace);
			for(size_t i = 0; i < m_vPostProcess.size(); ++i)
				op->set_dirichlet_post_process(*m_vPostProcess[i]);
			op->set_restriction_damping(m_dampRes);
			return op;
		}

		~P1ProlongationOperator()
		{
		}



	protected:
		matrix_type m_matrix;

		std::vector<IPostProcess<typename dof_distribution_type::implementation_type, algebra_type>*> m_vPostProcess;
		TApproximationSpace* m_pApproxSpace;
		size_t m_fineLevel;
		size_t m_coarseLevel;

		std::vector<bool> m_vIsRestricted;

		bool m_bInit;
		number m_dampRes;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__ */
