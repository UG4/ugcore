/*
 * linear_transfer.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__TRANSFER_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__TRANSFER_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/operator/operator_interface.h"

// library intern headers

namespace ug{


// TODO: This function should be put to an util file
/** AssembleVertexInterpolation
 *
 * This functions assembles the interpolation matrix between to
 * grid levels using only the Vertex degrees of freedom.
 *
 * \param[in] 	uCoarse			Grid function on coarse level
 * \param[in] 	uFine 			Grid function on fine level
 * \param[out]	mat 			Assembled interpolation matrix that interpolates u -> v
 *
 */
template <typename TApproximationSpace, typename TAlgebra>
bool AssembleVertexProlongation(typename TAlgebra::matrix_type& mat,
								IAssemble<typename TApproximationSpace::dof_distribution_type, TAlgebra>& ass,
								TApproximationSpace& approxSpace, size_t coarseLevel, size_t fineLevel)
{
	// TODO: Currently the dirichlet values are given by IAssemble, but not even each position is checked
	//       In case,that not the whole domain is dirichlet, this is a problem
	//       One should replace IAssemble by just a Dirichlet object, indicating dirichlet bnd more flexible

//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& coarseDoFDistr = approxSpace.get_level_dof_distribution(coarseLevel);
	const typename TApproximationSpace::dof_distribution_type& fineDoFDistr = approxSpace.get_level_dof_distribution(fineLevel);

// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDoFDistr.num_fct(); ++fct)
		if(fineDoFDistr.local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
			{UG_LOG("Interpolation only implemented for Lagrange P1 functions.\n"); return false;}

	// get subsethandler and grid
	MultiGrid& grid = approxSpace.get_domain().get_grid();
	MGSubsetHandler& sh = approxSpace.get_domain().get_subset_handler();

	// get number of dofs on different levels
	const size_t numFineDoFs = fineDoFDistr.num_dofs();
	const size_t numCoarseDoFs = coarseDoFDistr.num_dofs();

	// create matrix
	if(!mat.destroy())
		{UG_LOG("Cannot destroy Interpolation Matrix.\n"); return false;}
	if(!mat.create(numFineDoFs, numCoarseDoFs))
		{UG_LOG("Cannot create Interpolation Matrix.\n"); return false;}

	typename TApproximationSpace::dof_distribution_type::multi_index_vector_type fineMultInd;
	typename TApproximationSpace::dof_distribution_type::multi_index_vector_type coarseMultInd;

	// iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

	// loop subsets on fine level
	for(int si = 0; si < fineDoFDistr.num_subsets(); ++si)
	{
		iterBegin = fineDoFDistr.template begin<Vertex>(si);
		iterEnd = fineDoFDistr.template end<Vertex>(si);

		// loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);
			Edge* edge = dynamic_cast<Edge*>(geomObj);
			Quadrilateral* quad = dynamic_cast<Quadrilateral*>(geomObj);
			Hexahedron* hexaeder = dynamic_cast<Hexahedron*>(geomObj);

			for(size_t fct = 0; fct < fineDoFDistr.num_fct(); fct++)
			{
				if(ass.is_dirichlet(si, fct)) continue;
				if(!fineDoFDistr.is_def_in_subset(fct, si)) continue;

				// get global indices
				if(fineDoFDistr.get_inner_multi_indices(*iter, fct, fineMultInd) != 1)
					return false;

				// Check if father is Vertex
				if(vert != NULL)
				{
					// get global indices
					if(coarseDoFDistr.get_inner_multi_indices(vert, fct, coarseMultInd) != 1)
						return false;

					// skip boundary nodes
					int si = sh.get_subset_index(vert);
					if(ass.is_dirichlet(si, fct)) continue;

					BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
								fineMultInd[0][1], coarseMultInd[0][1]) = 1.0;
					continue;
				}

				// Check if father is Edge
				if(edge != NULL)
				{
					for(int i = 0; i < 2; ++i)
					{
						vert = edge->vertex(i);

						// get global indices
						if(coarseDoFDistr.get_inner_multi_indices(vert, fct, coarseMultInd) != 1)
							return false;

						// skip boundary nodes
						int si = sh.get_subset_index(vert);
						if(ass.is_dirichlet(si, fct)) continue;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.5;
					}
					continue;
				}

				// Check if father is Quad
				if(quad != NULL)
				{
					for(int i = 0; i < 4; ++i)
					{
						vert = quad->vertex(i);

						// get global indices
						if(coarseDoFDistr.get_inner_multi_indices(vert, fct, coarseMultInd) != 1)
							return false;

						// skip boundary nodes
						int si = sh.get_subset_index(vert);
						if(ass.is_dirichlet(si, fct)) continue;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.25;
					}
					continue;
				}

				// Check if father is Hexaeder
				if(hexaeder != NULL)
				{
					for(int i = 0; i < 8; ++i)
					{
						vert = hexaeder->vertex(i);

						// get global indices
						if(coarseDoFDistr.get_inner_multi_indices(vert, fct, coarseMultInd) != 1)
							return false;

						// skip boundary nodes
						int si = sh.get_subset_index(vert);
						if(ass.is_dirichlet(si, fct)) continue;

						BlockRef(mat(fineMultInd[0][0], coarseMultInd[0][0]),
									fineMultInd[0][1], coarseMultInd[0][1]) = 0.125;
					}
					continue;
				}

				UG_LOG("ERROR in assemble_interpolation: Element Father not detected." << std::endl);
				return false;
			}
		}
	}
	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
class ProlongationOperator :
	virtual public ILinearOperator<	typename TAlgebra::vector_type,
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
		ProlongationOperator(IAssemble<dof_distribution_type, algebra_type>& ass) :
			m_ass(ass), m_bInit(false)
		{};

	//	Set approximation level
		bool set_approximation_levels(approximation_space_type& approxSpace, size_t coarseLevel, size_t fineLevel)
		{
			m_pApproximationSpace = &approxSpace;
			m_fineLevel = fineLevel;
			m_coarseLevel = coarseLevel;
			return true;
		}

	public:
		virtual bool init(const vector_type& u)
		{
			return init();
		}

		virtual bool init()
		{
			if(m_pApproximationSpace == NULL)
			{
				UG_LOG("ERROR in 'ProjectionOperator::init': Approximation Space not set. Cannot init Projection.\n");
				return false;
			}

			if(m_fineLevel - m_coarseLevel != 1)
			{
				UG_LOG("ERROR in ProjectionOperator::set_approximation_levels:"
						" Can only project between successiv level.\n");
				return false;
			}

			if(!AssembleVertexProlongation(m_matrix, m_ass, *m_pApproximationSpace, m_coarseLevel, m_fineLevel))
				{UG_LOG("ERROR in 'TransferOperator::prepare(u,v)': Cannot assemble interpolation matrix.\n"); return false;}

			m_bInit = true;

			return true;
		}

		// apply Operator, interpolate function
		virtual bool apply(vector_type& uFineOut, const vector_type& uCoarseIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
			{
				UG_LOG("ERROR in 'ProjectionOperator::apply':Operator not initialized.\n");
				return false;
			}

		//	Some Assertions
			UG_ASSERT(uFineOut.size() == m_matrix.num_rows(),"Vector and Row sizes have to match!");
			UG_ASSERT(uCoarseIn.size() == m_matrix.num_cols(),"Vector and Column sizes have to match!");

			m_matrix.apply(uFineOut, uCoarseIn);
#ifdef UG_PARALLEL
			uFineOut.copy_storage_type(uCoarseIn);
#endif
			return true;
		}

		// apply transposed Operator, restrict function
		bool apply_transposed(vector_type& uCoarseOut, const vector_type& uFineIn)
		{
		//	Check, that operator is initiallized
			if(!m_bInit)
			{
				UG_LOG("ERROR in 'ProjectionOperator::apply':Operator not initialized.\n");
				return false;
			}

		//	Some Assertions
			UG_ASSERT(uFineIn.size() == m_matrix.num_rows(),"Vector and Row sizes have to match!");
			UG_ASSERT(uCoarseOut.size() == m_matrix.num_cols(),"Vector and Column sizes have to match!");

			m_matrix.apply_transposed(uCoarseOut, uFineIn);
#ifdef UG_PARALLEL
			uCoarseOut.copy_storage_type(uFineIn);
#endif
			return true;
		}

		// apply Operator, i.e. v := v - L(u);
		virtual bool apply_sub(vector_type& u, const vector_type& v)
		{
			UG_ASSERT(0, "Not Implemented.");
			return true;
		}

		~ProlongationOperator()
		{
			m_matrix.destroy();
		}



	protected:
		IAssemble<dof_distribution_type, algebra_type>& m_ass;
		matrix_type m_matrix;

		TApproximationSpace* m_pApproximationSpace;
		size_t m_fineLevel;
		size_t m_coarseLevel;

		bool m_bInit;
};


}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__TRANSFER_OPERATOR__ */
