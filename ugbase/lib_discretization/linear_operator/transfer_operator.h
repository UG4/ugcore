/*
 * linear_transfer.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LINEAR_OPERATOR__TRANSFER_OPERATOR__
#define __H__LIBDISCRETIZATION__LINEAR_OPERATOR__TRANSFER_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers

namespace ug{

template <typename TDomain, typename TAlgebra, typename TDoFManager>
class TransferOperator : public LinearOperator<ApproximationSpace<TDomain, TAlgebra, TDoFManager>, ApproximationSpace<TDomain, TAlgebra, TDoFManager> > {
	public:
		// domain space
		typedef ApproximationSpace<TDomain,TAlgebra,TDoFManager> domain_type;

		// range space
		typedef ApproximationSpace<TDomain,TAlgebra,TDoFManager>  codomain_type;

		// domain function type
		typedef typename domain_type::function_type domain_function_type;

		// codomain function type
		typedef typename codomain_type::function_type codomain_function_type;

	public:
		TransferOperator(ApproximationSpace<TDomain, TAlgebra, TDoFManager>& approxSpace) : m_approxSpace(approxSpace), m_construct(true) {};

		// Transfer Operator acts on level -> level + 1
		bool init(uint level)
		{
			if(level > m_I.size()) return false;

			m_level = level;
			return true;
		}

		// prepare Operator
		bool prepare()
		{
			if(m_construct == false) return true;

			const uint num_level = m_approxSpace.num_levels();
			m_I.resize(num_level);
			for(uint l=0; l < num_level; ++l)
			{
				if(this->assemble_interpolation(*m_I[l], l) != true) return false;
			}

			return true;
		}

		// apply Operator, i.e. v = L(u);
		bool apply(const domain_function_type& u, codomain_function_type& v)
		{
			m_I[m_level]->apply(v.get_vector(m_level + 1), u.get_vector(m_level));
			return true;
		}

		// apply Operator
		bool applyTransposed(const domain_function_type& u, codomain_function_type& v)
		{
			m_I[m_level]->applyTransposed(v.get_vector(m_level), u.get_vector(m_level + 1));
			return true;
		}

	protected:

		bool assemble_interpolation(typename TAlgebra::matrix_type& mat, uint coarseLevel)
		{
			typename TDomain::grid_type& grid = m_approxSpace.get_domain().get_grid();

			uint fineLevel = coarseLevel + 1;

			const uint num_dofs_fineLevel = m_approxSpace.num_dofs(fineLevel);
			const uint num_dofs_coarseLevel = m_approxSpace.num_dofs(coarseLevel);

			if(mat.create(num_dofs_fineLevel, num_dofs_coarseLevel) != true)
			{
				UG_LOG("Cannot create Interpolation Matrix.\n");
				return false;
			}

			uint num_fct = m_approxSpace.num_fct();

			for(uint i = 0; i < num_fct; ++i)
			{
				if(m_approxSpace.get_local_shape_function_set_id(i) != LSFS_LAGRANGEP1)
				{
					UG_LOG("Interpolation only implemented for Lagrange P1 functions.");
					return false;
				}
			}

			typename TAlgebra::matrix_type::local_matrix_type val(1,1);

			typename TAlgebra::matrix_type::local_index_type coarse_ind(1);
			typename TAlgebra::matrix_type::local_index_type fine_ind(1);

			geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

			for(int subsetIndex = 0; subsetIndex < m_approxSpace.num_subsets(); ++subsetIndex)
			{
				iterBegin = m_approxSpace.template begin<Vertex>(fineLevel, subsetIndex);
				iterEnd = m_approxSpace.template end<Vertex>(fineLevel, subsetIndex);

				for(iter = iterBegin; iter != iterEnd; ++iter)
				{
					// skip boundary nodes
					if(IsBoundaryVertex2D(grid, *iter)) continue;

					// get father
					GeometricObject* geomObj = grid.get_parent(*iter);
					VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);
					Edge* edge = dynamic_cast<Edge*>(geomObj);
					Quadrilateral* quad = dynamic_cast<Quadrilateral*>(geomObj);

					for(uint fct = 0; fct < num_fct; fct++)
					{
						if(m_approxSpace.fct_is_def_in_subset(fct, subsetIndex) != true) continue;

						if(m_approxSpace.get_multi_indices_of_geom_obj(*iter, fct, fine_ind) != 1)
						{
							UG_LOG("Cannot determine fine index of node."); return false;
						}


						// Check if father is Vertex
						if(vert != NULL)
						{
							// skip boundary nodes
							if(IsBoundaryVertex2D(grid, vert)) continue;

							val(0,0) = 1.0;

							if(m_approxSpace.get_multi_indices_of_geom_obj(*iter, fct, coarse_ind) != 1)
							{
								UG_LOG("Cannot determine fine index of node."); return false;
							}

							mat.add(val, fine_ind, coarse_ind);
							continue;
						}

						if(edge != NULL)
						{
							for(int i = 0; i < 2; ++i)
							{
								vert = edge->vertex(i);
								if(IsBoundaryVertex2D(grid, vert)) continue;

								val(0,0) = 0.5;

								if(m_approxSpace.get_multi_indices_of_geom_obj(*iter, fct, coarse_ind) != 1)
								{
									UG_LOG("Cannot determine fine index of node."); return false;
								}
								mat.add(val, fine_ind, coarse_ind);
							}
							continue;
						}

						if(quad != NULL)
						{
							val(0,0) = 0.25;

							for(int i = 0; i < 4; ++i)
							{
								vert = quad->vertex(i);
								if(IsBoundaryVertex2D(grid, vert)) continue;

								if(m_approxSpace.get_multi_indices_of_geom_obj(*iter, fct, coarse_ind) != 1)
								{
									UG_LOG("Cannot determine fine index of node."); return false;
								}
								mat.add(val, fine_ind, coarse_ind);
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


	protected:
		ApproximationSpace<TDomain, TAlgebra, TDoFManager>& m_approxSpace;
		uint m_level;
		bool m_construct;
		typedef typename TAlgebra::matrix_type matrix_type;
		std::vector<matrix_type*> m_I;
};


}

#endif /* __H__LIBDISCRETIZATION__LINEAR_OPERATOR__TRANSFER_OPERATOR__ */
