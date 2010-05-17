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
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers

namespace ug{

// TODO: Assert, that discrete function is a level grid function or define a senseful prolongation for surface functions

template <typename TDiscreteFunction>
class ProlongationOperator : public IDiscreteLinearOperator<TDiscreteFunction, TDiscreteFunction> {
	public:
		// domain space
		typedef TDiscreteFunction domain_function_type;

		// range space
		typedef TDiscreteFunction  codomain_function_type;

	protected:
		typedef typename TDiscreteFunction::approximation_space_type approximation_space_type;
		typedef typename TDiscreteFunction::domain_type phys_domain_type;
		typedef typename TDiscreteFunction::algebra_type algebra_type;
		typedef typename algebra_type::matrix_type matrix_type;


	public:
		// Transfer Operator acts on level -> level + 1
		ProlongationOperator(approximation_space_type& approxSpace, IAssemble<algebra_type, domain_function_type>& ass, uint level) :
			m_approxSpace(approxSpace), m_ass(ass), m_level(level)
		{
			UG_ASSERT(level < m_approxSpace.num_levels() - 1, "Interpolation is from (level -> level+1). Requested (" << level <<" -> " << level + 1<< "), but only " << m_approxSpace.num_levels() << " levels in multigrid.\n");

		};

		bool init()
		{
			return true;
		}

		// prepare Operator (u=coarse, v = fine)
		bool prepare(domain_function_type& u, codomain_function_type& v)
		{
			UG_DLOG(LIB_DISC_TRANSFER, 3, " ---- START: 'TransferOperator::prepare'\n");
			UG_DLOG(LIB_DISC_TRANSFER, 3, " Constructing Interpolation for levels (" << m_level << " -> "<<m_level+1<<").\n");
			if(this->assemble_interpolation(m_matrix, m_level, u, v) != true)
			{
				UG_LOG("ERROR in 'TransferOperator::prepare(u,v)': Cannot assemble interpolation matrix for level transfer (" << m_level << " -> "<<m_level+1<<").\n");
				return false;
			}
			UG_DLOG(LIB_DISC_TRANSFER, 9, "Interpolation Matrix: \n" << m_matrix);
			UG_DLOG(LIB_DISC_TRANSFER, 3, " ---- END: 'TransferOperator::prepare'\n");
			return true;
		}

		// apply Operator, i.e. v = L(u);
		bool apply(domain_function_type& u, codomain_function_type& v)
		{
			UG_DLOG(LIB_DISC_TRANSFER, 5, " Matrix: (" << m_matrix.row_size() << " x " << m_matrix.col_size() << ").\n");
			UG_DLOG(LIB_DISC_TRANSFER, 5, " Fine Vector: (" <<  v.get_vector().size() << ").\n");
			UG_DLOG(LIB_DISC_TRANSFER, 5, " Coarse Vector: (" <<  u.get_vector().size() << ").\n");

			UG_DLOG(LIB_DISC_TRANSFER, 10, "Interpolation Matrix: \n" << m_matrix);

			// v = coarse, u = fine
			m_matrix.apply(u.get_vector(), v.get_vector());
			return true;
		}

		// apply Operator
		bool applyTransposed(codomain_function_type& v, domain_function_type& u)
		{
			UG_DLOG(LIB_DISC_TRANSFER, 5, " Matrix: (" << m_matrix.row_size() << " x " << m_matrix.col_size() << ").\n");
			UG_DLOG(LIB_DISC_TRANSFER, 5, " Coarse Vector: (" <<  v.get_vector().size() << ").\n");
			UG_DLOG(LIB_DISC_TRANSFER, 5, " Fine Vector: (" <<  u.get_vector().size() << ").\n");

			UG_DLOG(LIB_DISC_TRANSFER, 10, "Interpolation Matrix: \n" << m_matrix);

			// v = coarse, u = fine
			m_matrix.applyTransposed(v.get_vector(), u.get_vector());
			return true;
		}

		// apply Operator, i.e. v := v - L(u);
		bool apply_sub(domain_function_type& u, codomain_function_type& v)
		{
			UG_ASSERT(0, "Not Implemented.");
			return true;
		}

		~ProlongationOperator()
		{
			m_matrix.destroy();
		}

	protected:
		// u = coarse level, v = fine level
		bool assemble_interpolation(matrix_type& mat, uint coarseLevel, domain_function_type& u, codomain_function_type& v)
		{
			UG_DLOG(LIB_DISC_TRANSFER, 10, " Assembling Interplation by assemble_interpolation.\n");

			typename phys_domain_type::position_accessor_type aaPos = m_approxSpace.get_domain().get_position_accessor();
			typename phys_domain_type::grid_type& grid = m_approxSpace.get_domain().get_grid();
			typename phys_domain_type::position_type corner;
			number dummy_val;
			number time = 0.0;

			const uint num_dofs_fineLevel = v.num_dofs();
			const uint num_dofs_coarseLevel = u.num_dofs();

			if(mat.create(num_dofs_fineLevel, num_dofs_coarseLevel) != true)
			{
				UG_LOG("Cannot create Interpolation Matrix.\n");
				return false;
			}

			// TODO: Handle this
			uint num_fct = m_ass.num_fct();

			for(uint i = 0; i < num_fct; ++i)
			{
				if(m_approxSpace.get_local_shape_function_set_id(i) != LSFS_LAGRANGEP1)
				{
					UG_LOG("Interpolation only implemented for Lagrange P1 functions.");
					return false;
				}
			}

			typename algebra_type::matrix_type::local_matrix_type val(1,1);

			typename algebra_type::matrix_type::local_index_type coarse_ind(1);
			typename algebra_type::matrix_type::local_index_type fine_ind(1);

			geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

			// loop fine level
			for(int subsetIndex = 0; subsetIndex < m_approxSpace.num_subsets(); ++subsetIndex)
			{
				iterBegin = v.template begin<Vertex>(subsetIndex);
				iterEnd = v.template end<Vertex>(subsetIndex);

				UG_DLOG(LIB_DISC_TRANSFER, 10, " Start loop over nodes on subset "<< subsetIndex <<".\n");
				UG_DLOG(LIB_DISC_TRANSFER, 10, " Start loop over "<< num_fct <<" fundamental function.\n");

				for(iter = iterBegin; iter != iterEnd; ++iter)
				{
					// get father
					GeometricObject* geomObj = grid.get_parent(*iter);
					VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);
					Edge* edge = dynamic_cast<Edge*>(geomObj);
					Quadrilateral* quad = dynamic_cast<Quadrilateral*>(geomObj);

					for(uint fct = 0; fct < num_fct; fct++)
					{
						if(v.fct_is_def_in_subset(fct, subsetIndex) != true) continue;

						// skip boundary nodes
						corner = aaPos[*iter];
						if(IsBoundaryVertex2D(grid, *iter))
							if(m_ass.boundary_value(dummy_val, corner, fct, time))
								continue;

						if(v.get_multi_indices_of_geom_obj(*iter, fct, fine_ind) != 1)
						{
							UG_LOG("Cannot determine fine index of node."); return false;
						}

						UG_DLOG(LIB_DISC_TRANSFER, 10, " Fine Node: "<< fine_ind << " will be interpolated.\n");

						// Check if father is Vertex
						if(vert != NULL)
						{
							val(0,0) = 1.0;

							if(u.get_multi_indices_of_geom_obj(vert, fct, coarse_ind) != 1)
							{
								UG_LOG("Cannot determine fine index of node."); return false;
							}

							// skip boundary nodes
							corner = aaPos[vert];
							if(IsBoundaryVertex2D(grid, vert))
								if(m_ass.boundary_value(dummy_val, corner, fct, time))
									continue;

							UG_DLOG(LIB_DISC_TRANSFER, 10, " Coarse Node: "<< coarse_ind << " will be used.\n");

							mat.add(val, fine_ind, coarse_ind);
							continue;
						}

						if(edge != NULL)
						{
							val(0,0) = 0.5;
							for(int i = 0; i < 2; ++i)
							{
								vert = edge->vertex(i);

								if(u.get_multi_indices_of_geom_obj(vert, fct, coarse_ind) != 1)
								{
									UG_LOG("Cannot determine fine index of node."); return false;
								}
								// skip boundary nodes
								corner = aaPos[vert];
								if(IsBoundaryVertex2D(grid, vert))
									if(m_ass.boundary_value(dummy_val, corner, fct, time))
										continue;

								UG_DLOG(LIB_DISC_TRANSFER, 10, " Coarse Node: "<< coarse_ind << " will be used.\n");
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

								if(u.get_multi_indices_of_geom_obj(vert, fct, coarse_ind) != 1)
								{
									UG_LOG("Cannot determine fine index of node."); return false;
								}
								// skip boundary nodes
								corner = aaPos[vert];
								if(IsBoundaryVertex2D(grid, vert))
									if(m_ass.boundary_value(dummy_val, corner, fct, time))
										continue;

								UG_DLOG(LIB_DISC_TRANSFER, 10, " Coarse Node: "<< coarse_ind << " will be used.\n");
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
		approximation_space_type& m_approxSpace;
		IAssemble<algebra_type, domain_function_type>& m_ass;
		uint m_level;
		matrix_type m_matrix;
};


}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__TRANSFER_OPERATOR__ */
