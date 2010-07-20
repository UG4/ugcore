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
#include "lib_algebra/lib_algebra.h"

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
template <typename TFunction, typename TAlgebra>
bool AssembleVertexProlongation(typename TAlgebra::matrix_type& mat, IAssemble<TFunction, TAlgebra>& ass, TFunction& uCoarse, TFunction& uFine)
{
	// TODO: Currently the dirichlet values are given by IAssemble, but not even each position is checked
	//       In case,that not the whole domain is dirichlet, this is a problem
	//       One should replace IAssemble by just a Dirichlet object, indicating dirichlet bnd more flexible


	// check, that u and v are from same approximatio space
	if(&uFine.get_approximation_space() != &uCoarse.get_approximation_space())
		{UG_LOG("Interpolation only implemented for equal approximation space.\n"); return false;}

	// allow only lagrange P1 functions
	for(size_t fct = 0; fct < uFine.num_fct(); ++fct)
		if(uFine.local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
			{UG_LOG("Interpolation only implemented for Lagrange P1 functions.\n"); return false;}

	// get subsethandler and grid
	MultiGrid& grid = uFine.get_approximation_space().get_domain().get_grid();
	MGSubsetHandler& sh = uFine.get_approximation_space().get_domain().get_subset_handler();

	// get number of dofs on different levels
	const size_t numFineDoFs = uFine.num_dofs();
	const size_t numCoarseDoFs = uCoarse.num_dofs();

	// create matrix
	if(!mat.destroy())
		{UG_LOG("Cannot destroy Interpolation Matrix.\n"); return false;}
	if(!mat.create(numFineDoFs, numCoarseDoFs))
		{UG_LOG("Cannot create Interpolation Matrix.\n"); return false;}

	typename TAlgebra::matrix_type::local_matrix_type val(1,1);
	typename TAlgebra::matrix_type::local_index_type coarse_ind(1);
	typename TAlgebra::matrix_type::local_index_type fine_ind(1);

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop subsets on fine level
	for(int si = 0; si < uFine.num_subsets(); ++si)
	{
		iterBegin = uFine.template begin<Vertex>(si);
		iterEnd = uFine.template end<Vertex>(si);

		// loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);
			Edge* edge = dynamic_cast<Edge*>(geomObj);
			Quadrilateral* quad = dynamic_cast<Quadrilateral*>(geomObj);

			for(size_t fct = 0; fct < uFine.num_fct(); fct++)
			{
				if(ass.is_dirichlet(si, fct)) continue;
				if(!uFine.is_def_in_subset(fct, si)) continue;

				if(uFine.get_multi_indices_of_geom_obj(*iter, fct, fine_ind) != 1)
					{UG_LOG("Cannot determine fine index of node."); return false;}

				// Check if father is Vertex
				if(vert != NULL)
				{
					val(0,0) = 1.0;

					if(uCoarse.get_multi_indices_of_geom_obj(vert, fct, coarse_ind) != 1)
						{UG_LOG("Cannot determine coarse index of node."); return false;}

					// skip boundary nodes
					int si = sh.get_subset_index(vert);
					if(ass.is_dirichlet(si, fct)) continue;
					mat.add(val, fine_ind, coarse_ind);
					continue;
				}

				// Check if father is Edge
				if(edge != NULL)
				{
					val(0,0) = 0.5;
					for(int i = 0; i < 2; ++i)
					{
						vert = edge->vertex(i);

						if(uCoarse.get_multi_indices_of_geom_obj(vert, fct, coarse_ind) != 1)
							{UG_LOG("Cannot determine coarse index of edge."); return false;}

						// skip boundary nodes
						int si = sh.get_subset_index(vert);
						if(ass.is_dirichlet(si, fct)) continue;
						mat.add(val, fine_ind, coarse_ind);
					}
					continue;
				}

				// Check if father is Quad
				if(quad != NULL)
				{
					val(0,0) = 0.25;
					for(int i = 0; i < 4; ++i)
					{
						vert = quad->vertex(i);

						if(uCoarse.get_multi_indices_of_geom_obj(vert, fct, coarse_ind) != 1)
							{UG_LOG("Cannot determine coarse index of edge."); return false;}

						// skip boundary nodes
						int si = sh.get_subset_index(vert);
						if(ass.is_dirichlet(si, fct)) continue;

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


// TODO: Assert, that discrete function is a level grid function or define a senseful prolongation for surface functions

template <typename TDiscreteFunction>
class ProlongationOperator : public ILinearOperator<TDiscreteFunction, TDiscreteFunction> {
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
		ProlongationOperator(IAssemble<domain_function_type, algebra_type>& ass) :
			m_ass(ass)
		{};

		bool init()
		{
			return true;
		}

		// prepare Operator (u=coarse, v = fine)
		bool prepare(domain_function_type& uCoarse, codomain_function_type& uFine)
		{
			if(!AssembleVertexProlongation(m_matrix, m_ass, uCoarse, uFine))
				{UG_LOG("ERROR in 'TransferOperator::prepare(u,v)': Cannot assemble interpolation matrix.\n"); return false;}
			return true;
		}

		// apply Operator, interpolate function
		bool apply(domain_function_type& uCoarse, codomain_function_type& uFine)
		{
			// v = coarse, u = fine
			m_matrix.apply(uFine.get_vector(), uCoarse.get_vector());
#ifdef UG_PARALLEL
			uFine.copy_storage_type(uCoarse);
#endif
			return true;
		}

		// apply transposed Operator, restrict function
		bool apply_transposed(codomain_function_type& uFine, domain_function_type& uCoarse)
		{
			// v = coarse, u = fine
			m_matrix.apply_transposed(uCoarse.get_vector(), uFine.get_vector());
#ifdef UG_PARALLEL
			uCoarse.copy_storage_type(uFine);
#endif
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
		IAssemble<domain_function_type, algebra_type>& m_ass;
		matrix_type m_matrix;
};


}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__TRANSFER_OPERATOR__ */
