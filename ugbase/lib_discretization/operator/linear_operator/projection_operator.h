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
#include "lib_algebra/lib_algebra.h"

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
template <typename TFunction, typename TMatrix>
bool AssembleVertexProjection(TMatrix& mat, TFunction& uFine, TFunction& uCoarse)
{
	// check, that u and v are from same approximatio space
	if(&uFine.get_approximation_space() != &uCoarse.get_approximation_space())
		{UG_LOG("Interpolation only implemented for equal approximation space.\n"); return false;}

	// allow only lagrange P1 functions
	for(size_t fct = 0; fct < uFine.num_fct(); ++fct)
		if(uFine.local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
			{UG_LOG("Interpolation only implemented for Lagrange P1 functions.\n"); return false;}

	// get MultiGrid
	MultiGrid& grid = uFine.get_approximation_space().get_domain().get_grid();

	// get number of dofs on different levels
	const size_t numFineDoFs = uFine.num_dofs();
	const size_t numCoarseDoFs = uCoarse.num_dofs();

	// create matrix
	if(!mat.destroy())
		{UG_LOG("Cannot destroy Interpolation Matrix.\n"); return false;}
	if(!mat.create(numCoarseDoFs, numFineDoFs))
		{UG_LOG("Cannot create Interpolation Matrix.\n"); return false;}

	LocalMatrix<typename TMatrix::entry_type> loc_mat;
	typename TFunction::local_index_type coarse_ind;
	typename TFunction::local_index_type fine_ind;

	// Vertex iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop subsets of fine level
	for(int si = 0; si < uFine.num_subsets(); ++si)
	{
		// prepare local indices for elem type
		if(!uCoarse.prepare_indices(ROID_VERTEX, si, fine_ind))
			{UG_LOG("Cannot prepare indices.\n"); return false;}
		// prepare local indices for elem type
		if(!uFine.prepare_indices(ROID_VERTEX, si, coarse_ind))
			{UG_LOG("Cannot prepare indices.\n"); return false;}

		iterBegin = uFine.template begin<Vertex>(si);
		iterEnd = uFine.template end<Vertex>(si);

		// loop nodes of fine subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);

			// get global indices
			uFine.update_indices(*iter, fine_ind);

			// Check if father is Vertex
			if(vert != NULL)
			{
				// get global indices
				uCoarse.update_indices(vert, coarse_ind);
			}
			else
			{UG_LOG("Vertex found, that has no father. This is impossible.\n"); return false;}

			// adjust local matrix
			loc_mat.set_indices(coarse_ind, fine_ind);
			loc_mat.set(0.0);

			// loop functions on subset
			for(size_t fct = 0; fct < uFine.num_fct(); ++fct)
			{
				if(!uFine.is_def_in_subset(fct, si)) continue;
				if(!uCoarse.is_def_in_subset(fct, si)) continue;

				// This is for nodal elements only
				loc_mat(fct, 0, fct, 0) = 1.0;

				// Add to global matrix
				mat.add(loc_mat);
			}
		}
	}
	return true;
}



template <typename TFunction>
class ProjectionOperator : public ILinearOperator<TFunction, TFunction> {
	public:
		// domain space
		typedef TFunction domain_function_type;

		// range space
		typedef TFunction  codomain_function_type;

	protected:
		typedef typename TFunction::algebra_type algebra_type;
		typedef typename algebra_type::matrix_type matrix_type;


	public:
		bool init(){return true;}

		// prepare Operator (u=coarse, v = fine)
		bool prepare(domain_function_type& uFine, codomain_function_type& uCoarse)
		{
			if(!AssembleVertexProjection(m_matrix, uFine, uCoarse))
				{UG_LOG("ERROR in 'TransferOperator::prepare(u,v)': Cannot assemble interpolation matrix.\n"); return false;}
			return true;
		}

		// Project uFine to uCoarse; uCoarse = P(uFine);
		bool apply(domain_function_type& uFine, codomain_function_type& uCoarse)
		{
			m_matrix.apply(uCoarse.get_vector(), uFine.get_vector());
#ifdef UG_PARALLEL
			uCoarse.copy_storage_type(uFine);
#endif
			return true;
		}

		/// Project uCoarse to uFine; uFine = P^{-1}(uCoarse);
		// ATTENTION: This will only affect fine nodes, that are also coarse node, thus
		//            this operation is not very senseful
		bool apply_transposed(codomain_function_type& uCoarse, domain_function_type& uFine)
		{
			m_matrix.apply_transposed(uFine.get_vector(), uCoarse.get_vector());
#ifdef UG_PARALLEL
			uFine.copy_storage_type(uCoarse);
#endif
			return true;
		}

		// apply sub not implemented
		bool apply_sub(domain_function_type& u, codomain_function_type& v)
		{
			UG_ASSERT(0, "Not Implemented.");
			return false;
		}

		~ProjectionOperator()
		{
			m_matrix.destroy();
		}

	protected:
		matrix_type m_matrix;
};


}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__ */
