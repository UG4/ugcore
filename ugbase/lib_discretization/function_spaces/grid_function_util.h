/*
 * grid_function_util.h
 *
 *  Created on: 17.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_UTIL__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_UTIL__

#include "./grid_function.h"
#include "lib_algebra/lib_algebra.h"
#include <vector>

namespace ug{


template <typename TFunction>
void ExtractPositions(const TFunction &u, std::vector<MathVector<TFunction::domain_type::dim> > &positions)
{
	typedef typename TFunction::domain_type domain_type;
	const typename domain_type::position_accessor_type& aaPos = u.get_approximation_space().get_domain().get_position_accessor();

	int nr = u.num_dofs();
	positions.resize(nr);

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		for(geometry_traits<VertexBase>::const_iterator iter = u.template begin<VertexBase>(si); iter != u.template end<VertexBase>(si); ++iter)
		{
			VertexBase* v = *iter;

			typename TFunction::algebra_index_vector_type ind;

			u.get_inner_algebra_indices(v, ind);

			for(size_t i = 0; i < ind.size(); ++i)
			{
				const size_t index = ind[i];
				positions[index] = aaPos[v];
			}
		}
	}
}

template <class TFunction>
void WriteMatrixToConnectionViewer(const char *filename, const typename TFunction::algebra_type::matrix_type &A, const TFunction &u)
{
	const static int dim = TFunction::domain_type::dim;
	std::vector<MathVector<dim> > positions;

	// get positions of vertices
	ExtractPositions(u, positions);

	// write matrix
	WriteMatrixToConnectionViewer(filename, A, &positions[0], dim);
}

template <class TFunction>
void WriteVectorToConnectionViewer(const char *filename, const typename TFunction::algebra_type::vector_type &b, const TFunction &u)
{
	const static int dim = TFunction::domain_type::dim;
	std::vector<MathVector<dim> > positions;

	// get positions of vertices
	ExtractPositions(u, positions);

	// write matrix
	WriteVectorToConnectionViewer(filename, b, &positions[0], dim);
}
} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_UTIL__ */
