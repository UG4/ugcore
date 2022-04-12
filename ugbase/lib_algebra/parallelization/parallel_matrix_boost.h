
#ifndef PARALLEL_MATRIX_BOOST_H
#define PARALLEL_MATRIX_BOOST_H
#include "lib_algebra/parallelization/parallel_matrix.h"
#include "lib_algebra/small_algebra/storage/fixed_array.h"


namespace boost{

template <class T> struct graph_traits<ug::ParallelMatrix<ug::SparseMatrix<T>>>{
	typedef int vertex_descriptor;
	typedef SM_edge<T> edge_descriptor;
	typedef directed_tag directed_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef SM_out_edge_iterator<T> out_edge_iterator;
	typedef SM_adjacency_iterator<T> adjacency_iterator;
	//typedef typename ug::SparseMatrix<T>::const_row_iterator adjacency_iterator;
};

template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::ParallelMatrix<ug::SparseMatrix<T>> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

} // boost

namespace ug {

using boost::counting_iterator;

// duplicate? why does it have to be here?
template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::ParallelMatrix<ug::SparseMatrix<T>> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<class T>
inline std::pair<boost::SM_out_edge_iterator<T>, boost::SM_out_edge_iterator<T>>
					out_edges(size_t v, ug::ParallelMatrix<ug::SparseMatrix<T>> const& g)
{
	typedef boost::SM_out_edge_iterator<T> Iter;
   auto a = boost::adjacent_vertices(v, g);
	return std::make_pair(Iter(v, a.first), Iter(v, a.second));
}

} // ug

#endif // guard
