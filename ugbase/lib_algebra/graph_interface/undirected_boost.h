#ifndef UG_GRAPH_UNDIRECTED_BOOST_H
#define UG_GRAPH_UNDIRECTED_BOOST_H

#include "undirected.h"
//#include "trace.h"
#include <boost/iterator/counting_iterator.hpp>

namespace boost{

#if 0
template<class T>
class UM_adjacency_iterator : public iterator_facade<
	UM_adjacency_iterator<T>,
	typename ug::SparseMatrix<T>::const_row_iterator,
	std::input_iterator_tag,
	size_t, // <= reference
	std::intmax_t // difference_type
	 >{ //
	typedef ug::SparseMatrix<T> M;
	typedef typename M::const_row_iterator iter_t;

public:
	typedef iter_t* value_type;
	typedef intmax_t difference_type;
	typedef size_t reference;

public:
	UM_adjacency_iterator() : _base(nullptr), _end(nullptr){}
	UM_adjacency_iterator(UM_adjacency_iterator&& p) = delete;
	UM_adjacency_iterator(UM_adjacency_iterator const& p)
		: _base(p._base?(new iter_t(*p._base)) : nullptr)
		, _end(p._end?(new iter_t(*p._end)) : nullptr){
	}
	UM_adjacency_iterator(value_type p, value_type e){
		assert(p);
		_base = new iter_t(*p);
		_end = new iter_t(*e);
		skip_zeroes();
	}
	~UM_adjacency_iterator(){
		delete _base;
		delete _end;
		_base = nullptr;
		_end = nullptr;
	}
	T const& value() const {
		assert(_base);
		return _base->value();
	}
	int row() const {
		assert(_base);
		return _base->row();
	}
	size_t idx() const { untested();
		assert(_base);
		return _base->idx();
	}
	reference dereference() const {
		assert(_base);
		return _base->index(); // row?
	}
	bool operator==(const UM_adjacency_iterator& other) const {
		assert(_base);
		assert(other._base);
		return *_base == *other._base;
	}
	bool operator!=(const UM_adjacency_iterator& other) const {
		assert(_base);
		assert(other._base);
		return *_base != *other._base;
	}
	UM_adjacency_iterator operator=(const UM_adjacency_iterator& other) {
		if(other._base){
			_base = new iter_t(*other._base);
		}else{ untested();
			_base = nullptr;
		}
		if(other._end){
			_end = new iter_t(*other._end);
		}else{ untested();
			_end = nullptr;
		}
		return *this;
	}

private:
	// could use boost::filter_iterator, but does not work yet.
	void skip_zeroes(){
		while((*_base) != (*_end)){
			if(iszero(_base->value())){
				++(*_base);
			}else{
				break;
			}
		}
	}
	bool equal(UM_adjacency_iterator const& other) const { untested();
		assert(_base);
		assert(other._base);
		return *_base == *other._base;
	}
	void increment() {
		assert(_base);
		assert(_end);
		++(*_base);
		skip_zeroes();
	}
	void decrement() { untested();
		// incomplete(); // don't use. too complicated.
		assert(_base);
		--(*_base);
	}

private:
	value_type _base;
	value_type _end;
	friend class iterator_core_access;
}; // UM_adjacency_iterator
#else
template<class T>
class UM_out_edge_iterator : public iterator_facade<
	UM_out_edge_iterator<T>,
	typename ug::UndirectedMatrix<T>::adjacency_iterator,
	// bidirectional_traversal_tag, // breaks InputIterator (why?)
	std::input_iterator_tag,
	typename ug::UndirectedMatrix<T>::edge, // <= reference
	std::intmax_t // difference_type
	 >{ //
private: // types
	typedef iterator_facade<
	UM_out_edge_iterator<T>,
	typename ug::UndirectedMatrix<T>::adjacency_iterator,
	// bidirectional_traversal_tag, // breaks InputIterator (why?)
	std::input_iterator_tag,
	typename ug::UndirectedMatrix<T>::edge, // <= reference
	std::intmax_t // difference_type
	 > base_class;
public: // types
	typedef ug::SparseMatrix<T> M;
	typedef typename ug::UndirectedMatrix<T>::adjacency_iterator value_type;
	typedef intmax_t difference_type;
	typedef typename ug::UndirectedMatrix<T>::edge reference;
	typedef typename ug::UndirectedMatrix<T>::edge edge_type;

public: // construct
	explicit UM_out_edge_iterator() : base_class(), _base() {
	}
	explicit UM_out_edge_iterator(int v, value_type w)
	    : base_class(), _v(v), _base(w) {
	}
	/* explicit */ UM_out_edge_iterator(UM_out_edge_iterator const& p)
	    : base_class(p), _v(p._v), _base(p._base){
	}
	UM_out_edge_iterator(UM_out_edge_iterator&& p) = delete;
	~UM_out_edge_iterator(){
	}
	UM_out_edge_iterator& operator=(UM_out_edge_iterator const& p) {
		_v = p._v;
		_base = p._base;
		return *this;
	}
	UM_out_edge_iterator& operator=(UM_out_edge_iterator&& p) = delete;
#if 0
public: // op
	UM_out_edge_iterator operator+(int a) const{ untested();
		UM_out_edge_iterator ret(*this);
		ret.base.first += a;
		return ret;
	}
	UM_out_edge_iterator operator-(int a) const{ untested();
		UM_out_edge_iterator ret(*this);
		ret.base.first -= a;
		return ret;
	}
	difference_type operator-(UM_out_edge_iterator const& other) const{ untested();
		UM_out_edge_iterator ret(*this);
		return ret.base.first - other.base.first;
	}
#endif
private:
	reference dereference() const {
		return edge_type(_v, *_base);
	}
	bool equal(UM_out_edge_iterator const& other) const { untested();
		assert(_base.first == other._base.first);
		return _base.second == other._base.second;
	}
	void increment() { untested();
		++_base;
	}
	void decrement() { untested();
		--_base;
	}
	//			bool operator==(const UM_out_edge_iterator& other) const
	//			{ incomplete();
	//				return false;
	//			}
public:
	UM_out_edge_iterator& operator++() {
		++_base;
		return *this;
	}
	UM_out_edge_iterator operator++(int) { untested();
		UM_out_edge_iterator copy(*this);
		++*this;
		return copy;
	}
	bool operator==(const UM_out_edge_iterator& other) const {
		return _base == other._base;
	}
	bool operator!=(const UM_out_edge_iterator& other) const {
		return _base != other._base;
	}

private:
	int _v;
	value_type _base;
	friend class iterator_core_access;
}; // UM_out_edge_iterator
#endif

template <class T> struct graph_traits<ug::UndirectedMatrix<T>>{
	typedef ug::UndirectedMatrix<T> G;
	typedef int vertex_descriptor;
	typedef typename G::edge edge_descriptor;
	typedef undirected_tag directed_category;
	typedef disallow_parallel_edge_tag edge_parallel_category;
	typedef SM_traversal_tag traversal_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef UM_out_edge_iterator<T> out_edge_iterator;
	typedef typename G::adjacency_iterator adjacency_iterator;
	typedef int degree_size_type;
	typedef int vertices_size_type;
};

template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::UndirectedMatrix<T> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<class T>
int degree(int v, ug::UndirectedMatrix<T> const& M)
{
	return M.degree(v);
}

template<class T>
int in_degree(int v, ug::UndirectedMatrix<T> const& M)
{
	return degree(v, M);
}

template<class T>
int out_degree(int v, ug::UndirectedMatrix<T> const& M)
{
	return degree(v, M);
}

template<class T>
std::pair<UM_out_edge_iterator<T>, UM_out_edge_iterator<T> >
					out_edges(int v, ug::UndirectedMatrix<T> const& M)
{
	typedef typename ug::UndirectedMatrix<T>::adjacency_iterator ai;
	typedef UM_out_edge_iterator<T> ei;

	ai b(M.begin_row(v));
	ai e(M.end_row(v));

	return std::make_pair(ei(v, b), ei(v, e));
}

template<class T>
std::pair<typename ug::UndirectedMatrix<T>::adjacency_iterator,
          typename ug::UndirectedMatrix<T>::adjacency_iterator >
					adjacent_vertices(size_t v, ug::UndirectedMatrix<T> const& M)
{
	typedef typename ug::UndirectedMatrix<T>::adjacency_iterator ei;

	ei b(M.begin_row(v));
	ei e(M.end_row(v));

	return std::make_pair(b, e);
}

template<class T>
int source(typename T::edge const& e, T const&)
{
	return e.first;
}

template<class T>
int target(typename T::edge const& e, T const&)
{
	return e.second;
}

template<class T>
struct property_map<ug::UndirectedMatrix<T>, vertex_index_t>{
	typedef sparse_matrix_index_map<T> type;
	typedef type const_type;
};

template<class T>
inline typename property_map<ug::UndirectedMatrix<T>, vertex_index_t>::const_type
get(vertex_index_t, ug::UndirectedMatrix<T> const& m)
{
	return sparse_matrix_index_map<typename T::value_type>(m);
}

template<class T>
int num_vertices(ug::UndirectedMatrix<T> const& M)
{
	return M.num_rows();
}

template <class T>
class degree_property_map< ug::UndirectedMatrix<T> >
: public put_get_helper< typename graph_traits< ug::UndirectedMatrix<T> >::degree_size_type,
      degree_property_map< ug::UndirectedMatrix<T> > >
{
public:
	typedef ug::UndirectedMatrix<T> Graph;
	typedef typename graph_traits< Graph >::vertex_descriptor key_type;
	typedef typename graph_traits< Graph >::degree_size_type value_type;
	typedef value_type reference;
	typedef readable_property_map_tag category;
	degree_property_map(const Graph& g) : m_g(g) {}
	value_type operator[](const key_type& v) const { return degree(v, m_g); }

private:
	Graph const& m_g;
};

template<class T>
degree_property_map< ug::UndirectedMatrix<T> >
make_degree_map(const ug::UndirectedMatrix<T>& g)
{
	return degree_property_map< ug::UndirectedMatrix<T> >(g);
}

} // boost

namespace ug{

// adl hack/workaround?
// (required by adjacency list concept
using boost::vertices;
using boost::adjacent_vertices;
// required by vertexListGraph concept
using boost::num_vertices;
using boost::out_edges; // used in IndicenceGraph concept
using boost::out_degree; // used in IndicenceGraph concept

} // ug
#endif
