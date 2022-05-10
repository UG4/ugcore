
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_UTIL__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_UTIL__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include "common/error.h" // UG_THROW

namespace ug{

#ifndef HAVE_BOOL
#define HAVE_BOOL

class BOOL{
public:
	BOOL() : value_(bool()){}
	BOOL(bool const& t): value_(t) {}
	operator bool() const { return value_; }
private:
	char value_;
};

#endif


#ifndef HAVE_IS_PERMUTATION
#define HAVE_IS_PERMUTATION

template <typename O_t>
bool is_permutation(O_t &o){
	std::vector<BOOL> container(o.size(), false);
	for(unsigned i = 0; i < o.size(); ++i){
		if(!container[o[i]]){
			container[o[i]] = true;
		}
		else{
			UG_THROW("is_permutation: multiple occurence of index, i=" << i << ", o[i]=" << o[i] << "\n");
			return false; //no doubles allowed
		}
	}

	for(unsigned i = 0; i < o.size(); ++i){
		if(!container[i]){
			UG_THROW("is_permutation: no occurence of index " << i);
			return false;
		}
	}

	return true;
}

#endif

template <typename G_t, typename M_t>
void induced_subgraph(G_t& ind_g, M_t* A, const std::vector<size_t>& inv_map){
	size_t n = A->num_rows();
	size_t k = inv_map.size();
	ind_g = G_t(k);

	std::vector<int> ind_map(n, -1);
	for(unsigned i = 0; i < k; ++i){
		ind_map[inv_map[i]] = i;
	}

	typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
	for(unsigned i = 0; i < inv_map.size(); ++i){
		for(typename M_t::row_iterator conn = A->begin_row(inv_map[i]); conn != A->end_row(inv_map[i]); ++conn){
			if(conn.value() != 0.0 && conn.index() != i){
				int idx = ind_map[conn.index()];
				if(idx >= 0){
					//boost::add_edge(i, idx, ind_g);
					boost::add_edge(idx, i, ind_g);
				}
			}
		}
	}
}

} //namespace

#endif //guard
