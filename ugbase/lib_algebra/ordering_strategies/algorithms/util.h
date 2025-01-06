
#ifndef UG_BASE_LIB_ALGEBRA_ORDERING_STRATEGIES_ALGORITHMS_UTIL_H
#define UG_BASE_LIB_ALGEBRA_ORDERING_STRATEGIES_ALGORITHMS_UTIL_H

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

}

#endif
