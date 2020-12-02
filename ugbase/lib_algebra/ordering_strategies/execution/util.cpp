
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_EXECUTION_UTIL__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_EXECUTION_UTIL__

#include "../typedefs.h" //BOOL

namespace ug{

template <typename O_t>
bool is_permutation(O_t &o){
	std::vector<BOOL> container(o.size(), false);
	for(unsigned i = 0; i < o.size(); ++i){
		if(!container[o[i]]){
			container[o[i]] = true;
		}
		else{
			return false; //no doubles allowed
		}
	}
	
	for(unsigned i = 0; i < o.size(); ++i){
		if(!container[i]){
			return false;
		}
	}

	return true;
}


} //namespace

#endif //guard
