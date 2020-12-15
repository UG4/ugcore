
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_EXECUTION_UTIL__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_EXECUTION_UTIL__

namespace ug{

#ifndef HAVE_BOOL
#define HAVE_BOOL

class BOOL{
public:
	BOOL() : value_(bool()){}
	/* explicit */ BOOL(bool const& t): value_(t) {}
	// /* explicit */ operator bool&() { return value_; }
	/* explicit */ operator bool() const { return value_; }
private:
	char value_;
};

#endif


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
