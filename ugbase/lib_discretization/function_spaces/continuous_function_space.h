
#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__CONTINUOUS_FUNCTION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__CONTINUOUS_FUNCTION_SPACE__

#include "common/common.h"

template <typename TDomain>
class ContinuousFunctionSpace {
		typedef typename TDomain::position_type position_type;

	public:
		typedef bool (function_type)(const position_type& x, number& val);

};


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__CONTINUOUS_FUNCTION_SPACE__ */
