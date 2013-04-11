//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d09

#ifndef __H__LIB_GRID__TYPES__
#define __H__LIB_GRID__TYPES__

#include <cstddef>

//	stdint.h is included with most compilers.
//	If stdint.h is not available (e.g. in Microsoft Visual C++.)
//	you may either download pstint.h or use boost/cstdint.h
#include <stdint.h>
//#include "pstdint.h"

/// \addtogroup ugbase_common
/// \{

///brief Include C99 int-types into the lg namespace.
namespace ugtypes
{
	using ::int8_t;
	using ::int16_t;
	using ::int32_t;
	using ::int64_t;
	using ::int_fast8_t;
	using ::int_fast16_t;
	using ::int_fast32_t;
	using ::int_fast64_t;
	using ::int_least8_t;
	using ::int_least16_t;
	using ::int_least32_t;
	using ::int_least64_t;
	using ::intmax_t;
	using ::uint8_t;
	using ::uint16_t;
	using ::uint32_t;
	using ::uint64_t;
	using ::uint_fast8_t;
	using ::uint_fast16_t;
	using ::uint_fast32_t;
	using ::uint_fast64_t;
	using ::uint_least8_t;
	using ::uint_least16_t;
	using ::uint_least32_t;
	using ::uint_least64_t;
	using ::uintmax_t;
};

/*
#include "boost/cstdint.hpp"
namespace ugtypes
{
	using boost::int8_t;
	using boost::int16_t;
	using boost::int32_t;
	using boost::int64_t;
	using boost::int_fast8_t;
	using boost::int_fast16_t;
	using boost::int_fast32_t;
	using boost::int_fast64_t;
	using boost::int_least8_t;
	using boost::int_least16_t;
	using boost::int_least32_t;
	using boost::int_least64_t;
	using boost::intmax_t;
	using boost::uint8_t;
	using boost::uint16_t;
	using boost::uint32_t;
	using boost::uint64_t;
	using boost::uint_fast8_t;
	using boost::uint_fast16_t;
	using boost::uint_fast32_t;
	using boost::uint_fast64_t;
	using boost::uint_least8_t;
	using boost::uint_least16_t;
	using boost::uint_least32_t;
	using boost::uint_least64_t;
	using boost::uintmax_t;
};
*/

using std::size_t;

typedef unsigned char byte;
typedef unsigned int uint;

typedef ugtypes::uint32_t uint32;
typedef ugtypes::uint64_t uint64;
typedef ugtypes::int32_t int32;
typedef ugtypes::int64_t int64;

#ifdef UG_SINGLE_PRECISION
	typedef float number;
#else
	typedef double number;
#endif

// end group ugbase_common
/// \}

#endif
