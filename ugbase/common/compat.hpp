#ifndef UG_COMMON_COMPATIBILITY_HPP
#define UG_COMMON_COMPATIBILITY_HPP


#include <algorithm>

#if defined(__cpp_lib_ranges)
  #include <ranges>
  #define HAS_RANGES 1
#else
  #define HAS_RANGES 0
#endif



#if defined(__cpp_lib_ranges_contains)
  #define HAS_RANGES_CONTAINS 1
#else
  #define HAS_RANGES_CONTAINS 0
#endif




#if HAS_RANGES_CONTAINS
#define ug_compat_contains(c, v) std::ranges::contains((c), (v))
#elif HAS_RANGES
#define ug_compat_contains(c, v) (std::ranges::find((c), (v)) != std::ranges::end(c))
#else
#define ug_compat_contains(c, v) (std::find((c).begin(), (c).end(), (v)) != (c).end())
#endif

#endif

// erase , remove
// find
// transform