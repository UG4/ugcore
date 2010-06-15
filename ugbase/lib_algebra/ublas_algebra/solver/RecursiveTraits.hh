#ifndef __RECURSIVE_TRAITS_HH__
#define __RECURSIVE_TRAITS_HH__

// empty dummy policy
template <typename VE, bool k>
struct recursive_traits_dummy{
  enum{depth = 0};
};

//! general type to define the level of recursion
template <typename T> 
struct recursive_traits : public
recursive_traits_dummy<T, boost::is_arithmetic<typename T::value_type>::value > {};

//dummy: recursive definition
template <typename VE>
struct recursive_traits_dummy<VE, false> 
{
  enum{depth = (recursive_traits<typename VE::value_type>::depth) + 1};
};

//dummy: specialization for leave
template <typename VE>
struct recursive_traits_dummy<VE, true> 
{
  enum{depth = 1};
};

//! general type to define the level of recursion
template <typename T> 
struct recursive_depth {
  enum {value = 
	recursive_traits_dummy<T, (boost::is_arithmetic<typename T::value_type>::value) >::depth };
};
#endif
