#ifndef __H__FUNCTION_CAST__
#define __H__FUNCTION_CAST__

/// \addtogroup ugbase_common_util
/// \{

///	Casts a function of one type to a function of another type
/**	ATTENTION: Use with care!
 * This method should only be used on functions which have similar parameters,
 * and where different parameters in TTarget and TSrc are pointers or references
 * to types which are connected through inheritance.
 *
 * The main purpose of this method is to make code readable and to allow
 * more secure compiler dependent implementations.
 *
 * \remark	The cast is currently performed through a c-style cast, since a
 *			static_cast or reinterpret_cast for the described purpose is not
 *			raises a compile error in Microsoft Visual Studio 2010.*/
template <typename TTarget, typename TSrc>
TTarget function_cast(TSrc src)
{
	return (TTarget)src;
}

// end group ugbase_common_util
/// \}

#endif
