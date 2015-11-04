// created by Sebastian Reiter
// s.b.reiter@gmail.com
// oct 2013

#ifndef __H__UG__hash_function__
#define __H__UG__hash_function__

#include "common/types.h"

namespace ug{

/// \addtogroup ugbase_common_util
/// \{

///	The hashing method can be specialized for different types.
/**	A default implementation exists, which casts each key
 * to a unsigned long.
 * \{
 */
template <typename TKey> size_t hash_key(const TKey& key);

template <typename TKey> size_t hash_key(const TKey& key)
{
	return static_cast<size_t>(key);
}

/** \} */
/** \} */

}// end of namespace

#endif
