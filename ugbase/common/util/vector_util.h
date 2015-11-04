// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.03.2012 (m,d,y)

#ifndef __H__UG__vector_util__
#define __H__UG__vector_util__

#include <vector>

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

///	Returns a pointer to the array which is managed by the std::vector.
/**	Note that this pointer may be invalidated once new elements are added
 * to the vector or old ones are removed.
 *
 * The method returns NULL if the vector is empty.
 * \{ */
template <class T>
T* GetDataPtr(std::vector<T>& v)
{
	if(v.empty())
		return NULL;
	return &v.front();
}

template <class T>
const T* GetDataPtr(const std::vector<T>& v)
{
	if(v.empty())
		return NULL;
	return &v.front();
}
/**	\} */

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif
