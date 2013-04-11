// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)

#ifndef __H__UG__os_info__
#define __H__UG__os_info__

#include "common/ug_config.h"

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

///	returns the standard prefix of static and dynamic libraries on this os
UG_API const char* GetDynamicLibraryPrefix();

///	returns the standard suffix of dynamic libraries on this os
UG_API const char* GetDynamicLibrarySuffix();

///	returns a string containing the path-separator for the current os
UG_API const char* GetPathSeparator();

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif
