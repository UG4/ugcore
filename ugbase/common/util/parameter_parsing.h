// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)

#ifndef __H__UG__parameter_parsing__
#define __H__UG__parameter_parsing__

#include "common/ug_config.h"

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns its position in argv.
 *	If the parameter is not contained in argv, -1 is returned.
 */
UG_API int GetParamIndex(const char* param, int argc, const char * const * argv);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns true if it is found.
 */
UG_API bool FindParam(const char* param, int argc, const char * const * argv);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter, and converts the
 *	associated value to an integer. Returns true if the parameter was
 *	found, false if not.
 */
UG_API bool ParamToInt(int& iOut, const char* param, int argc, const char * const * argv);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter, and converts the
 *	associated value to a double value. Returns true if the parameter was
 *	found, false if not.
 */
UG_API bool ParamToDouble(double &dOut, const char *param, int argc, const char * const * argv);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter, and returns the
 *	associated string (the argv directly following param).
 *	associated Returns true if the parameter was found, false if not.
 *
 *	Please note that spaces may not be contained in the associated string.
 */
UG_API bool ParamToString(const char** strOut, const char* param, int argc, const char * const * argv);



UG_API double ParamToDouble(const char *param, int argc, const char * const * argv, double dDefault);
UG_API int ParamToInt(const char *param, int argc, const char * const * argv, int iDefault);

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif
