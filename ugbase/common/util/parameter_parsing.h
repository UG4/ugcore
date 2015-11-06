/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
