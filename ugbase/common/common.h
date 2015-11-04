//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#ifndef __H__UG__COMMON__COMMON__
#define __H__UG__COMMON__COMMON__

#include <iostream>
#include <fstream>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

/**
 * \defgroup ugbase_common Common
 * \ingroup ugbase
 * \brief common utilities for ug4
 * \{
 */

#include "ug_config.h"
#include "types.h"
#include "log.h"
#include "assert.h"
#include "error.h"
#include "static_assert.h"
#include "util/metaprogramming_util.h"

// depreciated, currently here for backward compatibility
#define LOG(msg) UG_LOG(msg)
#define STATIC_ASSERT(expr, msg) UG_STATIC_ASSERT(expr, msg)

////////////////////////////////////////////////////////////////////////////////
// save pointer (de-)allocation
////////////////////////////////////////////////////////////////////////////////

#define SAFE_DELETE(a)		{if(a){ delete a; a = NULL;}}
#define SAFE_RELEASE(p)		{if(p) { (p)->Release(); (p)=NULL;}}

// end group ugbase_common
/// \}

#endif /* __H__UG__COMMON__COMMON__ */
