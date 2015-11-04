#ifndef __H__UG__UGBASE__
#define __H__UG__UGBASE__

/**
 * \defgroup ugbase ugbase
 * \brief the core functionality of ug4
 * \{
 */

/////////////////////////////////////////////////////
// This file is used to include all parts of the
// ugbase - library at once
/////////////////////////////////////////////////////

// common
#include "common/common.h"

// lib_algebra
#include "lib_algebra/lib_algebra.h"

// lib_grid
#include "lib_grid/lib_grid.h"

// lib_discretization
#include "lib_disc/lib_disc.h"

// node_tree
#include "common/node_tree/node_tree.h"

// ug_bridge
#include "bridge/bridge.h"

// ug_script
#include "bindings/lua/lua_util.h"

// pcl (parallel only)
#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif

// end group ugbase
/// \}

#endif
