/**
 * \file ugbase.h
 */

#ifndef __H__UG__UGBASE__
#define __H__UG__UGBASE__


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
#include "lib_disc/lib_discretization.h"

// node_tree
#include "common/node_tree/node_tree.h"

// ug_bridge
#include "ug_bridge/ug_bridge.h"

// ug_script
#include "bindings/lua/ug_script.h"

// pcl (parallel only)
#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif

#endif
