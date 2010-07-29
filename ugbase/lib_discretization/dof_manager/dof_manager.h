/*
 * dof_manager.h
 *
 *  Created on: 05.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER__

#include "./mg_dof_manager.h"
#include "./dof_distribution.h"
#include "./mg_dof_manager.h"
#include "./function_pattern.h"

#include "./p1conform/p1conform.h"

namespace ug{

enum DoFManagerGroupStrategy {
	DMGS_INVALID = -1,
	DMGS_NO_GROUPING = 0,
	DMGS_GROUP_ALL,
	DMGS_NUM_STRATEGIES
};

}

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER__ */
