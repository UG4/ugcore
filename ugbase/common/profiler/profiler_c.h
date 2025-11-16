/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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


#ifndef __H__UG__COMMON__PROFILER_C__
#define __H__UG__COMMON__PROFILER_C__

#ifdef __H__UG__COMMON__PROFILER__
#error "can't use profiler_c.h and profiler.h"
#endif

#ifdef UG_PROFILER


#ifdef UG_PROFILER_SHINY

//	this is just a wrapper-include for the shiny-profiler by Aidin Abedi
	#define SHINY_PROFILER TRUE
	#include "src/ShinyManager.h"
	#include "src/ShinyNode.h"


	#define C_PROFILE_BEGIN_PARAMS(id, name, group, file, line)			\
															\
		static Shiny::ProfileZone __ShinyZone_##id = {		\
			nullptr, Shiny::ProfileZone::STATE_HIDDEN, name, \
			group, file, line,	\
			{ { 0, 0 }, { 0, 0 }, { 0, 0 } }				\
		};													\
		{													\
			static Shiny::ProfileNodeCache cache =			\
				&Shiny::ProfileNode::_dummy;				\
															\
			Shiny::ProfileManager::instance._beginNode(&cache, &__ShinyZone_##id);\
		}


	#define C_PROFILE_BEGIN_GROUP(name, group)						\
		C_PROFILE_BEGIN_PARAMS(apn_##name, #name, group, __FILE__, __LINE__)

	#define C_PROFILE_BEGIN(name)						\
		C_PROFILE_BEGIN_GROUP(name, nullptr)

	#define C_PROFILE_END()														\
		Shiny::ProfileManager::instance._endCurNode()

	#define C_PROFILE_FUNC_BEGIN() \
		C_PROFILE_BEGIN(__FUNCTION__)
	#define PROFILE_FUNC_GROUP_BEGIN(groups) \
		C_PROFILE_BEGIN_GROUP(__FUNCTION__, groups)

#else
#error "not defined for C"
#endif // UG_PROFILER_SHINY

#else
	#define C_PROFILE_BEGIN(name)
	#define C_PROFILE_BEGIN_GROUP(name, groups)
	#define C_PROFILE_END()
	#define C_PROFILE_FUNC_BEGIN()
	#define C_PROFILE_FUNC_GROUP_BEGIN(groups)

#endif // UG_PROFILER

#endif	// __H__UG__COMMON__PROFILER__
