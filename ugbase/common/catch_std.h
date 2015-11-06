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

#ifndef __H_COMMON_CATCH_STD_H_
#define __H_COMMON_CATCH_STD_H_
#include "log.h"
#include "assert.h"

#ifdef NDEBUG
#define UG_LOG_CATCH(expr)\
		UG_LOG_ALL_PROCS(	"\nA std::exception has been thrown:\n"\
				"Description:   " << expr << "\n"\
				"File:        " << __FILE__ << "\n"\
				"Line:        " << __LINE__ << "\n\n"); \
				assert(0);
#else
	#define UG_LOG_CATCH(expr) UG_ASSERT(false, "A std::exception has been thrown:\n" << expr)
#endif

#define CATCH_STD_EXCEPTIONS()\
	catch(std::bad_alloc& ex)	{	UG_LOG_CATCH("bad_alloc caught: " << ex.what() << "\nThis is mostly caused by an OUT OF MEMORY - condition.\n"\
			 	 << "You might have to reduce your problem size, go parallel or increase memory.")	} \
	catch(std::bad_cast& ex)	{	UG_LOG_CATCH("bad_cast caught: " << ex.what() << "\nThis is caused by a Casting error in classes.")	} \
	catch(std::exception& ex)	{	UG_LOG_CATCH("std::exception caught: " << ex.what() << "\n.")	} \



#endif /* __H_COMMON_CATCH_STD_H_ */
