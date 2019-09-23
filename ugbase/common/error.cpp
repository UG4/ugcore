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

#include "error.h"
#include <typeinfo>
#include <new>
#include <stdexcept>
namespace ug{

std::string ErrorStringFromStdException(const std::exception *pex)
{
	std::stringstream ss;
	if(dynamic_cast<const std::bad_alloc*>(pex))
		ss << "bad_alloc exception: " << pex->what() << "\nThis is mostly caused by an OUT OF MEMORY - condition.\n"\
			<< "You might have to reduce your problem size, go parallel or increase memory.";
	if(dynamic_cast<const std::bad_typeid*>(pex))
			ss << "bad_typeid exception: " << pex->what();
	else if(dynamic_cast<const std::bad_cast*>(pex))
		ss << "bad_cast exception: " << pex->what() << "\nThis is caused by a casting error in classes: "
			"An object is expected to be be or derive from class A, but is not.";
	else if(dynamic_cast<const std::bad_exception*>(pex))
		ss << "bad_exception exception: " << pex->what() << "\nException thrown by unexpected handler.";
	else if(dynamic_cast<const std::runtime_error*>(pex))
		ss << "runtime_error exception: " << pex->what();
	else if(dynamic_cast<const std::logic_error*>(pex))
	{
		if(dynamic_cast<const std::out_of_range*>(pex))
			ss << "out_of_range exception: " << pex->what();
		else if(dynamic_cast<const std::length_error*>(pex))
			ss << "length_error exception: " << pex->what();
		else
			ss << "logic_error exception: " << pex->what();
	}
	else
	{
		ss << "Unknown std::exception: " << pex->what();
	}
	return ss.str();
}

UGError::UGError(const std::string &msg, const std::exception &ex, const char *file, const unsigned long line)
{
	push_msg(ErrorStringFromStdException(&ex), file, line);
	push_msg(msg, file, line);
}
}
