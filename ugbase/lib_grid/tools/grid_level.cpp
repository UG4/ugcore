/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "grid_level.h"

#include <sstream>
#include <iomanip>     // std::setw

#include "common/common.h"

namespace ug {


std::ostream& operator << (std::ostream& out, const GridLevel& v)
{
	if(v.is_surface()) out << "(surf, ";
	else if(v.is_level()) out << "(lev,  ";
	else UG_THROW("GridLevel: type of GridLevel not found.");

	if(v.top()) out << "top";
	else out << std::setw(3) << v.level();

	if(v.ghosts() == true) out << ", g)";
	else out << ")";

	return out;
}

std::ostream& operator << (std::ostream& stream, const GridLevel::ViewType & type) {
	if ( type == GridLevel::ViewType::LEVEL) {
		stream << "ViewType::Level";
	} else {
		stream << "ViewType::Surface";
	}
	return stream;
}

std::string GridLevelAppendix(const GridLevel& gl, int minfill)
{
	std::stringstream ss; ss << std::setfill('0') << "_";

	if(gl.is_level()){
		if(gl.ghosts()) ss << "gl" << std::setw(minfill) << gl.level();
		else 			ss << "l" << std::setw(minfill) << gl.level();
	} else if (gl.is_surface()){
		if(gl.ghosts()) ss << "gsl";
		else 			ss << "sl";
		if(gl.top()){
			ss << "top";
		} else {
			ss << std::setw(minfill) << gl.level();
		}
	} else {
		UG_THROW("GMG: GridLevel not supported.")
	}

	return ss.str();

}


} // end namespace ug
