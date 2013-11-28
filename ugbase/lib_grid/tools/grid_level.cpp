/*
 * grid_level.cpp
 *
 *  Created on: 28.11.2013
 *      Author: andreasvogel
 */

#include <sstream>
#include <iomanip>     // std::setw

#include "common/common.h"
#include "grid_level.h"

namespace ug{

std::ostream& operator<<(std::ostream& out,	const GridLevel& v)
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
