/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
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

#ifndef __H__UG__COMMON__LOG_IMPL__
#define __H__UG__COMMON__LOG_IMPL__

#include <iostream>
#include <string>  // added for 'string'
#include <sstream>  // added for 'stringstream'
#include <iomanip> // added for 'setprecision()'
#include <cmath>

namespace ug{

inline std::ostream& LogAssistant::
debug_logger()
{
	return std::cout;
}

inline std::ostream& LogAssistant::
logger()
{
	return std::cout;
}

inline std::ostream& LogAssistant::
error_logger()
{
	return m_errStream;
}

inline int LogAssistant::
get_output_process()
{
	return m_outputProc;
}

inline LogAssistant&
GetLogAssistant()
{
	return LogAssistant::instance();
}

inline unsigned int
GetNumberOfDigits (uint64_t i)
{
    return i > 0 ? static_cast<unsigned int>(log10((double) i)) + 1 : 1;
}

inline std::string
ConvertNumber (uint64_t size, unsigned int width, unsigned int numDisplayedDigits)
{
	std::stringstream ss;

	if (GetNumberOfDigits(size) > width) {
		if (size >= UNIT_EXA) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_EXA)
			   << " Ei";

		} else if (size >= UNIT_PETA) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_PETA)
			   << " Pi";

		} else if (size >= UNIT_TERA) {
			
			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_TERA)
			   << " Ti";

		} else if (size >= UNIT_GIGA) {
			
			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_GIGA)
			   << " Gi";

		} else if (size >= UNIT_MEGA) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_MEGA)
			   << " Mi";

		} else if (size >= UNIT_KILO) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_KILO)
			   << " Ki";

		}
	} else {

		ss << size;
	}

	return(ss.str());
}

inline std::string
ConvertNumberSI (uint64_t size, unsigned int width, unsigned int numDisplayedDigits)
{
	std::stringstream ss;

	if (GetNumberOfDigits(size) > width) {
		if (size >= UNIT_EXA_SI) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_EXA_SI)
			   << " E";

		} else if (size >= UNIT_PETA_SI) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_PETA_SI)
			   << " P";

		} else if (size >= UNIT_TERA_SI) {
			
			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_TERA_SI)
			   << " T";

		} else if (size >= UNIT_GIGA_SI) {
			
			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_GIGA_SI)
			   << " G";

		} else if (size >= UNIT_MEGA_SI) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_MEGA_SI)
			   << " M";

		} else if (size >= UNIT_KILO_SI) {

			ss << std::setprecision(numDisplayedDigits)
			   << (((double)size) / UNIT_KILO_SI)
			   << " K";

		}
	} else {

		ss << size;
	}

	return(ss.str());
}

} // end namespace ug

#endif