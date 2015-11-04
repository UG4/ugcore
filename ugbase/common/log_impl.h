/*
 * log_impl.h
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel, sebastianreiter
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
    return i > 0 ? (unsigned int) log10 ((double) i) + 1 : 1;
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

#endif /* __H__UG__COMMON__LOG_IMPL__ */
