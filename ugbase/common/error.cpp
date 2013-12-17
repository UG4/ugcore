/*
 * error.cpp
 *
 *  Created on: 09.12.2013
 *      Author: mrupp
 */

#include "error.h"
namespace ug{

UGError::UGError(const std::string &msg, std::bad_alloc &ex, const char *file, const unsigned long line)
{
	std::stringstream ss;
	ss << "bad_alloc caught: " << ex.what() << "\nThis is mostly caused by an OUT OF MEMORY - condition.\n"\
			 	 << "You might have to reduce your problem size, go parallel or increase memory.";
	push_msg(ss.str(), file, line);
	push_msg(msg, file, line);

}

UGError::UGError(const std::string &msg, std::bad_cast &ex, const char *file, const unsigned long line)
{
	std::stringstream ss;
	ss << "bad_cast caught: " << ex.what() << "\nThis is caused by a Casting error in classes.";
	push_msg(ss.str(), file, line);
	push_msg(msg, file, line);
}

UGError::UGError(const std::string &msg, std::exception &ex, const char *file, const unsigned long line)
{
	std::stringstream ss;
	ss << "std::exception caught: " << ex.what();
	push_msg(ss.str(), file, line);
	push_msg(msg, file, line);
}
}
