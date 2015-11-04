
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
		ss << "bad_cast exception: " << pex->what() << "\nThis is caused by a Casting error in classes:"
				"An object is expected to be be or derive from Class A, but is not.";
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
