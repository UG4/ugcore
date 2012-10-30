// 23.11.2011 (m,d,y)
 
#include "registry_util.h"

namespace ug{

/*
 * 
 * For compatibility to other languages like VRL we need some mappings of
 * classes/interfaces and functions
 * 
 * F_ function
 * C_ class
 * I_ interface
 * const__ for const methods
 * 
 * Otherwise we will get name collisions
 * At the moment, __ is reserved for further extensions of this scheme,
 * but not if an identifier is starting with __, like LUAs
 * __index and __tostring methods.
 * 
 */
	
bool IsValidRegistryIdentifier(const std::string& name) {
	return StartsWith(name,"__") || (!Contains(name,"__") &&
			!StartsWith(name, "F_") &&
			!StartsWith(name, "C_") &&
			!StartsWith(name, "I_") &&
			name!="constructor");
}

std::string GetRegistryIdentifierMessage() {
	return "Identifier names must not start with"
			" 'F_', 'C_' or 'I_' and must not contain"
			" '__' (double underscore) nor be equal to 'constructor'.";
}

}// end of namespace
