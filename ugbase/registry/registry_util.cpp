// 23.11.2011 (m,d,y)
 
#include "registry_util.h"

namespace ug{

bool IsValidRegistryIdentifier(std::string name) {
	return !Contains(name,"__") &&
			!StartsWith(name, "F_") &&
			!StartsWith(name, "C_") &&
			!StartsWith(name, "I_") &&
			name!="constructor";
}

std::string GetRegistryIdentifierMessage() {
	return "Identifier names must not start with"
			" 'F_', 'C_' or 'I_' and must not contain"
			" '__' (double underscore) nor be equal to 'constructor'.";
}

}// end of namespace
