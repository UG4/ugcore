#ifndef __H__UG__registry_util__
#define __H__UG__registry_util__

#include "common/util/string_util.h"

namespace ug
{

/// \addtogroup registry
/// \{

/**
 * Checks whether the specified name is a valid registry identifier name.
 * <p>
 * <b>Note:</b> identifiers starting with <code>F_</code>, <code>C_</code>,
 * <code>I_</code> or containing <code>__</code>
 *  (double underscore) or being equal to 'constructor' are invalid.
 * </p>
 * @param name name to check
 * @return  <code>true</code> if the specified name is valid;
 *          <code>false</code> otherwise
 */
UG_API bool IsValidRegistryIdentifier(const std::string& name);

/**
 * Returns a message describing which registry identifiers are valid and which are not.
 * @return message string
 */
UG_API std::string GetRegistryIdentifierMessage();

// end group registry
/// \}

}//	end of namespace

#endif
