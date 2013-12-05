//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#ifndef __H__UG_BRIDGE__UG_BRIDGE__
#define __H__UG_BRIDGE__UG_BRIDGE__

#include <string>
#include <sstream>
#include "registry/registry.h"
#include "lib_algebra/algebra_type.h"
#include "common/ug_config.h"

namespace ug{
namespace bridge{

/**
 * \defgroup bridge Bridge
 * \ingroup ugbase
 * \{
 */

/// string for ug4 group
extern const char* UG4_GRP;

///	returns the default registry used in ug
UG_API Registry & GetUGRegistry();

///	Sets the default classes of class-groups based on a tags using default DoFManager
UG_API void InitUG(int dim, const AlgebraType& algebraType);


/// calls RegisterStandardInterfaces and LoadPlugins if UG_PLUGINS is defined
UG_API void InitBridge();

///	registers all standard interfaces.
/**	This method is called by the constructor of Registry automatically.
 *	You don't have to call it yourself!
 */
UG_API void RegisterStandardBridges(Registry& reg, std::string grp = UG4_GRP);

// end group bridge
/// \}

}//	end bridge
}//	end ug

#include "suffix_tag.h"

#endif /* __H__UG_BRIDGE__UG_BRIDGE__ */
