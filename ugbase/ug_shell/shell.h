
#ifndef __UG__UG_SHELL_H_
#define __UG__UG_SHELL_H_

#include "bindings/lua/lua_debug.h"


namespace ug{
namespace bridge
{

void InitShell();
int RunShell(const char *prompt=NULL);
script::debug_return DebugShell();

}
}
#endif /* __UG__UG_SHELL_H_ */
