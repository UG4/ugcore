

#ifndef __H__VRL_BRIDGE__
#define __H__VRL_BRIDGE__

#include <string>
#include "registry/registry.h"

namespace ug {
namespace vrl {

void RegisterVRLFunctionality(ug::bridge::Registry& reg, std::string grp);

}}

#endif

