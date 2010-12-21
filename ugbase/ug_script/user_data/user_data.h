
#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__

namespace ug
{
namespace bridge
{

void RegisterLuaUserData(Registry& reg, const char* parentGroup);

void RegisterLuaBoundaryNumber(Registry& reg, const char* parentGroup);

} // end namepace bridge
} // end namespace ug

#endif /* __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__ */
