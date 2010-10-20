
#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__

namespace ug
{
namespace bridge
{

void RegisterUserNumber(Registry& reg, const char* parentGroup);
void RegisterUserVector(Registry& reg, const char* parentGroup);
void RegisterUserMatrix(Registry& reg, const char* parentGroup);

void RegisterBoundaryNumber(Registry& reg, const char* parentGroup);


void RegisterElderUserFunctions(Registry& reg, const char* parentGroup);

} // end namepace bridge
} // end namespace ug

#endif /* __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__ */
