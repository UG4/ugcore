
#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__

namespace ug
{
namespace bridge
{

void RegisterUserNumber(Registry& reg);
void RegisterUserVector(Registry& reg);
void RegisterUserMatrix(Registry& reg);
void RegisterElderUserFunctions(Registry& reg);

} // end namepace bridge
} // end namespace ug

#endif /* __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__ */
