/* 
 * File:   const_user_data.h
 * Author: miho
 *
 * Created on 15. April 2011, 11:11
 */

#ifndef USER_DATA_H
#define	USER_DATA_H

namespace ug {
namespace vrl {

void RegisterUserData(ug::bridge::Registry& reg, const char* parentGroup);
//void RegisterBoundaryNumber(Registry& reg, const char* parentGroup);

} // vrl::
} // ug::

#endif	/* USER_DATA_H */

