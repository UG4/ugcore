/*
 * endian.h
 *
 *  Created on: 31.08.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__COMMON__UTIL__ENDIAN__
#define __H__UG__COMMON__UTIL__ENDIAN__

namespace ug{

/// \addtogroup ugbase_common_util
/// \{

/**
 * Returns true if the system endianess is little endian, i.e. the least
 * significant byte is stored first. If returned false the system endianess is
 * big endian, i.e. most significant byte first. Please note, that this is a
 * runtime check.
 *
 * @return	true if little endian, false if big endian
 */
inline bool IsLittleEndian()
{
    short int number = 0x1;
    char *numPtr = (char*)&number;
    return (numPtr[0] == 1);
}

/**
 * Returns false if the system endianess is little endian, i.e. the least
 * significant byte is stored first. If returned true the system endianess is
 * big endian, i.e. most significant byte first. Please note, that this is a
 * runtime check.
 *
 * @return	false if little endian, true if big endian
 */
inline bool IsBigEndian()
{
	return !IsLittleEndian();
}

// end group ugbase_common_util
/// \}

} // end namespace ug

#endif /* __H__UG__COMMON__UTIL__ENDIAN__ */
