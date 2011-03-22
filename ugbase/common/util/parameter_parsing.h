// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)

#ifndef __H__UG__parameter_parsing__
#define __H__UG__parameter_parsing__

namespace ug
{
////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns its position in argv.
 *	If the parameter is not contained in argv, -1 is returned.
 */
int GetParamIndex(const char* param, int argc, char* argv[]);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns true if it is found.
 */
bool FindParam(const char* param, int argc, char* argv[]);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter, and converts the
 *	associated value to an integer. Returns true if the parameter was
 *	found, false if not.
 */
bool ParamToInt(int& iOut, const char* param, int argc, char* argv[]);

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter, and returns the
 *	associated string (the argv directly following param).
 *	associated Returns true if the parameter was found, false if not.
 *
 *	Please note that spaces may not be contained in the associated string.
 */
bool ParamToString(char** strOut, const char* param, int argc, char* argv[]);

}//	end of namespace

#endif
