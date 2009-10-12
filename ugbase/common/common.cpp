//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#include <fstream>
#include "common.h"

#ifdef LG_DEF__ENABLE_LOGGING
std::ofstream& lib_grid_logger()
{
	static std::ofstream fLog("libGrid_LOG.txt");
	return fLog;
}
#endif
