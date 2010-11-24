//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d21

//	this is just a wrapper-include for the shiny-profiler by Aidin Abedi

//	To enable or disable the profiler the following define can be set.
//	It is preferable to do this via cmake options or similar.
//#define UG_PROFILER
#ifdef UG_PROFILER
	#define SHINY_PROFILER TRUE
#else
	#define SHINY_PROFILER FALSE
#endif

#include "src/Shiny.h"
