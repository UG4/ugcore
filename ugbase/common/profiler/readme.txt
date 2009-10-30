The profiler source code has been taken from the shiny-profiler project by Aidin Abedi.
See http://shinyprofiler.sourceforge.net/ and have a look at the "Shiny  C++ Profiler.htm" in this folder.
Some files have been slightly altered to ensure compatibility with ug and to remove compilation errors.
altered files are:
	ShinyMacros.h
	ShinyPrereqs.h
	ShinyNodePool.cpp
	ShinyTools.h

ShinyTools.h defines a hash-function. The original implementation which simply casted the pointer to and uint32_t did not compile on 64-bit machines. I added a small workaround. Sadly this will most likely introduce errors in the hashing on 64 bit platforms.
If you are using a 64 bit platform and experience strange profiler-behaviour you should definitively think about improving hashing in shiny.

The profiler can be disabled by adding
#define SHINY_PROFILER FALSE
before you include profiler.h

Sebastian	(s.b.reiter@googlemail.com)
