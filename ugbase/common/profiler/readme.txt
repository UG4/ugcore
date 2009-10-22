The profiler source code has been taken from the shiny-profiler project by Aidin Abedi.
See http://shinyprofiler.sourceforge.net/ and have a look at the "Shiny  C++ Profiler.htm" in this folder.
Some files have been slightly altered to ensure compatibility with ug and to remove compilation errors.
altered files are:
	ShinyMacros.h
	ShinyPrereqs.h
	ShinyNodePool.cpp

The profiler can be disabled by adding
#define SHINY_PROFILER = FALSE
before you include profiler.h

Sebastian	(s.b.reiter@googlemail.com)
