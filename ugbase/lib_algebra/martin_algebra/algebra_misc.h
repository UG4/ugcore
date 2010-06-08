#pragma once

////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef NDEBUG

//!
//! use this to force the creation of print routines or similar for use in gdb.
#define FORCE_CREATION volatile int ___never_happens___ = 0; if(___never_happens___)
#else
#define FORCE_CREATION if(0)
#endif
