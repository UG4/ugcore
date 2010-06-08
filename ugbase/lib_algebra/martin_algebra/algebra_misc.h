#pragma once

////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef NDEBUG
static int ___never_happens___ = 0;
//!
//! use this to force the creation of print routines or similar for use in gdb.
#define FORCE_CREATION if(___never_happens___)
#else
#define FORCE_CREATION if(0)
#endif
