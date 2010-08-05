#ifndef __H__UG__MARTIN_ALGEBRA__ALGEBRA_MISC__
#define __H__UG__MARTIN_ALGEBRA__ALGEBRA_MISC__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{
#ifndef NDEBUG
//!
//! use this to force the creation of prsize_t routines or similar for use in gdb.
#define FORCE_CREATION volatile size_t ___never_happens___ = 0; if(___never_happens___)
#else
#define FORCE_CREATION if(0)
#endif

//! prevent unused variable-warnings
#define UNUSED_VARIABLE(var) ((void) var);
}


#endif // __H__UG__MARTIN_ALGEBRA__ALGEBRA_MISC__
