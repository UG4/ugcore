#ifndef __H__UG__CPU_ALGEBRA__ALGEBRA_MISC__
#define __H__UG__CPU_ALGEBRA__ALGEBRA_MISC__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{

/// \addtogroup lib_algebra
/// \{

#ifndef NDEBUG
//!
//! use this to force the creation of prsize_t routines or similar for use in gdb.
#define FORCE_CREATION volatile size_t ___never_happens___ = 0; if(___never_happens___)
#else
#define FORCE_CREATION if(0)
#endif

//! prevent unused variable-warnings
#define UNUSED_VARIABLE(var) ((void) var);


//!
//! template struct for sorting some keys after values
//! for example, sorting a vector of ints and know original pos
template<typename T>
struct sortStruct
{
	size_t index; // for example "original" position.
	T sortValue;

	bool operator < (const sortStruct<T> &other) const
	{
		return sortValue < other.sortValue;
	}
};

// end group lib_algebra
/// \}

}


#endif // __H__UG__CPU_ALGEBRA__ALGEBRA_MISC__
