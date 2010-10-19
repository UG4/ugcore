/*
 * parallel_storage_type.h
 *
 *  Created on: 3.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_STORAGE_TYPE__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_STORAGE_TYPE__

namespace ug
{


/**
 * The storage type of a vector is used in parallel applications.
 * We assume that the dofs are distributed to the processes in the way that
 * each dof is master on exactly one process and can be a slave (i.e. a local
 * copy) on several other processes. Given the real values of the dofs the different
 * storage type are defined as follows:
 *  - PST_UNDEFINED: no information given
 *  - PST_CONSISTENT: The real value is saved in the master and every slave
 *  - PST_ADDITIVE: The sum over the values in the master and all slaves gives the exact value
 *  - PST_UNIQUE: Same as PST_ADDITIV, but value is zero in all slaves (i.e. master has exact value)
 *
 *  Note, that a Vector can have more than one type. E.g. every unique Vector is additive. Moreover,
 *  the Vector being zero everywhere is consistent, additive and unique at the same time. Therefore,
 *  the information is given bitwise.
 *  To check, whether a function is in a type one may ask:
 *  	 -- v.has_storage_type(PST_CONSISTENT), etc.
 *  To change into another type
 *       -- v.change_storage_type(PST_ADDITIVE)
 */
enum ParallelStorageType {
	PST_UNDEFINED = 0,
	PST_CONSISTENT = 1 << 0,
	PST_ADDITIVE = 1 << 1,
	PST_UNIQUE = 1 << 2
};

inline ParallelStorageType operator & (const ParallelStorageType &a, const ParallelStorageType &b)
{
	return (ParallelStorageType) ((int)a&(int)b);
}

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_STORAGE_TYPE__ */
