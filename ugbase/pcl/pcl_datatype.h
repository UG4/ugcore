/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <mpi.h>
#include "common/types.h"
#include "common/error.h"
#ifndef __H__PCL__pcl_datatype__
#define __H__PCL__pcl_datatype__

namespace pcl
{

/// \addtogroup pcl
/// \{

//	DataType
#define PCL_DT_NULL				MPI_DATATYPE_NULL
#define PCL_DT_BYTE 				MPI_BYTE
#define PCL_DT_PACKED 				MPI_PACKED
#define PCL_DT_CHAR 				MPI_CHAR
#define PCL_DT_SHORT 				MPI_SHORT
#define PCL_DT_INT 				MPI_INT
#define PCL_DT_LONG 				MPI_LONG
#define PCL_DT_UNSIGNED_LONG		MPI_UNSIGNED_LONG
#define PCL_DT_LONG_LONG_INT		MPI_LONG_LONG_INT
#define PCL_DT_UNSIGNED_LONG_LONG	MPI_UNSIGNED_LONG_LONG
#define PCL_DT_FLOAT 				MPI_FLOAT
#define PCL_DT_DOUBLE 				MPI_DOUBLE
#define PCL_DT_LONG_DOUBLE 		MPI_LONG_DOUBLE
#define PCL_DT_UNSIGNED_CHAR 		MPI_UNSIGNED_CHAR

typedef MPI_Datatype DataType;

class DataTypeDirectlySupported {};
class DataTypeIndirectlySupported {};

template<typename T>
class DataTypeTraits
{
public:
	static DataType get_data_type() {return T::PCL_DATATYPE_NOT_SUPPORTED(); }
	typedef DataTypeIndirectlySupported supported;
	enum { directlySupported = false };
};

template<>
class DataTypeTraits<unsigned long>
{
public:
	static DataType get_data_type() {return PCL_DT_UNSIGNED_LONG; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};

template<>
class DataTypeTraits<long>
{
public:
	static DataType get_data_type() {return PCL_DT_LONG; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};
template<>
class DataTypeTraits<int>
{
public:
	static DataType get_data_type() {return PCL_DT_INT; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};
template<>
class DataTypeTraits<float>
{
public:
	static DataType get_data_type() {return PCL_DT_FLOAT; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};
template<>
class DataTypeTraits<double>
{
public:
	static DataType get_data_type() {return PCL_DT_DOUBLE; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};

template<>
class DataTypeTraits<char>
{
public:
	static DataType get_data_type() {return PCL_DT_CHAR; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};

template<>
class DataTypeTraits<unsigned char>
{
public:
	static DataType get_data_type() {return PCL_DT_UNSIGNED_CHAR; }
	typedef DataTypeDirectlySupported supported;
	enum { directlySupported = true };
};

inline size_t GetSize(const DataType &t)
{
	if(t == PCL_DT_UNSIGNED_CHAR) return sizeof(unsigned char);
	else if(t == PCL_DT_CHAR) return sizeof(char);
	else if(t == PCL_DT_DOUBLE) return sizeof(double);
	else if(t == PCL_DT_FLOAT) return sizeof(float);
	else if(t == PCL_DT_INT) return sizeof(int);
	else if(t == PCL_DT_LONG) return sizeof(long);
	else if(t == PCL_DT_UNSIGNED_LONG) return sizeof(unsigned long);
	UG_THROW("Datatype not supported: " << t << " ???");
	return 1;
}

// end group pcl
/// \}

}

#endif
