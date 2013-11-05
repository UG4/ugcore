// created by Martin Rupp
// martin.rupp@gcsc.uni-frankfurt.de
// y2012m04d27
#include "mpi.h"
#include "common/types.h"

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

// end group pcl
/// \}

}
