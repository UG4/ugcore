/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__
#define __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__

#include "algebra_type.h"

// vector and sparse_matrix
#include "cpu_algebra/vector.h"
#include "cpu_algebra/sparsematrix.h"

#ifdef UG_GPU
//#include "gpu_algebra/gpuvector.h"
//#include "gpu_algebra/gpusparsematrix.h"
#endif

// parallel support
#ifdef UG_PARALLEL
	#include "parallelization/parallel_vector.h"
	#include "parallelization/parallel_matrix.h"
#endif

#include "small_algebra/small_algebra.h"
namespace ug{

////////////////////////////////////////////////////////////////////////////////
//   CPU Algebra
////////////////////////////////////////////////////////////////////////////////

/**  This module defines different algebra types.
 *
 *  An Algebra should export the following typedef:
 *  - matrix_type (for representation of sparse matrices)
 *  - vector_type (for representation of DoF vectors)
 *  - flag_unit_type (integer type for the elementary units of the DoF flag storage)
 *
 *  Furthermore, Algebra should define the following constants:
 *  - blockSize (size of the elementary blocks, same as in AlgebraType)
 *  - flagBlocksPerUnit (number of flag blocks stored in one elementary unit)
 *
 *  The type 'vector_type' describes the storage of the values of the DoFs,
 *  so that 'vector_type' is a type of a "grid function without the grid".
 *  This type declares 'vector_type::value_type', the type of the values
 *  of the considered grid functions. Note that 'value_type' may be a scalar
 *  'double', or a small array (vector) of 'double's, or any other data
 *  structure. In any case, it is refered to as a 'DoF block'. 'Vector_type'
 *  can be considered as an array of these DoF blocks, and the indices
 *  identify the entire block (not a particular element of it). These
 *  indices are refered to as 'algebra indices'.
 *
 *  The 'matrix_type' describes representations of sparse matrices (or
 *  linear operators) on these 'vectors' (of type 'vector_type'). Note that
 *  the entries of these matrices are generally small matrices, too, and
 *  have types 'matrix_type::value_type'. Only in the case of scalar
 *  'vector_type::value_type', the type 'matrix_type::value_type' is scalar,
 *  too.
 *
 *  Operations with 'vector_type::value_type' and 'matrix_type::value_type'
 *  belong to the s.c. "small algebra", see lib_algebra/small_algebra. Using
 *  those template functions, the implementation can be made independent
 *  on the type.
 *
 *  Flags are assigned to DoF blocks. All flags assigned to the same DoF
 *  block are stored in one or several flag units having the same index
 *  as the DoF block itself. A mapping assigning the flag units to DoF block
 *  indices is called a flag set. (E.i., one can consider a flag set as
 *  an array of the flag units indexed in the same way as the DoF blocks
 *  themseves.)
 *
 *  Inside the flag units, the flags are organized in blocks. The flag blocks
 *  have the same size as the DoF blocks. Every flag in a flag block marks
 *  the corresponding DoF in the DoF block and has for it the same meeting
 *  as the other flags in the same block. Every entire flag block is stored
 *  completely in one flag unit. How the blocks are placed in different units
 *  (i.e. the allocation of the flag blocks in the flag sets), is regulated
 *  by the flag manager.
 */

////////////////////////////////////////////////////////////////////////////////
//   CPU Algebra (Block 1x1 Algebra)
////////////////////////////////////////////////////////////////////////////////

/**
 * \defgroup cpu_algebra CPU Algebra
 * \ingroup lib_algebra
 * \{
 */

struct CPUAlgebra
{
#ifdef UG_PARALLEL
		typedef ParallelMatrix<SparseMatrix<double> > matrix_type;
		typedef ParallelVector<Vector<double> > vector_type;
#else
		typedef SparseMatrix<double> matrix_type;
		typedef Vector<double> vector_type;
#endif

	static const int blockSize = 1; ///< i.e., storing scalar DoFs
	
	typedef typename vector_type::flag_unit_type flag_unit_type; ///< elementary type to store DoF flags
	static const int flagBlocksPerUnit = sizeof(flag_unit_type) * 8; ///< number of flag blocks per unit
	
	/// returns the algebra type
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, 1);
	}
};

// end group cpu_algebra
/// \}

/**
 * \defgroup crs_algebra CRS Algebra
 * \ingroup lib_algebra
 * \{
 */
/*#ifdef UG_GPU
struct GPUAlgebra
{
#ifdef UG_PARALLEL
		typedef ParallelMatrix<GPUSparseMatrix<double> > matrix_type;
		typedef ParallelVector<GPUVector<double> > vector_type;
#else
		typedef GPUSparseMatrix<double> matrix_type;
		typedef GPUVector<double> vector_type;
#endif

	static const int blockSize = 1;
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::GPU, 1);
	}
};
#endif
// end group crs_algebra
/// \}
*/
////////////////////////////////////////////////////////////////////////////////
//   CPU Fixed Block Algebra
////////////////////////////////////////////////////////////////////////////////
/// \addtogroup cpu_algebra
template<int TBlockSize>
struct CPUBlockAlgebra
{
#ifdef UG_PARALLEL
	typedef ParallelMatrix<SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > > matrix_type;
	typedef ParallelVector<Vector<DenseVector<FixedArray1<double, TBlockSize> > > > vector_type;
#else
	typedef SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > matrix_type;
	typedef Vector<DenseVector<FixedArray1<double, TBlockSize> > > vector_type;
#endif

	static const int blockSize = TBlockSize; ///< the algebra is operating with small blocks
	
	typedef typename vector_type::flag_unit_type flag_unit_type; ///< elementary type to store DoF flags
	static const int flagBlocksPerUnit = (sizeof(flag_unit_type) * 8) / blockSize; ///< number of flag blocks per unit
	
	/// returns the algebra type
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, TBlockSize);
	}
};

////////////////////////////////////////////////////////////////////////////////
//   CPU Variable Block Algebra
////////////////////////////////////////////////////////////////////////////////

/// \addtogroup cpu_algebra
struct CPUVariableBlockAlgebra
{
#ifdef UG_PARALLEL
	typedef ParallelMatrix<SparseMatrix<DenseMatrix<VariableArray2<double> > > > matrix_type;
	typedef ParallelVector<Vector<DenseVector<VariableArray1<double> > > > vector_type;
#else
	typedef SparseMatrix<DenseMatrix<VariableArray2<double> > > matrix_type;
	typedef Vector<DenseVector<VariableArray1<double> > > vector_type;
#endif

	static const int blockSize = AlgebraType::VariableBlockSize; ///< the algebra is operating with blocks of variable type
	
	typedef typename vector_type::flag_unit_type flag_unit_type; ///< elementary type to store DoF flags
	static const int flagBlocksPerUnit = 1; ///< only one flag block per unit as we cannot compute the offset
	/*ToDo: The size of the block is then restricted by sizeof(flag_unit_type)*8 */
	
	/// returns the algebra type
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, AlgebraType::VariableBlockSize);
	}
};

}
#endif /* __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__ */
