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

/*  Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
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

	static const int blockSize = 1;
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
	typedef  SparseMatrix<DenseMatrix<FixedArray2<double, TBlockSize, TBlockSize> > > matrix_type;
	typedef Vector<DenseVector<FixedArray1<double, TBlockSize> > > vector_type;
#endif

	static const int blockSize = TBlockSize;
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
	typedef  SparseMatrix<DenseMatrix<VariableArray2<double> > > matrix_type;
	typedef Vector<DenseVector<VariableArray1<double> > > vector_type;
#endif

	static const int blockSize = AlgebraType::VariableBlockSize;
	static AlgebraType get_type()
	{
		return AlgebraType(AlgebraType::CPU, AlgebraType::VariableBlockSize);
	}
};

}
#endif /* __H__UG__LIB_ALGEBRA__CPU_ALGEBRA_TYPES__ */
