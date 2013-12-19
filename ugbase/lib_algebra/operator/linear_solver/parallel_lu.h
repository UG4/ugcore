/*
 * parallel_lu.h
 *
 *  Created on: 03.12.2013
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_PARALLEL_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_PARALLEL_LU_OPERATOR__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioner.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/collect_matrix.h"
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
#endif

namespace ug{


template <typename TBase, typename TAlgebra>
class AgglomeratingBase : public TBase
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	// 	Destructor
		virtual ~AgglomeratingBase() {};



	std::vector<size_t> originalToA;
	std::vector<size_t> AtoOriginal;
	size_t borderStart=k;

	void init(matrix_type &oldMat)
	{
		matrix_type mat;

		IndexLayout totalMaster, totalSlave;
		std::vector<IndexLayout> vMasterLayouts, vSlaveLayouts;
		// braucht man nciht
		GenerateOverlapClass<matrix_type> c(oldMat, mat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts);

		c.m_overlapDepthMaster = 1;
		c.m_overlapDepthSlave = 1;
		c.m_masterDirichletLast = false;
		c.m_slaveDirichletLast = false;
		bool b = c.calculate();
		overlapSize = c.m_overlapSize;
		return b;

		/////

		size_t N = mat.num_rows();
		std::vector<size_t> originalToA(N), originalToInterface(N);
		std::vector<size_t> AtoOriginal, InterfaceToOriginal;

		std::vector<bool> border(N, false);

		ConstSmartPtr<AlgebraLayouts> layouts = mat.layouts();

		MarkAllFromLayout(mark, totalMasterLayout);
		MarkAllFromLayout(mark, totalSlaveLayout);

		size_t k=0;
		for(size_t i=0; i<N; i++)
		{
			if(border[i]==false)
			{
				originalToA[i] = k;
				AtoOriginal.push_back(i);
				k++;
			}
		}
		size_t Asize = k;
		k = 0;
		for(size_t i=0; i<N; i++)
		{
			if(border[i]==true)
			{
				originalToInterface[i] = k;
				InterfaceToOriginal.push_back(i);
				k++;
			}
		}
		size_t interfaceSize = k;


		////////////////////////

		matrix_type A_ii, A_ib, A_bi, A_bb;

		CreateMatrix(A_ii, A, originalToInterface, interfaceSize, originalToInterface, interfaceSize);
		CreateMatrix(A_ai, A, originalToA, Asize, originalToInterface, interfaceSize);
		CreateMatrix(A_ia, A, originalToInterface, interfaceSize, originalToA, Asize);
		CreateMatrix(A_aa, A, originalToA, Asize, originalToA, Asize);

		lu.init(A_ii);

		S = A_aa;
		for(size_t i=0; i<interfaceSize; i++)
		{
			vector_type v;
			A_ia.get_column_vector(i, v);
			//A_ai * lu^{-1} * A_ia[i, .]
			lu.apply(x, v);
//			S_1.set_col(i, x);


			// A (v1 v2 v3) = (Av1 Av2 Av3)

			// S = A_aa + A_ai * S_1;
			v = A_ai * x;
			S.add_col(i, v)
		}

//		S = A_aa + A_ai * S_1;

		// grid
		////////////////////
		/// aaaaaibbbbb
		/// aaaaaibbbbb
		/// aaaaaibbbbb
		/// aaaaaibbbbb

		// P1
		/// aaaaai
		/// aaaaai
		/// aaaaai
		/// aaaaai

		// P2
		///      ibbbbb
		///      ibbbbb
		///      ibbbbb
		///      ibbbbb

		// A_11 A_12 A_13
		// A_21 A_22 A_23
		// A_31 A_32 A_33

		// A_11 = _aa : a -> a
		// A_31 = _ia : a -> i
		// A_13 = _ai : i -> a

		// A_22 = _bb
		// A_21 = _ib
		// A_12 = _bi

		// A_33 = _ii


		// A_11		0		A_13
		// 0		A_22	A_23
		// A_31		A_32	A_33

		// S = A_{33} + A_{13} A_11^{-1} A_{31} + A_{23} A_22^{-1} A_{21}

	}

};




} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
