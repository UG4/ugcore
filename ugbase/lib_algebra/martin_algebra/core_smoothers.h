/*
 *  core_smoothers.h
 *
 *  Created by Martin Rupp on 21.07.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__MARTIN_ALGEBRA__CORE_SMOOTHERS__
#define __H__UG__MARTIN_ALGEBRA__CORE_SMOOTHERS__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{
template<typename Matrix_type, typename Vector_type>
bool gs_step_LL(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// gs LL has preconditioning matrix
	// (D-L)^{-1}

	typename Vector_type::entry_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(typename Matrix_type::cLowerLeftIterator it = A.beginLowerLeftRow(i); !it.isEnd(); ++it)
			s -= (*it).dValue * x[(*it).iIndex];
		x[i] = s/A(i,i);
	}

	return true;
}

template<typename Matrix_type, typename Vector_type>
bool gs_step_UR(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// gs UR has preconditioning matrix
	// (D-U)^{-1}
	typename Vector_type::entry_type s;
	if(x.size() == 0) return true;
	for(size_t i = x.size()-1; ; --i )
	{
		s = b[i];
		for(typename Matrix_type::cUpperRightIterator it = A.beginUpperRightRow(i); !it.isEnd(); ++it)
		//for(typename Matrix_type::cRowIterator it = A.beginRow(i); !it.isEnd(); ++it)
		{
		//	if((*it).iIndex > i)
				s -= (*it).dValue * x[(*it).iIndex];
		}
		x[i] = s/A.get_diag(i);
		if(i==0) break;
	}

	return true;
}

template<typename Matrix_type, typename Vector_type>
bool sgs_step(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// sgs has preconditioning matrix
	// (D-U)^{-1} D (D-L)^{-1}

	// x1 = (D-L)^{-1} b
	gs_step_LL(A, x, b);

	// x2 = D x1
	for(size_t i = 0; i<x.size(); i++)
		x[i] *= A.get_diag(i);

	// x3 = (D-U)^{-1} x2
	gs_step_UR(A, x, x);

	return true;
}

template<typename Matrix_type, typename Vector_type>
bool diag_step(const Matrix_type& A, Vector_type& x, const Vector_type& b, number damp)
{
	//exit(3);
	UG_ASSERT(x.size() == b.size() && x.size() == A.num_rows(), x << ", " << b << " and " << A << " need to have same size.");

	for(size_t j=0; j < x.size(); j++)
		x[j] = b[j] / A.get_diag(j);

	return true;
}
}
#endif // __H__UG__MARTIN_ALGEBRA__CORE_SMOOTHERS__
