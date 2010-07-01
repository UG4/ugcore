#ifndef __H__UG__MARTIN_ALGEBRA__ALGEBRA_MISC__
#define __H__UG__MARTIN_ALGEBRA__ALGEBRA_MISC__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug{
#ifndef NDEBUG
//!
//! use this to force the creation of prsize_t routines or similar for use in gdb.
#define FORCE_CREATION volatile size_t ___never_happens___ = 0; if(___never_happens___)
#else
#define FORCE_CREATION if(0)
#endif

template<typename Matrix_type, typename Vector_type>
bool gs_step_LL(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// gs LL has preconditioning matrix
	// (D-L)^{-1}

	typename Vector_type::entry_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(typename Matrix_type::cLowerLeftIterator it = A.beginRow(i); !it.isEnd(); ++it)
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
	for(size_t i = x.size()-1; i >= 0; i--)
	{
		s = b[i];
		for(typename Matrix_type::cUpperRightIterator it = A.beginRow(i); !it.isEnd(); ++it)
			s -= (*it).dValue * x[(*it).iIndex];
		x[i] = s/A(i,i);
	}

	return true;
}

template<typename Matrix_type, typename Vector_type>
bool sgs_step(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// sgs has preconditioning matrix
	// (D-U)^{-1} D (D-L)^{-1}

	// x1 = (D-L)^{-1} b
	typename Vector_type::entry_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(typename Matrix_type::cLowerLeftIterator it = A.beginRow(i); !it.isEnd(); ++it)
			s -= (*it).dValue * x[(*it).iIndex];
		x[i] = s/A(i,i);
	}

	// x2 = D x1
	for(size_t i = 0; i<x.size(); i++)
		x[i] *= A(i,i);

	// x3 = (D-U)^{-1} x2
	for(size_t i = x.size()-1; i >= 0; i--)
	{
		s = x[i];
		for(typename Matrix_type::cUpperRightIterator it = A.beginRow(i); !it.isEnd(); ++it)
			s -= (*it).dValue * x[(*it).iIndex];
		x[i] = s/A(i,i);
	}

	return true;
}

template<typename Matrix_type, typename Vector_type>
bool diag_step(const Matrix_type& A, Vector_type& x, const Vector_type& b, number damp)
{
	//exit(3);
	UG_ASSERT(x.size() == b.size() && x.size() == A.row_size(), x << ", " << b << " and " << A << " need to have same size.");

	for(size_t j=0; j < x.size(); j++)
		x[j] = b[j] / A.getDiag(j);

	return true;
}

}

#endif
