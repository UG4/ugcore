// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// Somewhen back in 2005-2007

#include "common/math/ug_math.h"
#include "sparse_matrix.h"
#include "nvec_math_templates.h"

#ifndef __UG__LIB_GRID_LEAST_SQUARES_SOLVER__
#define __UG__LIB_GRID_LEAST_SQUARES_SOLVER__

namespace ug{
namespace libgrid_simplealg{

///	This method finds v so that the error |Av-d| for a linear system Av = d is minimized.
/**	The method is implemented by a sort of cg-method.*/
template <class Number>
int LeastSquaresSolver(libMesh::SparseMatrix<Number>& A, Number* v, Number* d,
						int MaxIter, int NumRows, int NumCols, Number rhs = SMALL)
{
	Number* r = new Number[NumCols];
	Number* p = new Number[NumCols];
	Number* q = new Number[NumCols];
	Number* tmp = new Number[NumRows];

	A.vmultT(p, d);
	A.vmult(tmp, v);
	A.vmultT(r, tmp);

	NVecSub(r, p, r, (Number)1., NumCols);

	Number fpo, fpn;
	int i;
	for(i = 0; i < MaxIter; ++i)
	{
		fpn = NVecDot(r, r, NumCols);

		if(fpn < rhs)
			break;

		if(i == 0)
			NVecCopy(p, r, NumCols);
		else
		{
			Number beta = fpn / fpo;
			NVecAdd(p, r, p, beta, NumCols);
		}
		A.vmult(tmp, p);
		A.vmultT(q, tmp);
		Number alpha = fpn / NVecDot(p, q, NumCols);
		NVecAdd(v, v, p, alpha, NumCols);
		NVecSub(r, r, q, alpha, NumCols);

		fpo = fpn;
	}
	delete[] r;
	delete[] p;
	delete[] q;
	delete[] tmp;

	return i;
}

}//	end of namespace
}//	end of namespace
#endif
