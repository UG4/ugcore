
#include "arnelinearsolver.h"

namespace ug{


/// jacobi step
bool diag_step(const ArneMatrix& A, ArneVector& c, ArneVector& d, number damp)
{
	typedef ublas::vector<double> ScalarVector;
	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;

	const ScalarMatrix& Amat = const_cast<ScalarMatrix&>(*A.getStorage());
	ScalarVector& cVec = *c.getStorage();
	ScalarVector& dVec = *d.getStorage();

	// invert approx(A):  c = diag(A)^{-1} * d
	mv_dsolve(Amat, cVec, dVec);

	// damp correction
	c *= damp;

	// update defect
	mvsub(Amat, cVec, dVec);

	return true;
}


ArneJacobi::ArneJacobi(int maxIter, number tol, bool verbose)
{
	_maxIter = maxIter;
	_tol = tol;
	_verbose = verbose;
}

bool ArneJacobi::solve(ArneMatrix& A, ArneVector& x, ArneVector& b)
{
	typedef ublas::vector<double> ScalarVector;
	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;

	typedef MatrixLinOp<ScalarMatrix, ScalarVector, ScalarVector> ScalarLinop;
	typedef PJacobi<ScalarMatrix, ScalarVector, ScalarVector> ScalarPrecond;
	typedef StdConvCheck ScalarConvCheck;
	typedef LinOpInverse<ScalarLinop, ScalarPrecond, ScalarConvCheck> ScalarSolver;

	ScalarLinop linop_s(*A.getStorage());
	ScalarPrecond pjac_s(*A.getStorage());
	ScalarConvCheck conv_s(_maxIter, _tol, _verbose);
	ScalarSolver solver_s(linop_s, pjac_s, conv_s);

	solver_s.apply(*x.getStorage(),*b.getStorage());

	return true;
}

bool ArneJacobi::step(ArneMatrix& A, ArneVector& c, ArneVector& d, number damp)
{
	typedef ublas::vector<double> ScalarVector;
	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
	typedef PJacobi<ScalarMatrix, ScalarVector, ScalarVector> ScalarPrecond;

	mv_dsolve(*A.getStorage(), *c.getStorage(), *d.getStorage());

	*c.getStorage() *= damp;

	mvsub(*A.getStorage(), *c.getStorage(), *d.getStorage());

	return true;
}

ArneJacobi::~ArneJacobi()
{

}


}
