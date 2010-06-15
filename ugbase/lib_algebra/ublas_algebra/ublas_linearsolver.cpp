
#include "ublas_linearsolver.h"

namespace ug{

using namespace std;
/// jacobi step
bool diag_step(const UblasMatrix& A, UblasVector& c, UblasVector& d, number damp)
{
	typedef ublas::vector<double> ScalarVector;
	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;

	const ScalarMatrix& Amat = const_cast<ScalarMatrix&>(*A.getStorage());
	ScalarVector& cVec = *c.getStorage();
	ScalarVector& dVec = *d.getStorage();

	// invert approx(A):  c = diag(A)^{-1} * d
	mv_dsolve(Amat, cVec, dVec);

/* THIS IS NOW MOVED TO JACOBI SOLVER, only STEP left
	// damp correction
	c *= damp;

	// update defect
	// dVec -= Amat*cVec
	mvsub(Amat, cVec, dVec);
*/
	return true;
}


UblasJacobi::UblasJacobi(int maxIter, number tol, bool verbose)
{
	m_maxIter = maxIter;
	m_tol = tol;
	m_verbose = verbose;
}

bool UblasJacobi::solve(UblasMatrix& A, UblasVector& x, UblasVector& b)
{
	typedef ublas::vector<double> ScalarVector;
	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;

	typedef MatrixLinOp<ScalarMatrix, ScalarVector, ScalarVector> ScalarLinop;
	typedef PJacobi<ScalarMatrix, ScalarVector, ScalarVector> ScalarPrecond;
	typedef StdConvCheck ScalarConvCheck;
	typedef LinOpInverse<ScalarLinop, ScalarPrecond, ScalarConvCheck> ScalarSolver;

	ScalarLinop linop_s(*A.getStorage());
	ScalarPrecond pjac_s(*A.getStorage());
	ScalarConvCheck conv_s(m_maxIter, m_tol, m_verbose);
	ScalarSolver solver_s(linop_s, pjac_s, conv_s);

	solver_s.apply(*x.getStorage(),*b.getStorage());

	return true;
}

bool UblasJacobi::step(UblasMatrix& A, UblasVector& c, UblasVector& d, number damp)
{
	typedef ublas::vector<double> ScalarVector;
	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
	typedef PJacobi<ScalarMatrix, ScalarVector, ScalarVector> ScalarPrecond;

	mv_dsolve(*A.getStorage(), *c.getStorage(), *d.getStorage());

	*c.getStorage() *= damp;

	mvsub(*A.getStorage(), *c.getStorage(), *d.getStorage());

	return true;
}

UblasJacobi::~UblasJacobi()
{

}


}
