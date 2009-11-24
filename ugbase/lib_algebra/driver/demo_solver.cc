#include "BoostBlock.hh"
#include "ExtAlgSolver.hh"


const int myBlockSize=2;
const int myGridSize = 3;
const int myNumDof = myGridSize*myGridSize;

typedef block_vector<double,myBlockSize> FixedSizeVector;
typedef ublas::vector<FixedSizeVector> BlockVector;


typedef block_matrix<double,myBlockSize,myBlockSize> FixedSizeMatrix;

//typedef ublas::generalized_vector_of_vector<FixedSizeMatrix, ublas::row_major,
//					    ublas::vector<ublas::coordinate_vector<FixedSizeMatrix> > > TempMatrix;
typedef ublas::coordinate_matrix<FixedSizeMatrix, ublas::row_major> TempMatrix;
typedef ublas::compressed_matrix<FixedSizeMatrix, ublas::row_major> BlockMatrix;
typedef ublas::compressed_matrix<FixedSizeMatrix, ublas::row_major> BlockMatrix;


//typedef matrix_slice<block_matrix> FixedSizeMatrixSlice;
//typedef ublas::compressed_matrix<FixedSizeMatrixSlice, ublas::row_major> BlockMatrixSlice;

// unknown-approach
typedef ublas::vector<double> ScalarVector;
typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
typedef block_vector<ScalarVector*,2> UnknownRefVector;
typedef block_matrix<ScalarMatrix*,2,2> UnknownRefMatrix;



template <typename M>
void AssembleMatrix(M & Mat, int N)
{
  typedef typename M::value_type S;
  const int bsize = myBlockSize;
 
  S Id =ublas::identity_matrix<double>(bsize);
  S Aii(Id); 
  Aii*=4.0;
  
  S Aij(Id); 
  Aij*=-1.0;
  
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      {
	// diagonal element
	int ind = i*N+j;
	
	if (ind-N >= 0)	
	  Mat.push_back(ind, ind-N, Aij);

	if (ind-1 >= 0)	
	  Mat.push_back(ind, ind-1, Aij);
	
	Mat.push_back(ind, ind, Aii);

	if (ind+1 < Mat.size2())
	  Mat.push_back(ind, ind+1, Aij);

	if (ind+N < Mat.size2())
	  Mat.push_back(ind, ind+N, Aij);

	// more offdiag elements
	// ...
      }
}


template <typename M>
void AssembleSparseMatrix(M & Mat, int N)
{
  typedef typename M::value_type S;
  S Aii=4.0;
  S Aij=-1.0;
  
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      {
	// diagonal element
	int ind = i*N+j;
	
	if (ind-N >= 0)	
	  Mat.push_back(ind, ind-N, Aij);

	if (ind-1 >= 0)	
	  Mat.push_back(ind, ind-1, Aij);
	
	Mat.push_back(ind, ind, Aii);

	if (ind+1 < Mat.size2())
	  Mat.push_back(ind, ind+1, Aij);

	if (ind+N < Mat.size2())
	  Mat.push_back(ind, ind+N, Aij);

	// more offdiag elements
	// ...
      }
}
template <typename T> 
void dummy(T t) {
  std::cout << "This non-int!" <<  std::endl;
};


void dummy(int i)
{
  std::cout << "This is int!" <<  std::endl; 
}





int main ()
{

  typedef MatrixLinOp<BlockMatrix, BlockVector, BlockVector> Linop;
  typedef PJacobi<BlockMatrix, BlockVector, BlockVector> Precond;
  typedef StdConvCheck ConvCheck;
  typedef LinOpInverse<Linop, Precond, ConvCheck> Solver;
  
  

  // A) initialize matrix
  TempMatrix Atmp(myNumDof, myNumDof);  
  //BlockMatrix Atmp(myNumDof, myNumDof);
  AssembleMatrix(Atmp, myGridSize);
  std::cout << Atmp << std::endl;
  
  


  BlockMatrix A(Atmp);  // fails for 
  //BlockMatrix A(myNumDof,myNumDof); 
  //A = Atmp;
  std::cout << A << std::endl;
 
  BlockVector sol(myNumDof), rhs(myNumDof);


  
  // diagonal solves
  mv_dsolve(A, sol, rhs); // two equivalent calls
  mv_dsolve_rec<2>(A, sol, rhs);
  mv_dsolve_rec<1>(A, sol, rhs);  // block level call (must fail)

  // lower part solve
  /*
GenLSolve(A, sol, rhs); // two equivalent calls
  GenLSolveRec<2>(A, sol, rhs);
  GenLSolveRec<1>(A, sol, rhs);  // block level call (must fail)
  */
  // B) initialize solver 
  
  Linop linop(A);
  Precond pjac(A);
  ConvCheck conv(10, 1.0e-3);

  Solver solver(linop, pjac, conv);
  solver.apply(sol, rhs);
  
  return 0;



  /*
  // defined in traits.hpp
  typedef boost::numeric::ublas::promote_traits<FixedSizeVector, FixedSizeVector>::promote_type ptype;
  ptype a = 4; // FixedSizeVector
  
  typedef boost::numeric::ublas::promote_traits<FixedSizeMatrix, FixedSizeMatrix>::promote_type pmtype;
  pmtype b = 4.0; // FixedSizeMatrix
    
  typedef FixedSizeMatrix::vector_temporary_type mtemp_type;
  mtemp_type c = 4;  // okay: compressed_vector<double>
  
  typedef BlockMatrix::vector_temporary_type  mtemp_type2;  
  mtemp_type2 d = 4;   // fails: compressed_vector<FixedSizeMatrix>, should be compressed_vector<compressed_vector<double>>
 
  typedef BlockMatrix::vector_temporary_type::value_type mtemp_type3;  
  mtemp_type3 e = 4; // FixedSizeMatrix
  
  typedef BlockMatrix::vector_temporary_type::value_type::vector_temporary_type mtemp_type4;  
  mtemp_type4 f = 4; //  okay: compressed_vector<double> a.k.a FixedSizeVector
  return 0;
  */



  /*
  // all crap!
  ScalarVector u1(myNumDof), u2(myNumDof), f1(myNumDof), f2(myNumDof);
  ScalarMatrix A11(myNumDof, myNumDof), A12(myNumDof, myNumDof), A21(myNumDof, myNumDof), A22(myNumDof, myNumDof);
  AssembleSparseMatrix(A11, 3);  
  AssembleSparseMatrix(A12, 3);  
  AssembleSparseMatrix(A21, 3);  
  AssembleSparseMatrix(A22, 3);  
  
  std::cout << A11 <<  std::endl;
  std::cout << A22 <<  std::endl;
  
  UnknownRefVector uVector, fVector;
  UnknownRefMatrix uMatrix;
  std::cout << uMatrix <<  std::endl;
  std::cout << uVector <<  std::endl;
  std::cout << fVector <<  std::endl;

  uMatrix.push_back(0,0, &A11); 
  uMatrix.push_back(0,1, &A12); 
  uMatrix.push_back(1,0, &A21); 
  uMatrix.push_back(1,1, &A22);
  
  uVector(0) = (&u1);  uVector(1) = (&u2);
  fVector(0) = (&f1);  fVector(1) = (&f2);

  std::cout << *uMatrix(0,0) <<  std::endl;
  std::cout << *uVector(0) <<  std::endl;
  std::cout << *fVector(0) <<  std::endl;
  */

  //  mv_prod(uMatrix, uVector, fVector);
}
