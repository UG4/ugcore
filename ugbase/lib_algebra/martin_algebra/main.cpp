#if 1
/*
 *  main.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#include "algebra.h"


void discretization(myMatrix &A, myVector &x, myVector &b, int NX, int NY, int NZ, bool full_star, int problemtype);
double getDiscErr(myVector &vx, int NX, int NY, int NZ, int problemtype);

extern std::vector<postype> positions;


/*#define VECLIB
#include <flens/flens.h>

using namespace flens;
using namespace std;
*/
int flexamg_dimensions;

// 16704ms
////////// main
////////////////////////////////////////////////////////////////



//! simple matrix multiplication test
void test3(myMatrix &A, myVector &x, myVector &b)
{
	stopwatch SW;
	SW.start();
	b.print();
	x += 0.3*x + 0.5*b - A*b;
	x.print();

	SW.printTimeDiff();				
}

//! simple sgs/LinearSolver test
void sgstest(myMatrix &A, myVector &x, myVector &b)
{
	stopwatch SW;
	// CG / sgs test
	SW.start();
	sgs<myMatrix, myVector> prec;
	LinearSolver(x, A, b, prec, 10, 1e-8);
	cout << "took " << SW << endl;
}

//! AMG test: create one AMG level (do not solve)
void onelevel_amg(myMatrix &A, int aggressive_coarsening)
{
	amg<myMatrix, myVector> AMG;
	if(aggressive_coarsening)
		AMG.setAggressiveCoarsening_A_2();
		
	AMG.onlyOneLevel(A);
	if(A.getLength() < 10003)
	{
		writeToPosFile("/tmp/grid");
		AMG.writeMatrices("/Users/mrupp/matrices/SparseMatrix");
	}
}

//! AMG test: AMG hierachy analysis
int amg_test(myMatrix &A, myVector &x, myVector &b, int aggressive_coarsening)
{	
	amg<myMatrix, myVector> AMG;
	if(aggressive_coarsening)
		AMG.setAggressiveCoarsening_A_2();
	
	AMG.amgTest(A,x,b);
	return 0;
}

//! AMG test: create full AMG hierachy, solve
int solve_flexamg(myMatrix &A, myVector &x, myVector &b, int aggressive_coarsening, bool bCG)
{	
	amg<myMatrix, myVector> AMG;
	if(aggressive_coarsening)
		AMG.setAggressiveCoarsening_A_2();
	
	if(bCG)
		CG(x, A, b, AMG, 10);	// with preconditioner
	else
		LinearSolver(x, A, b, AMG, 40, 1e-12); // without preconditioner

	if(A.getLength() < 10003)
	{
		writeToPosFile("/tmp/grid");
		//AMG.writeMatrices("/Users/mrupp/matrices/SparseMatrix");
		x.printtofile(OUTPUT_DIR "SparseMatrix/vx.dat"); 
	}	

    return 0;
}

#ifndef USE_COMMAND_LINE
int main2(int argc, char **argv) ;
int main(int argc, char **argv) 
{
	int argc2 = 8;
	char *argv2[] = { "flexamg",
		"33", // NX
		"33", // NY
		"9", // NZ
		"2", // dim
		"0", // fullstar
		"0", // problemtype
		"0", // test
		NULL};
	return main2(argc2, argv2);
}

/*
template<typename T>
int inverse(T &A)
{
	DenseVector<Array<int> > P(A.numRows());
	int info = trf(A, P);
	if(info) return info;
	info = tri(A, P);
	return info;
}*/

int main2(int argc, char **argv) 
#else
int main(int argc, char **argv) 
#endif
{
	cout << "                           FLEXAMG" << endl;
	cout << "==================================================================" << endl;
	
	if(argc < 8)
	{
		cout << "Not enough parameters!" << endl;
		cout << "[1] : NX" << endl;
		cout << "[2] : NY" << endl;
		cout << "[3] : NZ" << endl;
		cout << "[4] : 2/3 (dimensions)" << endl;
		cout << "[4] : 0/1 (full_star)" << endl;		
		cout << "[5] : 0/1 (problemtype)" << endl;
		cout << "[6] : 0-4 (testtype: 0 = full AMG test, 1 = sgs test, 2 = multiplication test, 3 = one level AMG init test, 4 = amg analysis)" << endl;		
		return 1;
	}
	
	int NX = atoi(argv[1]); cout << "NX = " << NX;
	int NY = atoi(argv[2]); cout << ", NY = " << NY;
	int NZ = atoi(argv[3]); cout << ", NZ = " << NZ << endl;
	::flexamg_dimensions = atoi(argv[4]);
	
	bool full_star = atoi(argv[5]) != 0; if(full_star) cout << "full star, ";
	int problemtype = atoi(argv[6]); cout << "problem " << problemtype; 
	int test = atoi(argv[7]); cout << ", test " << test << endl << endl;
	
	
	myVector b; b.name = "b";
	myVector vx; vx.name = "x";	
	myMatrix A; A.name = "A";
	
	stopwatch SW;
	SW.start();
	

	// run discretization (see disc.cpp)
	discretization(A, vx, b, NX, NY, NZ, full_star, problemtype);
	
	
	/*myVector a1; a1.create(vx.getLength());
	myVector a2; a2.create(vx.getLength());*/
	
	int n = A.getLength();
	cout << "n = " << n << ", " << UNKNOWN_NR << " unknowns ";
#if UNKNOWN_NR > 1
	cout << "(" << myStorageType::getType() << " Storage)";
#endif
	cout << endl;
	cout << "==================================================================" << endl << endl;
	
	/*
	vector<int> indices;
	A.getNeighborhood(403, 3, indices);
	for(int i; i<indices.size(); i++)
		cout << indices[i] << endl;
	submatrix<myBlockMat> SM(&indices[0], indices.size());
	
	A.get(SM);
	//cout << SM << endl;
	SM.print();
	//////
	__CLPK_integer rows = SM.getRows();
	__CLPK_integer cols = SM.getCols();
	
	variableArray2<double> densemat;
	densemat.setSize(rows, cols);
	for(int r=0; r < rows; r++)
		for(int c=0; c < cols; c++)
			densemat[c + r*cols] = SM(r, c);

	__CLPK_integer interchange[20];	
	__CLPK_integer info = 0;
	
	dgetrf_(&rows, &cols, &SM(0,0), &rows, &interchange[0], &info);
	ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	
	double worksize; __CLPK_integer iWorksize = -1; 
	dgetri_(&rows, &SM(0,0), &rows, &interchange[0], &worksize, &iWorksize, &info);
	ASSERT1(info == 0);
	iWorksize = worksize;
	double *work = new double[iWorksize];
	
	dgetri_(&rows, &SM(0,0), &rows, &interchange[0], work, &iWorksize, &info);
	ASSERT1(info == 0);
	
	int ind = SM.getLocalColIndex(403);
	double alpha = SM(ind, ind);
	for(int i=0; i<SM.getRows(); i++)
		for(int j=0; j<SM.getCols(); j++)
			SM(i,j)/=alpha;
	
	//delete[] work;
	cout << endl << endl;
	SM.print();
	
	for(int i=0; i<SM.getRows(); i++)
		cout << endl << SM.getRowIndex(i) << ": " <<  SM(ind, i);
	cout << endl;
	
	//////
	
	
	A.get(SM);
	 

	//return 0;
	/*GeMatrix<FullStorage<double, ColMajor> > fSM( SM.getRows(), SM.getCols(), 0, 0);
	for(int i=0; i<SM.getRows(); i++)
		for(int j=0; j<SM.getCols(); j++)
		{
			fSM(i,j) = SM(i,j);
		}
	cout << "_____________________________________________________________" << endl << endl;
	
	cout << fSM;
	GeMatrix<FullStorage<double, ColMajor> > fSM2(SM.getRows(), SM.getCols(), 0, 0); fSM2 = fSM;
	inverse(fSM);
	
	cout << fSM;
	
	return 0;
	
	//A.print();
	/*b.print();
	vx.print();*/
	
	if(test == 0)
	{
		cout << "full AMG solving test ..." << endl << endl;
		solve_flexamg(A, vx, b, !full_star, false);	
		cout << endl << endl << "disc. err: " << getDiscErr(vx, NX, NY, NZ, problemtype) << endl;	
	}
	else if(test == 1)
	{
		cout << "sgs test ..." << endl << endl;
		sgstest(A, vx, b);
	}
	else if(test == 2)
	{
		cout << "multiplication test ..." << endl << endl;
		test3(A, vx, b);
	}
	else if(test == 3)
	{
		cout << "one level AMG init test ..." << endl << endl;
		onelevel_amg(A, !full_star);
	}
	else if(test == 4)
	{
		cout << "analysing AMG ..." << endl << endl;
		amg_test(A, vx, b, !full_star);
	}
	else
		cout << "unsupported test nr. " << test << endl;
	
	//vx.print();
	return 0;
}



#endif