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

using namespace ug;

void discretization(myMatrix &A, myVector &x, myVector &b, int NX, int NY, int NZ, bool full_stencil, int problemtype);
double getDiscErr(myVector &vx, int NX, int NY, int NZ, int problemtype);



namespace ug{
int flexamg_dimensions;
}


// 16704ms
////////// main
////////////////////////////////////////////////////////////////

//! simple matrix multiplication test
void test3(myMatrix &A, myVector &x, myVector &b)
{
	stopwatch SW;
	SW.start();
	x += 0.3*x + 0.5*b - A*b;

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
int main2(int argc, const char **argv) ;
int main(int argc, const char **argv) 
{
	int argc2 = 8;
	const char *argv2[] = { "flexamg",
		"513", // NX
		"513", // NY
		"3", // NZ
		"2", // dim
		"0", // full_stencil
		"0", // problemtype
		"4", // test
		NULL};
	return main2(argc2, argv2);
}

int main2(int argc, const char **argv) 
#else
int main(int argc, const char **argv) 
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
		cout << "[4] : 0/1 (full_stencil)" << endl;		
		cout << "[5] : 0/1 (problemtype)" << endl;
		cout << "[6] : 0-4 (testtype: 0 = full AMG test, 1 = sgs test, 2 = multiplication test, 3 = one level AMG init test, 4 = amg analysis)" << endl;		
		return 1;
	}
	
	int NX = atoi(argv[1]); cout << "NX = " << NX;
	int NY = atoi(argv[2]); cout << ", NY = " << NY;
	int NZ = atoi(argv[3]); cout << ", NZ = " << NZ << endl;
	::flexamg_dimensions = atoi(argv[4]);
	
	bool full_stencil = atoi(argv[5]) != 0; if(full_stencil) cout << "full stencil, ";
	int problemtype = atoi(argv[6]); cout << "problem " << problemtype; 
	int test = atoi(argv[7]); cout << ", test " << test << endl << endl;
	
	
	myVector b; b.name = "b";
	myVector vx; vx.name = "x";	
	myMatrix A; A.name = "A";
	
	stopwatch SW;
	SW.start();
	

	// run discretization (see disc.cpp)
	discretization(A, vx, b, NX, NY, NZ, full_stencil, problemtype);
	
	/*A.print();
	vx.print();
	b.print();*/
	
	
	int n = A.getLength();
	cout << "n = " << n << ", " << UNKNOWN_NR << " unknowns ";
#if UNKNOWN_NR > 1
	cout << "(" << myStorageType::getType() << " Storage)";
#endif
	cout << endl;
	cout << "==================================================================" << endl << endl;
	

	if(test == 0)
	{
		cout << "full AMG solving test ..." << endl << endl;
		solve_flexamg(A, vx, b, !full_stencil, false);	
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
		onelevel_amg(A, !full_stencil);
	}
	else if(test == 4)
	{
		cout << "analysing AMG ..." << endl << endl;
		amg_test(A, vx, b, !full_stencil);
	}
	else
		cout << "unsupported test nr. " << test << endl;
	
	//vx.print();
	return 0;
}



#endif