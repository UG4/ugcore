#if 1

/*
 *  ugInterface.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 06.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

/****************************************************************************/
/*									                                        */
/* include files							                                */
/*			  system include files				                            */
/*			  application include files 	                	     	    */
/*																			*/
/****************************************************************************/

/* usage:
 npcreate flexamg $c flexamg;
 npinit flexamg $A stiffMat $x sol $b rhs $m 20 $red 1.0E-8 $abslimit 1.0E-15 $display full;
 
 npexecute flexamg $i $d $r $s $p;
 
 */

// UG 
extern "C"
{

#include "gm.h"
#include "parallel.h"
#include "ls.h"

#include "ugInterface.h"
}


#include <iostream>

#ifndef NDEBUG
#define DEBUG
#endif


#include "algebra.h"
#include <time.h>

int *parentIndex[AMG_MAX_LEVELS];
int flexamg_dimensions=2;

std::vector<postype> positions;
postype GetPosForIndex(int i)
{
	return positions[i];
}

//START_UGDIM_NAMESPACE
USING_UG_NAMESPACES

double ugtime;


string nrstring(double d)
{
	char s[255];
	sprintf(s, "%g", d);
	return string(s);
}

string nrstring(int i)
{
	char s[255];
	sprintf(s, "%d", i);
	return string(s);
}


void writePosToStream(ostream &out)
{
	out << positions.size() << endl;
	for(int i=0; i< positions.size() ; i++)
		out << positions[i] << endl;
}


#define FLEXAMGLIB_TIMING


/****************************************************************************/
/*									    */
/* defines in the following order				       	    */
/*																			*/
/*		  compile time constants defining static data size (e.g. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

REP_ERR_FILE; 
/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef SparseMatrix<myBlockMat> myMatrix;
typedef Vector<myBlockVec> myVector;


typedef struct 
{
	NP_LINEAR_SOLVER ls;   //<! Base class
	
	INT baselevel;
	
	myMatrix A;
	
	myVector b;
	myVector x;
	
	INT offset;
	
	int xcomp;
	int bcomp;
	int Acomp;
} NP_FLEXAMG;

// TODO: change fixed array size to variable see src famg_grid.c TransferVector

#ifdef ModelP
struct s_exchange_struct
{
	int n;
	int ind[MAX_BORDER_CONNECTIONS];
	double val[MAX_BORDER_CONNECTIONS];
};
typedef struct s_exchange_struct exchange_struct;
#endif

/****************************************************************************/
/*																			*/
/* definition of exported global variables				    */
/*									    */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*									    */
/* definition of variables global to this source file only (static!)	    */
/*									    */
/****************************************************************************/

static NP_FLEXAMG *g_np;

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

INT GetOffset(INT nr);

/****************************************************************************/
/*																			*/
/* interface																*/
/*																			*/
/****************************************************************************/


/***************************************************************************/
/** \brief Preprocessing before Solve is executed
 */
/***************************************************************************/

static INT FLEXAMGPreProcess  (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
	NP_FLEXAMG *np = (NP_FLEXAMG*) theNP;
	HEAP *theHeap = MGHEAP(NP_MG(theNP));
	GRID *grid = NP_GRID(theNP,level);
	
	np->baselevel=level;
	
	int *rows, *ncols;
	int *cols;
	double *values;
#ifdef ModelP
	int mpierror, mpiinitialized;
#endif
	VECTOR *v;
	MATRIX *m;  
	
	SHORT *Acomp;
	SHORT *bcomp;
	
	int n;
	int nV;
	INT row,col;
	int hErr;
	int me, procs;
#ifdef ModelP
	MPI_Status status;
#endif
	
#ifdef FLEXAMG_TIMING
	clock_t t1, t2;
	t1 = cclock();
#endif
	
	//Debuggm = 5;   
	
	g_np = np;     
	
	PRINTDEBUG(gm, 2, ("\n FLEXAMGPreProcess {\n"));
	
	/* Count number of (master) vectors */     
	for(v=FIRSTVECTOR(grid), nV = 0; v != NULL; v= SUCCVC(v)) 
	{
		//	nV+=VD_NCMPS_IN_TYPE(b,VTYPE(v));
		nV++;		
	}
	PRINTDEBUG(gm, 2, ("%d Nodes\n", nV));
	
	np->offset=0;
	int i=0;
	for(v=FIRSTVECTOR(grid); v != NULL; v= SUCCVC(v))
	{
		VINDEX(v) = i;
		//i+=VD_NCMPS_IN_TYPE(b,VTYPE(v)); // instead of i++;
		i++;
	}
	
	
	
 	/* ALLOCATE MEMORY */
	
	PRINTDEBUG(gm, 2, ("Allocate Memory\n"));
	
	
	// FLEXAMG: CREATE A, b and x
	
	np->A.create(nV, nV);
	np->b.create(nV);
	np->x.create(nV);

	// UG A --> FLEXAMG A
	//-------------------
	PRINTDEBUG(gm, 2, ("UG A --> FLEXAMG A.\n"));
	
	
	/* Copy SparseMatrix row to FLEXAMG for all (master) nodes
	 (row-wise, always TMP_MAX_VALUES at once)
	 */
	
	PRINTDEBUG(gm, 2, (" get\n\n"));

	
	std::vector<myMatrix::connection> cons;
	positions.resize(nV);
		
	n=0;  // total number of entries 
	nV=0; // row counter 
	
	for(v=FIRSTVECTOR(grid); v != NULL; v= SUCCVC(v))
	{
	  
	  	INT vindex=VINDEX(v);
		short vtype=VTYPE(v);
		short vncmps=VD_NCMPS_IN_TYPE(x, vtype);
	
		cons.clear();
		
	    for(m=VSTART(v); m!=NULL; m=MNEXT(m))
		{ 
			// copy row of A_ij
			VECTOR *dest = MDEST(m);
			short dtype  = VTYPE(dest);
			short dncmps=VD_NCMPS_IN_TYPE(b, dtype);
			Acomp = MD_MCMPPTR_OF_RT_CT(A, vtype, dtype);  
			
			myMatrix::connection con;
			con.iIndex = VINDEX(dest);
			con.dValue = 0.0;
			
			
			ASSERT2(vncmps == UNKNOWN_NR, "vncmps = " << vncmps << ", but UNKNOWN_NR = " << UNKNOWN_NR << "!");
			ASSERT2(dncmps == UNKNOWN_NR, "dncmps = " << dncmps << ", but UNKNOWN_NR = " << UNKNOWN_NR << "!");
			
			for (int r=0; r<vncmps; ++r)
				for (int c=0; c<dncmps; ++c)
				  setAt(con.dValue, r, c, MVALUE(m, Acomp[r*dncmps + c]));					
			if(con.dValue != 0.0)
				cons.push_back(con);
		}
		
		np->A.setMatrixRow(vindex, &cons[0], cons.size());		
		positions[vindex].x = CVECT(MYVERTEX(VMYNODE(v)))[0];
		positions[vindex].y = CVECT(MYVERTEX(VMYNODE(v)))[1];
		
	}

	PRINTDEBUG(gm, 2, ("stored all in flexamg structures\n"));
	
	//np->A.finish();
	//np->A.print();
	//if(nV < 100000) 	np->A.writeToFile("/Users/mrupp/matrices/A.mat");

	PRINTDEBUG(gm, 2, ("} FLEXAMGPreProcess.\n"));
	
#ifdef FLEXAMG_TIMING
	t2 = clock();
	ugtime = ((double)t2-t1)/CLOCKS_PER_SEC;
#endif
	// return 
	return(0);
}

/***************************************************************************/
/** \brief Copies x and b to FLEXAMG and executes the solver
 */
/***************************************************************************/
static INT FLEXAMGSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
	printf("yes, its me. the FLEXAMG Solver.\n");
	
  NP_FLEXAMG *np = (NP_FLEXAMG*) theNP;
	GRID *grid = NP_GRID(theNP, level);

	VECTOR *v;
	int nV;   
	double *values;
	int *indices;
	INT i; int hErr, n;

	HEAP *theHeap = MGHEAP(NP_MG(theNP));

	VECDATA_DESC *temp;
#ifdef FLEXAMG_TIMING
	clock_t t1, t2; double hypretime;
	t1 = clock();
#endif
	DOUBLE predefect=0.0, postdefect=0.0;
	PRINTDEBUG(gm, 2, ("FLEXAMGSolver {\n"));

	g_np = np;     

	// create a unique copy of defect
	temp=NULL;
	if (AllocVDFromVD(NP_MG((NP_BASE*) theNP),level,level,x,&temp))
	   	REP_ERR_RETURN (1);

	if (dcopy(NP_MG((NP_BASE*)theNP),level,level,ALL_VECTORS,temp,b)!= NUM_OK)							
		REP_ERR_RETURN (1);

	if (dnrm2(NP_MG((NP_BASE*)theNP),level,level,ON_SURFACE,b,&predefect))  
		     return 1;
	//printf("sqrt<b,b>=%e.\n",predefect);
	

	/* Copy b : UG --> FLEXAMG */
	PRINTDEBUG(gm, 2, ("UG b --> FLEXAMG b\n"));
	nV=0;
	/* copy blocks of TMP_MAX_VALUES at once*/	
	for(v=FIRSTVECTOR(grid); v != NULL; v= SUCCVC(v))
	{
		short vtype = VTYPE(v);
		const short vncmps = VD_NCMPS_IN_TYPE(b, vtype);
		short c;
		for (c=0; c<vncmps; c++)
		{
			double d = VVALUE(v, VD_CMP_OF_TYPE(b, vtype, c));
			setAt(np->b[VINDEX(v)], c, d);
		}		 
		nV++;
	}


	IFDEBUG(gm, 2)
		//np->b.print();
	ENDDEBUG;

	cout.flush();


	PRINTDEBUG(gm, 2, ("Execute Solver:\n"));

	amg<myMatrix, myVector> AMG;
	//AMG.setAggressiveCoarsening_A_2();
	//AMG.onlyOneLevel(np->A);
	

	np->x = 0.1;

	AMG.amgTest(np->A,np->x,np->b);
	//AMG.onlyOneLevel(np->A);
	
	//LinearSolver(np->x, np->A, np->b, AMG, 40, 1e-12);
	//CG(np->x, np->A, np->b, AMG, 100);
	

	PRINTDEBUG(gm, 2, ("finished!!!!\n"));
	
	if(nV < 100000)
		AMG.writeMatrices("/Users/mrupp/matrices/SparseMatrix");


//	np->x.print();

	printf("__new:\n");
	/* Copy x : FLEXAMG --> UG */
	for(v=FIRSTVECTOR(grid); v != NULL; v= SUCCVC(v))
	{
		short vtype = VTYPE(v);
		const short vncmps = VD_NCMPS_IN_TYPE(b, vtype);
		short c;
		for (c=0; c<vncmps; c++)
		{
			double d;
			d = getAt(np->x[VINDEX(v)], c);

			//printf("%d.%d: %g + %g = %g\n", VINDEX(v), c, VVALUE(v, VD_CMP_OF_TYPE(x,vtype,c)), d,  VVALUE(v, VD_CMP_OF_TYPE(x,vtype,c))+d);
			VVALUE(v, VD_CMP_OF_TYPE(temp,vtype,c)) = d; 
			VVALUE(v, VD_CMP_OF_TYPE(x,vtype,c)) += d; // we are solving the defect equation!
			
		}		    
	}

			
	       
		
#ifdef FLEXAMG_TIMING
	t2 = clock();
	ugtime += ((double)t2-t1)/CLOCKS_PER_SEC;
#endif	
	return (NUM_OK);
}

/***************************************************************************/
/** \brief Cleanup (Destroy FLEXAMG Matrices)
 */
/***************************************************************************/

static INT FLEXAMGPostProcess (NP_LINEAR_SOLVER *theNP, INT level,VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{     
	NP_FLEXAMG *np = (NP_FLEXAMG*) theNP;
	PRINTDEBUG(gm, 2, ("FLEXAMGPostProcess { \n"));
	
	/* Free SparseMatrix/vector data structures */
	
	PRINTDEBUG(gm, 2, ("} FLEXAMGPostProcess \n"));
	
	return(0);
}

/****************************************************************************/
/** \brief Read command line integers
 
 \param name - name of the argument
 \param a - integer value 1
 \param n - integer value 2
 \param argc - argument counter 
 \param argv - argument vector
 
 This function reads command strings and returns an integer value in 'a' and 'b'.
 
 \return <ul>
 <li>  0 if the argument was found and a value could be read </li>
 <li>  1 else. </li>
 </ul>
 */
/****************************************************************************/




static INT FLEXAMGInit (NP_BASE *theNP, INT argc , char **argv)
{  
	NP_FLEXAMG *np = (NP_FLEXAMG *) theNP;
	
	return (NPLinearSolverInit(&np->ls,argc,argv));
	
}

/***************************************************************************/
/***************************************************************************/
//! Display routine for numproc
static INT FLEXAMGDisplay (NP_BASE *theNP)
{
	NP_FLEXAMG *np;
	np = (NP_FLEXAMG *) theNP;
	printf("FLEXAMGSolver!!! there can be only one :)))");
	
	return(0); // (FLEXAMG_Output(theNP));
}

/***************************************************************************/
/***************************************************************************/


INT FLEXAMGExecute (NP_BASE *theNP, INT argc , char **argv)
{
	return(1);
}

/***************************************************************************/
/***************************************************************************/

static INT LinearDefect (NP_LINEAR_SOLVER *theNP, INT level,
						 VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
						 INT *result)
{
	NP_FLEXAMG *np;	
	INT bl; // \todo needs to be introduced
	
	np = (NP_FLEXAMG *) theNP;
	bl = MIN(FULLREFINELEVEL(NP_MG(theNP)),MAX(0,np->baselevel));
	if (dmatmul_minus(NP_MG(theNP),bl,level,ON_SURFACE,b,A,x) != NUM_OK)  
		NP_RETURN(1,result[0]);
	
	return (*result);
}

/***************************************************************************/
/***************************************************************************/

static INT LinearResiduum (NP_LINEAR_SOLVER *theNP, INT bl, INT level,
						   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
						   LRESULT *lresult)
{
	NP_LINEAR_SOLVER  *np;	
	
	np = (NP_LINEAR_SOLVER  *) theNP;
#ifdef ModelP
	if (a_vector_collect(NP_MG(theNP),bl,level,b))  
		NP_RETURN(1,lresult->error_code);
#endif
	if (dnrm2x(NP_MG(theNP),bl,level,ON_SURFACE,b,lresult->last_defect))  
		NP_RETURN(1,lresult->error_code);
	
	return(0);	  
	
}

/***************************************************************************/
//! Numproc constructor
/***************************************************************************/
static INT FLEXAMGConstruct (NP_BASE *theNP)
{
	NP_LINEAR_SOLVER *np;
	
	theNP->Init     = FLEXAMGInit;
	theNP->Display  = FLEXAMGDisplay;
	theNP->Execute  = NPLinearSolverExecute;
	
	np = (NP_LINEAR_SOLVER *) theNP;
	np->PreProcess  = FLEXAMGPreProcess;
	np->Defect      = LinearDefect;
	np->Residuum    = LinearResiduum;
	np->Solver      = FLEXAMGSolver;
	np->PostProcess = FLEXAMGPostProcess;
	
	return(0);
}


/****************************************************************************/
//! InitFLEXAMG	- Init this file
/****************************************************************************/
extern "C" INT InitFLEXAMG ()
{
	if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".flexamg",sizeof(NP_FLEXAMG),FLEXAMGConstruct))
		REP_ERR_RETURN (__LINE__);
	return(0);	
}

//END_UGDIM_NAMESPACE

#endif
