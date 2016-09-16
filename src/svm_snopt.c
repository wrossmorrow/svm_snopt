/*
 *  mnl_snopt.c
 *
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Use these includes on wrmorrow-01,wrmorrow-02
//
#include <vecLib/cblas.h>
#include <vecLib/clapack.h>

// Use this include on wrmorrow-MBP
//
// #include <Accelerate/Accelerate.h>

#include <unuran.h>
#include <unuran_urng_rngstreams.h>

#include "f2c.h"
#include "snfilewrapper.h"
#include "snopt.h"
#include "snoextras.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int   A;
static int * L;
static int   AL;
static int * As;

static int   J;
static int * X;
static int * y;

static double C;

static double * xp0;

static char initcond;

static clock_t tmp_ticks;
static clock_t func_eval_ticks;

void printl( char * filename , char * message )
{
	FILE * fp  = fopen( filename , "a" );
	if( fp != NULL ) { fprintf( fp , message ); fclose( fp ); }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * PRIMAL PROBLEM  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ARGUMENTS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * RETURNS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * NOTES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void svm_primal_size( integer * N , integer * M , integer * Annz , integer * Gnnz )
{
    
	// number of variables (see formulation above for details)
	N[0] = AL + 1 + J; // AL + 1 parameters, J slack variables
    
	// number of constraints (see formulation above for details)
	M[0] = J + A; // one constraint for each observation, and A constraints for uniqueness
	
	// number of nonzeros in linear part:
    //
    //      + C sum_j xi(j)                             J nonzeros, objective is linear in slacks
    //      y(j) x(j,:)' w + y(j) b + xi(j) >= 1        (2+A)J nonzeros in observation constraints
    //      L w = 0                                     AL nonzeros in parameter uniqueness constraints
    //
	Annz[0] = ( A + 3 ) * J + AL;
	
	// number of nonlinear objective gradient and constraint Jacobian non-zeros
	Gnnz[0] = AL; // nonlinear objective is just the squared norm of the parameters
    
}

void svm_primal_bounds( integer N , integer M , double * xLoBnds , double * xUpBnds , double * FLoBnds , double * FUpBnds )
{
	int n;
	
	// no parameter bounds, but slacks are non-negative
	for( n = 0 ; n < AL+1 ; n++ ) { xLoBnds[n] = -1.0e20;   xUpBnds[n] = 1.0e20; }
    for( n = AL+1 ; n < N ; n++ ) { xLoBnds[n] = 0.0;       xUpBnds[n] = 1.0e20; }
	
	// objective has no bounds
	FLoBnds[0] = -1.0e20;	FUpBnds[0] = 1.0e20;
	
	// observation constraints are >= 1, uniqueness constraints are == 0
	for( n = 1   ; n <= J ; n++ ) { FLoBnds[n] = 1.0; FUpBnds[n] = 1.0e20; }
	for( n = J+1 ; n <= M ; n++ ) { FLoBnds[n] = 0.0; FUpBnds[n] = 0.0;    }
	
}

void svm_primal_initcond( integer N , double * x0 )
{
    switch( initcond ) {
            
        case 'p': // prescribed initial parameters
            cblas_dcopy( AL+1 , x0 , 1 , xp0 , 1 );
            cblas_dscal( J , 0.0 , x0 + AL+1 , 1 );
            break;
            
        default:
            cblas_dscal( (int)N , 0.0 , x0 , 1 );
            break;
            
    }
}

void svm_primal_linmap( integer N , integer M , integer Annz , integer * Arows , integer * Acols , double * Adata )
{
    int j, a, l, row, base;
    
    base = 0;
    
    // objective is linear in slack terms, with slope C (penalty)
    
    for( j = 0 ; j < J ; j++ ) {
        Arows[base] = 0;
        Acols[base] = AL+1+j;
        Adata[base] = C;
        base++;
    }
    
    // jth constraint is
    //
    //      y(j) sum_{a,l} x(j,(a,l)) w(a,l) + y(j) b + xi(j) >= 1
    //
    // because of aspect form, we have
    //
    //      sum_a y(j) w( a , l(x(j)) ) + y(j) b + xi(j) >= 1
    //
    // so the entries in A are just -1 and +1... in the right places
    
    for( j = 0 ; j < J ; j++ ) {
        
        row = j + 1;
        
        for( a = 0 ; a < A ; a++ ) {
            Arows[base] = row;
            Acols[base] = As[a] + X[A*j+a]-1; // C-style indexing here, but just here
            Adata[base] = y[j];
            base++;
        }
        
        Arows[base] = row;
        Acols[base] = AL;
        Adata[base] = y[j];
        base++;
        
        Arows[base] = row;
        Acols[base] = AL+1+j;
        Adata[base] = 1.0;
        base++;
        
    }
    
    // we also mandate sums of attribute-level parameters are zero,
    // because this won't affect screening
    printf("As: ");
    for( a = 0 ; a < A ; a++ ) {
        printf(" %i , ",As[a]);
        row = J + 1 + a;
        for( l = 0 ; l < L[a] ; l++ ) {
            Arows[base] = row;
            Acols[base] = As[a] + l;
            Adata[base] = 1.0;
            base++;
        }
    }
    printf("\n");
	
	// enforce FORTRAN-style indexing
	for( j = 0 ; j < Annz ; j++ ) { Arows[j] += 1; Acols[j] += 1; }
}

void svm_primal_structure( integer N , integer M , integer Gnnz , integer * Grows , integer * Gcols )
{
	int i;
	
    for( i = 0 ; i < AL ; i++ ) {
        Grows[i] = 1; Gcols[i] = i+1;
    }
}

int svm_primal_callback_(integer		*Status,	// SNOPT status code
                      integer		*N,			// number of variables
                      doublereal	*x,			// current variable values
                      integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
                      integer		*FN,		// length of the vector of objective and constraint values
                      doublereal	*Fvals,		// values (to be calculated) for objective and constraint values
                      integer		*needG,     // 0 if G(x) not needed, > 0 if it is
                      integer		*Gnnz,		// length of arrays iGvar and jGfun
                      doublereal	*Gvals,		// derivative values (MMF format)
                      char			*cu,		// character workspace
                      integer		*lencu,		// length of character workspace
                      integer		*iu,		// integer workspace
                      integer		*leniu,		// length of integer workspace
                      doublereal	*ru,		// double workspace
                      integer		*lenru )	// length of double workspace
{
	
	if( needF[0] > 0 ) { Fvals[0] = 0.5 * cblas_ddot( AL , x , 1 , x , 1 ); }
	if( needG[0] > 0 ) { cblas_dcopy( AL , x , 1 , Gvals , 1 ); }
    return 0;
	
}

// set-up, solve, and post-process SNOPT formulation
int svm_primal_snopt(double * w ,     // computed parameters: actual parameters
                     double * b ,     // computed parameters: constant term
                     char * optfn ,   // options file name
                     char * logfn )   // log file name
{
	integer			Start = 0; // Cold = 0, Basis = 1, Warm = 2;
	
	integer			l;
	
    integer          N, M, FN;
	
	integer			 Annz;
	double			*Adata;
	integer			*Arows, *Acols;
	
	integer			 Gnnz;
    double			*Gdata;
	integer			*Grows, *Gcols;
	
	integer			*xState;
	double			*x, *xLoBnds, *xUpBnds, *lambdaB;
	
	integer			*FState;
	double			*Fvals, *FLoBnds, *FUpBnds, *lambdaC;
	
	double			ObjAdd; // constant to add to objective
	integer			ObjRow; // row of objective in combined function
	integer			INFO;	//
	
	integer			minrw, miniw, mincw;
	
	// USER workspace
	
	// real (double) workspace
	integer			lenru = 500;
	double			ru[8*500];
	
	// integer workspace
	integer			leniu = 500;
	integer			iu[8*500];
	
	// char workspace
	integer			lencu = 500;
	char			cu[8*500];
	
	// SNOPT workspace
	// we initialize to 500, and re-allocate below
	
	// real (double) workspace
	integer			lenrw = 500;
	double			*rw;
	
	// integer workspace
	integer			leniw = 500;
	integer			*iw;
	
	// char workspace
	integer			lencw = 500;
	char			*cw;
	
	integer			nxName = 1; // do not use variable names
	char			xNames[1*8];
	
	integer			nFName = 1; // do not use constraint names
	char			FNames[1*8];
	
	integer			iPrint = 9; // "unit number for the print file"
	integer			iSumm  = 6; // "unit number for the Summary file"
	
	integer			prnt_len; //
	integer			iSpecs = 4,  spec_len;
	
	char			Prob[200];
	char			printname[200];
	char			specname[200];
	
	integer			nS, nInf, npname = 1;
	double			sInf;
	
	integer			iSum, iPrt, strOpt_len;
	char			strOpt[200];
	
	clock_t ticks;
	
    int a;
	
	FILE * resfp = NULL;
	FILE * logfp = NULL;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	svm_primal_size( &N , &M , &Annz , &Gnnz ); FN = M + 1;
	
	logfp = fopen( logfn , "a" );
	if( logfp != NULL ) {
		fprintf(logfp,"Problem Sizes:\n");
		fprintf(logfp,"  %i variables\n",(int)N);
		fprintf(logfp,"  %i constraints\n",(int)M);
		fprintf(logfp,"  %i linear part nonzeros (%0.1f%% dense)\n",(int)Annz,100.0*(double)Annz/((double)(N*M)));
		fprintf(logfp,"  %i Jacobian nonzeros (%0.1f%% dense)\n",(int)Gnnz,100.0*(double)Gnnz/((double)(N*M)));
		fclose( logfp );
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    xState		 = (integer *)calloc (N, sizeof(integer)); // initialized to zero for no information
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
    lambdaB      = (double *) calloc (N, sizeof(double)); // bounds multipliers
	
	// combined function (objective and constraints)
    Fvals		 = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FState       = (integer *)calloc (FN, sizeof(integer));
    FLoBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FUpBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    lambdaC      = (double *) calloc (FN, sizeof(double)); // constraint multipliers
	
	// linear part of the objective and constraints
    Adata		 = (double *) calloc (Annz, sizeof(double));
    Arows		 = (integer *)calloc (Annz, sizeof(integer));
    Acols		 = (integer *)calloc (Annz, sizeof(integer));
	
	// Jacobian of the nonlinear part of objective and constraints
    Gdata		 = (double *) calloc (Gnnz, sizeof(double));
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	// initial SNOPT workspace; resized below
	cw			 = (char*)   calloc(8*lencw,sizeof(char   ));
	iw			 = (integer*)calloc(  leniw,sizeof(integer));
	rw			 = (double*) calloc(  lenrw,sizeof(double ));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE FORTRAN-STYLE FILE REFERENCES  * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// label spec (options) file, using SNOPT's FORTRAN utilities
    if( optfn != NULL ) { sprintf(specname ,   "%s", optfn); spec_len = strlen(specname); }
    else { sprintf(specname ,   "%s", "opt/primal.spc"); spec_len = strlen(specname); }
	
	// open Print file, using SNOPT's FORTRAN utilities
	sprintf(printname,   "%s", "out/svm_primal_snopt.out"); prnt_len = strlen(printname);
	snopenappend_( &iPrint, printname,   &INFO, prnt_len );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE SNOPT  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "Starting SNOPT...\n" );
	
	sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * SNOPT MEMORY ALLOCATION * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	snmema_(&INFO,
			&FN,
			&N,
			&nxName, &nFName,
			&Annz, &Gnnz,
			&mincw, &miniw, &minrw,
			cw, &lencw,
			iw, &leniw,
			rw, &lenrw,
			8*500);
	
	// if memory was NOT sized successfully,
	if( INFO != 104 ) {
		
		printl( logfn , "WARNING:: SNOPT could not estimate memory requirements correctly.\n" );
		
	} else {
		
		printl( logfn , "SNOPT estimated memory requirements.\n" );
		// fprintf(logfp,"SNOPT estimated memory requirements: %i, %i, %i.\n",(int)mincw,(int)miniw,(int)minrw);
		
		// re-initializing SNOPT workspace, if needed
		
		if( lencw < mincw ) {
			lencw = mincw;
			cw = (char*)realloc(cw, 8*lencw*sizeof(char));
		}
		
		if( leniw < miniw ) {
			leniw = miniw;
			iw = (integer*)realloc(iw, leniw*sizeof(integer));
		}
		
		if( lenrw < minrw ) {
			lenrw = minrw;
			rw = (double*) realloc(rw, lenrw*sizeof(double));
		}
		
		printl( logfn , "Re-initializating SNOPT.\n" );
		// fprintf(logfp,"Re-initializating SNOPT with sizes (%li,%li,%li)\n",lencw,leniw,lenrw);
        
		sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE SNOPT OPTIONS  * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// options
	snfilewrapper_(specname,
				   &iSpecs,
				   &INFO,
				   cw, &lencw,
				   iw, &leniw,
				   rw, &lenrw,
				   spec_len,
				   8*lencw);
	
	if( INFO != 101 ) {
		printl( logfn , "WARNING: trouble reading specs file. Using default options.\n" );
    }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	Start = 0;													// cold start
	
	strcpy(Prob,"SVM");                                         // Problem name
	
	ObjRow = 1; ObjAdd = 0.0;									// objective information
	
	svm_primal_bounds( N , M , xLoBnds , xUpBnds , FLoBnds , FUpBnds ); // define bounds
	printl( logfn , "defined bounds...\n" );
	
	svm_primal_linmap( N , M , Annz , Arows , Acols , Adata );         // define linear part
	printl( logfn , "defined linear map...\n" );
	
	svm_primal_structure( N , M , Gnnz , Grows , Gcols );              // define nonlinear Jacobian structure
	printl( logfn , "defined Jacobian structure...\n" );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    
    int n;
    
    printf("Variables (%li):\n",N);
    for( n = 0 ; n < N ; n++ ) {
        printf("  %0.4f <= x_%i <= %0.4f\n",xLoBnds[n],n+1,xUpBnds[n]);
    }
    
    printf("Combined Map: \n");
    for( n = 0 ; n < FN ; n++ ) {
        printf("  %0.8f <= F(%i) <= %0.8f\n",FLoBnds[n],n+1,FUpBnds[n]);
    }
    
    printf("Linear Part: \n");
    for( n = 0 ; n < Annz ; n++ ) {
        printf("  A[%i] ~ A(%i,%i) = %0.8f\n",n+1,(int)Arows[n],(int)Acols[n],Adata[n]);
    }
    
    printf("Nonlinear Part Jacobian Structure: \n");
    for( n = 0 ; n < Gnnz ; n++ ) {
        printf("  Jdata(%i) ~ J(%i,%i)\n",n+1,(int)Grows[n],(int)Gcols[n]);
    }
    
     /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "\n\nSVM:: Starting...\n" );
	
    // get initial condition
    svm_primal_initcond( N , x );
    
    /*
    // optional derivative check
    if( checkders ) {
        snopta_eval_G_check(N, FN,
                            svm_primal_callback_,
                            Gnnz, Grows, Gcols,
                            x, xLoBnds, xUpBnds,
                            cu, &lencu,
                            iu, &leniu,
                            ru, &lenru);
        exit(1);
    }
     */
    
    // (timed) solve
    ticks = clock();
    
    snopta_(&Start,								// 0: Cold, 1: Basis, 2: Warm
            &FN,								// FN = M + 1
            &N,									// number of variables
            &nxName,							// 1 if no names used, otherwise == N
            &nFName,							// 1 if no names used, otherwise == FN
            &ObjAdd,							// scalar to add to objective (for printing)
            &ObjRow,							// row of the objective in F (FORTRAN-Style)
            Prob,								// problem name
            svm_primal_callback_,				// combined function F(x) needed by snopta_()
            Arows, Acols, &Annz, &Annz, Adata,	// sparse "A" matrix such that F(x) = userfun_(x) + Ax
            Grows, Gcols, &Gnnz, &Gnnz,			// Jacobian structure for G(x) = DF(x)
            xLoBnds, xUpBnds, xNames,			// lower bounds, upper bounds, and names for x variables
            FLoBnds, FUpBnds, FNames,			// lower bounds, upper bounds, and names for F values
            x, xState, lambdaB,					// x values, "states" (see pg. 18), and associated dual variables (multipliers)
            Fvals, FState, lambdaC,				// F values, "states" (see pg. 18?), and associated dual variables (multipliers)
            &INFO,								// result of call to snopta_(). See docs, pg. 19 for details
            &mincw, &miniw, &minrw,				// minimum values of SNOPT workspace sizes for snopta_()
            &nS, &nInf, &sInf,					// see docs, pg. 18 & 20
            cu, &lencu,							// character user workspace
            iu, &leniu,							// integer user workspace
            ru, &lenru,							// real user workspace
            cw, &lencw,							// character SNOPT workspace (at leat 500 + N + NF if names used)
            iw, &leniw,							// integer SNOPT workspace; minimum values given by snmema_()
            rw, &lenrw,							// real SNOPT workspace; minimum values given by snmema_()
            npname, 8*nxName, 8*nFName,
            8*500,
            8*500);
    
    ticks = clock() - ticks;
    
    // (linear) constraint check
    double * Ax;
    
    Ax = ( double * ) calloc ( M , sizeof( double ) );
    for( n = 0 ; n < Annz ; n++ ) {
        if( Arows[n] > 1 ) { // don't include objective
            Ax[ Arows[n]-2 ] += Adata[n] * x[ Acols[n]-1 ];
        }
    }
    
    printf("Constraint Values: \n");
    for( n = 0 ; n < M ; n++ ) {
        printf("  %0.8f <= %0.8f <= %0.8f\n",FLoBnds[n+1],Ax[n],FUpBnds[n+1]);
    }
    
    free( Ax );
    
    // pass back parameter data
    if( INFO == 1 && w != NULL ) {
        cblas_dcopy( AL , x , 1 , w , 1 );
        if( b != NULL ) { b[0] = x[AL]; }
    }
    
    // write out solve status
    logfp = fopen( logfn , "a" );
    if( logfp != NULL ) {
        
        if( INFO != 1 ) {
            fprintf(logfp,"SVM:: SNOPT failed to solve the problem, final status = %i\n", (int)INFO);
        } else {
            fprintf(logfp,"SVM:: SNOPT solved the problem!\n");
        }
        
        fclose( logfp );
    }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
    free(xState);
    free(xLoBnds);
    free(xUpBnds);
	free(lambdaB);
	
	free(Fvals);
    free(FState);
    free(FLoBnds);
    free(FUpBnds);
	free(lambdaC);
	
	free(Adata);
    free(Arows);
    free(Acols);
	
	free(Gdata);
    free(Grows);
    free(Gcols);
	
	// SNOPT workspace
	free(cw);
	free(iw);
	free(rw);
	
	// close print and specs files
	snclose_( &iPrint );
	snclose_( &iSpecs );

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return (int)INFO;
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * DUAL PROBLEM  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ARGUMENTS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * RETURNS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * NOTES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void svm_dual_size( integer * N , integer * M , integer * Annz , integer * Gnnz )
{
    
	// number of variables (see formulation above for details)
	N[0] = J + A + AL; // J + A constraint multipliers, AL "Lagrangian-optimal" parameters
    
	// number of constraints (see formulation above for details)
    //
    //      sum_j lambda(j) y(j) = 0            one constraint
    //      w - X diag(y) lambda - L' mu = 0    AL constraints
    //      L w = 0                             A constraints
    //
	M[0] = 1 + AL + A;
	
	// number of nonzeros in linear part:
    //
    //      - sum_j lambda(j)                   J nonzeros (objective linear in lambda)
    //      sum_j lambda(j) y(j) = 0            J nonzeros (y)
    //      w - X diag(y) lambda - L' mu = 0    2 AL + A J (X diag(y) has A J nonzeros)
    //      L w = 0                             A constraints
    //
	Annz[0] = (2+A) * J + 2 * AL + A;
	
	// number of nonlinear objective gradient and constraint Jacobian non-zeros
	Gnnz[0] = AL; // nonlinear objective is just the squared norm of the parameters
    
}

void svm_dual_bounds( integer N , integer M , double * xLoBnds , double * xUpBnds , double * FLoBnds , double * FUpBnds )
{
	int n;
	
    // lambda multipliers are bounded, everything else unbounded
	for( n = 0 ; n < J ; n++ ) { xLoBnds[n] = 0.0;     xUpBnds[n] = C;      }
    for( n = J ; n < N ; n++ ) { xLoBnds[n] = -1.0e20; xUpBnds[n] = 1.0e20; }
	
	// objective has no bounds
	FLoBnds[0] = -1.0e20;	FUpBnds[0] = 1.0e20;
	
	// all constraints are null-form equality constraints
	for( n = 1 ; n <= M ; n++ ) { FLoBnds[n] = 0.0; FUpBnds[n] = 0.0; }
	
}

void svm_dual_initcond( integer N , double * x0 )
{
    switch( initcond ) {
            
        case 'p': // prescribed initial parameters
            cblas_dcopy( AL+1 , x0 , 1 , xp0 , 1 );
            cblas_dscal( J , 0.0 , x0 + AL+1 , 1 );
            break;
            
        default:
            cblas_dscal( (int)N , 0.0 , x0 , 1 );
            break;
            
    }
}

void svm_dual_linmap( integer N , integer M , integer Annz , integer * Arows , integer * Acols , double * Adata )
{
    int j, a, l, row, base;
    
    base = 0;
    
    // objective is linear in lambda multipliers, with slope -1
    
    for( j = 0 ; j < J ; j++ ) {
        Arows[base] =  0;
        Acols[base] =  j;
        Adata[base] = -1;
        base++;
    }
    
    //      sum_j lambda(j) y(j) = 0            one constraint
    
    for( j = 0 ; j < J ; j++ ) {
        Arows[base] = 1;
        Acols[base] = j;
        Adata[base] = y[j];
        base++;
    }
    
    //      w - X diag(y) lambda - L' mu = 0    AL constraints
    
    for( a = 0 ; a < A ; a++ ) {
        
        for( l = 0 ; l < L[a] ; l++ ) {
            
            // w term
            Arows[base] = 2+a;
            Acols[base] = J + A + As[a] + l;
            Adata[base] = 1.0;
            base++;
            
            // mu term
            Arows[base] = 2+a;
            Acols[base] = J + a;
            Adata[base] = - 1.0;
            base++;
            
        }
        
        // lambda term for this attribute, each product
        for( j = 0 ; j < J ; j++ ) {
            Arows[base] = 2 + As[a] + X[A*j+a]-1;
            Acols[base] = j;
            Adata[base] = - y[j];
            base++;
        }
        
    }
    
    // we also mandate sums of attribute-level parameters are zero,
    // because this won't affect screening
    //
    //      L w = 0
    //
    for( a = 0 ; a < A ; a++ ) {
        row = 2 + AL + a;
        for( l = 0 ; l < L[a] ; l++ ) {
            Arows[base] = row;
            Acols[base] = J + A + As[a] + l;
            Adata[base] = 1.0;
            base++;
        }
    }
	
	// enforce FORTRAN-style indexing
	for( j = 0 ; j < Annz ; j++ ) { Arows[j] += 1; Acols[j] += 1; }
}

void svm_dual_structure( integer N , integer M , integer Gnnz , integer * Grows , integer * Gcols )
{
	int i;
	
    for( i = 0 ; i < AL ; i++ ) {
        Grows[i] = 1; Gcols[i] = J+A+i + 1;
    }
}

int svm_dual_callback_(integer		*Status,	// SNOPT status code
                         integer		*N,			// number of variables
                         doublereal	*x,			// current variable values
                         integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
                         integer		*FN,		// length of the vector of objective and constraint values
                         doublereal	*Fvals,		// values (to be calculated) for objective and constraint values
                         integer		*needG,     // 0 if G(x) not needed, > 0 if it is
                         integer		*Gnnz,		// length of arrays iGvar and jGfun
                         doublereal	*Gvals,		// derivative values (MMF format)
                         char			*cu,		// character workspace
                         integer		*lencu,		// length of character workspace
                         integer		*iu,		// integer workspace
                         integer		*leniu,		// length of integer workspace
                         doublereal	*ru,		// double workspace
                         integer		*lenru )	// length of double workspace
{
	
	if( needF[0] > 0 ) { Fvals[0] = 0.5 * cblas_ddot( AL , x+J+A , 1 , x+J+A , 1 ); }
	if( needG[0] > 0 ) { cblas_dcopy( AL , x+J+A , 1 , Gvals , 1 ); }
    return 0;
	
}

// set-up, solve, and post-process SNOPT formulation
int svm_dual_snopt(double * w ,     // computed parameters: actual parameters
                   double * b ,     // computed parameters: constant term
                   char * optfn ,   // options file name
                   char * logfn )   // log file name
{
	integer			Start = 0; // Cold = 0, Basis = 1, Warm = 2;
	
	integer			l;
	
    integer          N, M, FN;
	
	integer			 Annz;
	double			*Adata;
	integer			*Arows, *Acols;
	
	integer			 Gnnz;
    double			*Gdata;
	integer			*Grows, *Gcols;
	
	integer			*xState;
	double			*x, *xLoBnds, *xUpBnds, *lambdaB;
	
	integer			*FState;
	double			*Fvals, *FLoBnds, *FUpBnds, *lambdaC;
	
	double			ObjAdd; // constant to add to objective
	integer			ObjRow; // row of objective in combined function
	integer			INFO;	//
	
	integer			minrw, miniw, mincw;
	
	// USER workspace
	
	// real (double) workspace
	integer			lenru = 500;
	double			ru[8*500];
	
	// integer workspace
	integer			leniu = 500;
	integer			iu[8*500];
	
	// char workspace
	integer			lencu = 500;
	char			cu[8*500];
	
	// SNOPT workspace
	// we initialize to 500, and re-allocate below
	
	// real (double) workspace
	integer			lenrw = 500;
	double			*rw;
	
	// integer workspace
	integer			leniw = 500;
	integer			*iw;
	
	// char workspace
	integer			lencw = 500;
	char			*cw;
	
	integer			nxName = 1; // do not use variable names
	char			xNames[1*8];
	
	integer			nFName = 1; // do not use constraint names
	char			FNames[1*8];
	
	integer			iPrint = 9; // "unit number for the print file"
	integer			iSumm  = 6; // "unit number for the Summary file"
	
	integer			prnt_len; //
	integer			iSpecs = 4,  spec_len;
	
	char			Prob[200];
	char			printname[200];
	char			specname[200];
	
	integer			nS, nInf, npname = 1;
	double			sInf;
	
	integer			iSum, iPrt, strOpt_len;
	char			strOpt[200];
	
	clock_t ticks;
	
    int a;
	
	FILE * resfp = NULL;
	FILE * logfp = NULL;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	svm_dual_size( &N , &M , &Annz , &Gnnz ); FN = M + 1;
	
	logfp = fopen( logfn , "a" );
	if( logfp != NULL ) {
		fprintf(logfp,"Problem Sizes:\n");
		fprintf(logfp,"  %i variables\n",(int)N);
		fprintf(logfp,"  %i constraints\n",(int)M);
		fprintf(logfp,"  %i linear part nonzeros (%0.1f%% dense)\n",(int)Annz,100.0*(double)Annz/((double)(N*M)));
		fprintf(logfp,"  %i Jacobian nonzeros (%0.1f%% dense)\n",(int)Gnnz,100.0*(double)Gnnz/((double)(N*M)));
		fclose( logfp );
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    xState		 = (integer *)calloc (N, sizeof(integer)); // initialized to zero for no information
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
    lambdaB      = (double *) calloc (N, sizeof(double)); // bounds multipliers
	
	// combined function (objective and constraints)
    Fvals		 = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FState       = (integer *)calloc (FN, sizeof(integer));
    FLoBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FUpBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    lambdaC      = (double *) calloc (FN, sizeof(double)); // constraint multipliers
	
	// linear part of the objective and constraints
    Adata		 = (double *) calloc (Annz, sizeof(double));
    Arows		 = (integer *)calloc (Annz, sizeof(integer));
    Acols		 = (integer *)calloc (Annz, sizeof(integer));
	
	// Jacobian of the nonlinear part of objective and constraints
    Gdata		 = (double *) calloc (Gnnz, sizeof(double));
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	// initial SNOPT workspace; resized below
	cw			 = (char*)   calloc(8*lencw,sizeof(char   ));
	iw			 = (integer*)calloc(  leniw,sizeof(integer));
	rw			 = (double*) calloc(  lenrw,sizeof(double ));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE FORTRAN-STYLE FILE REFERENCES  * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// label spec (options) file, using SNOPT's FORTRAN utilities
    if( optfn != NULL ) { sprintf(specname ,   "%s", optfn); spec_len = strlen(specname); }
    else { sprintf(specname ,   "%s", "opt/dual.spc"); spec_len = strlen(specname); }
	
	// open Print file, using SNOPT's FORTRAN utilities
	sprintf(printname,   "%s", "out/svm_dual_snopt.out"); prnt_len = strlen(printname);
	snopenappend_( &iPrint, printname,   &INFO, prnt_len );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE SNOPT  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "Starting SNOPT...\n" );
	
	sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * SNOPT MEMORY ALLOCATION * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	snmema_(&INFO,
			&FN,
			&N,
			&nxName, &nFName,
			&Annz, &Gnnz,
			&mincw, &miniw, &minrw,
			cw, &lencw,
			iw, &leniw,
			rw, &lenrw,
			8*500);
	
	// if memory was NOT sized successfully,
	if( INFO != 104 ) {
		
		printl( logfn , "WARNING:: SNOPT could not estimate memory requirements correctly.\n" );
		
	} else {
		
		printl( logfn , "SNOPT estimated memory requirements.\n" );
		// fprintf(logfp,"SNOPT estimated memory requirements: %i, %i, %i.\n",(int)mincw,(int)miniw,(int)minrw);
		
		// re-initializing SNOPT workspace, if needed
		
		if( lencw < mincw ) {
			lencw = mincw;
			cw = (char*)realloc(cw, 8*lencw*sizeof(char));
		}
		
		if( leniw < miniw ) {
			leniw = miniw;
			iw = (integer*)realloc(iw, leniw*sizeof(integer));
		}
		
		if( lenrw < minrw ) {
			lenrw = minrw;
			rw = (double*) realloc(rw, lenrw*sizeof(double));
		}
		
		printl( logfn , "Re-initializating SNOPT.\n" );
		// fprintf(logfp,"Re-initializating SNOPT with sizes (%li,%li,%li)\n",lencw,leniw,lenrw);
        
		sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE SNOPT OPTIONS  * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// options
	snfilewrapper_(specname,
				   &iSpecs,
				   &INFO,
				   cw, &lencw,
				   iw, &leniw,
				   rw, &lenrw,
				   spec_len,
				   8*lencw);
	
	if( INFO != 101 ) {
		printl( logfn , "WARNING: trouble reading specs file. Using default options.\n" );
    }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	Start = 0;													// cold start
	
	strcpy(Prob,"SVM");                                         // Problem name
	
	ObjRow = 1; ObjAdd = 0.0;									// objective information
	
	svm_dual_bounds( N , M , xLoBnds , xUpBnds , FLoBnds , FUpBnds ); // define bounds
	printl( logfn , "defined bounds...\n" );
	
	svm_dual_linmap( N , M , Annz , Arows , Acols , Adata );         // define linear part
	printl( logfn , "defined linear map...\n" );
	
	svm_dual_structure( N , M , Gnnz , Grows , Gcols );              // define nonlinear Jacobian structure
	printl( logfn , "defined Jacobian structure...\n" );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     
     printf("Variables (%li):\n",N);
     for( n = 0 ; n < N ; n++ ) {
     printf("  %0.4f <= x_%i <= %0.4f\n",xLoBnds[n],n+1,xUpBnds[n]);
     }
     
     printf("Combined Map: \n");
     for( n = 0 ; n < FN ; n++ ) {
     printf("  %0.8f <= F(%i) <= %0.8f\n",FLoBnds[n],n+1,FUpBnds[n]);
     }
     
     printf("Linear Part: \n");
     for( n = 0 ; n < Annz ; n++ ) {
     printf("  A(%i,%i) = %0.8f\n",(int)Arows[n],(int)Acols[n],Adata[n]);
     }
     
     printf("Nonlinear Part Jacobian Structure: \n");
     for( n = 0 ; n < Gnnz ; n++ ) {
     printf("  Jdata(%i) ~ J(%i,%i)\n",n+1,(int)Grows[n],(int)Gcols[n]);
     }
     
     /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "\n\nSVM:: Starting...\n" );
	
    // get initial condition
    svm_dual_initcond( N , x );
    
    /*
     // optional derivative check
     if( checkders ) {
     snopta_eval_G_check(N, FN,
     svm_primal_callback_,
     Gnnz, Grows, Gcols,
     x, xLoBnds, xUpBnds,
     cu, &lencu,
     iu, &leniu,
     ru, &lenru);
     exit(1);
     }
     */
    
    // (timed) solve
    ticks = clock();
    
    snopta_(&Start,								// 0: Cold, 1: Basis, 2: Warm
            &FN,								// FN = M + 1
            &N,									// number of variables
            &nxName,							// 1 if no names used, otherwise == N
            &nFName,							// 1 if no names used, otherwise == FN
            &ObjAdd,							// scalar to add to objective (for printing)
            &ObjRow,							// row of the objective in F (FORTRAN-Style)
            Prob,								// problem name
            svm_dual_callback_,				// combined function F(x) needed by snopta_()
            Arows, Acols, &Annz, &Annz, Adata,	// sparse "A" matrix such that F(x) = userfun_(x) + Ax
            Grows, Gcols, &Gnnz, &Gnnz,			// Jacobian structure for G(x) = DF(x)
            xLoBnds, xUpBnds, xNames,			// lower bounds, upper bounds, and names for x variables
            FLoBnds, FUpBnds, FNames,			// lower bounds, upper bounds, and names for F values
            x, xState, lambdaB,					// x values, "states" (see pg. 18), and associated dual variables (multipliers)
            Fvals, FState, lambdaC,				// F values, "states" (see pg. 18?), and associated dual variables (multipliers)
            &INFO,								// result of call to snopta_(). See docs, pg. 19 for details
            &mincw, &miniw, &minrw,				// minimum values of SNOPT workspace sizes for snopta_()
            &nS, &nInf, &sInf,					// see docs, pg. 18 & 20
            cu, &lencu,							// character user workspace
            iu, &leniu,							// integer user workspace
            ru, &lenru,							// real user workspace
            cw, &lencw,							// character SNOPT workspace (at leat 500 + N + NF if names used)
            iw, &leniw,							// integer SNOPT workspace; minimum values given by snmema_()
            rw, &lenrw,							// real SNOPT workspace; minimum values given by snmema_()
            npname, 8*nxName, 8*nFName,
            8*500,
            8*500);
    
    ticks = clock() - ticks;
    
    // pass back parameter data
    if( INFO == 1 && w != NULL ) {
        
        cblas_dcopy( AL , x + J + A , 1 , w , 1 );
        
        if( b != NULL ) {
            
            
            
        }
    }
    
    // write out solve status
    logfp = fopen( logfn , "a" );
    if( logfp != NULL ) {
        
        if( INFO != 1 ) {
            fprintf(logfp,"SVM:: SNOPT failed to solve the problem, final status = %i\n", (int)INFO);
        } else {
            fprintf(logfp,"SVM:: SNOPT solved the problem!\n");
        }
        
        fclose( logfp );
    }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
    free(xState);
    free(xLoBnds);
    free(xUpBnds);
	free(lambdaB);
	
	free(Fvals);
    free(FState);
    free(FLoBnds);
    free(FUpBnds);
	free(lambdaC);
	
	free(Adata);
    free(Arows);
    free(Acols);
	
	free(Gdata);
    free(Grows);
    free(Gcols);
	
	// SNOPT workspace
	free(cw);
	free(iw);
	free(rw);
	
	// close print and specs files
	snclose_( &iPrint );
	snclose_( &iSpecs );
    
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return (int)INFO;
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * SVM PROBLEM * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ARGUMENTS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * RETURNS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * NOTES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int svm_snopt(int   extA ,      // number of attributes
              int * extL ,      // A-element vector of number of levels per attribute
              int   extJ ,      // number of profiles
              int * extX ,      // J x A element matrix of profiles (level indexing form)
              int * exty ,      // J-vector of classification responses (-1,1)
              int   extC ,      // soft-margin penalty term
              char  ic ,        // initial condition code
              double * x0 ,     // optional prescribed initial condition for parameters
              char prob ,       // problem code ('d' for "dual", anything else for primal)
              double * w ,      // computed parameters: actual parameters
              double * b ,      // computed parameters: constant term
              char * optfn ,    // options filename
              char * logfn)     // log file name
{
    
    int a, j, info, count;
    
    count = 0;
    for( j = 0 ; j < extJ ; j++ ) {
        if( exty[j] > 0 ) { count++; }
    }
    if( count == 0 || count == extJ ) {
        printf("SVM:: Cannot classify uniform data. Exiting.\n");
        return -1;
    }
    
    // store internal data
    
    A = extA; L = extL;
    As = ( int * ) calloc ( A , sizeof(int) );
    AL = L[0]; As[0] = 0; for( a = 1 ; a < A ; a++ ) { AL += L[a]; As[a] = As[a-1] + L[a-1]; }
    
    J = extJ; X = extX; y = exty;
    
    C = extC;
    
    initcond = 'n'; if( ic == 'p' ) { initcond = 'p'; xp0 = x0; }
    
    // solve problem
    
    switch( prob ) {
        case 'd': info = svm_dual_snopt(w,b,optfn,logfn);
        default:  info = svm_primal_snopt(w,b,optfn,logfn);
    }
    
    // clean up
    
    free( As );
    
    return info;
    
}