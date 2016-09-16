//
//  basictest.c
//  
//
//  Created by W. Ross Morrow on 7/21/14.
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "svm_snopt.h"

int getnext( int A , int * L , int * x )
{
	int a, b;
	
	for( a = 0 ; a < A ; a++ ) {
		if( x[a] < L[a] ) {
			x[a]++; for( b = 0 ; b < a ; b++ ) { x[b] = 1; } return 1;
		}
	}
	
	return 0;
}

int classify( int A , int * As , int * x , double * wT , double bT )
{
    int a;
    double s;
    
    s = wT[x[0]-1]; for( a = 1 ; a < A ; a++ ) { s += wT[As[a]+x[a]-1]; }
    s -= bT;
    
    return ( s >= 0 ? 1 : -1 );
    
}

void printclass( int A , int * x , int y )
{
    int a;
    printf("( %i ",x[0]);
    for( a = 1 ; a < A ; a++ ) {
        printf(", %i ",x[a]);
    }
    printf( " ) , %i \n" , y );
}

double urand()
{
    return ((double)rand()) / ((double)RAND_MAX);
}

int main( int argc , char * argv[] )
{
    int a , l , j , info;
    double W;
    
    int A = 4;
    int L[4] = { 3 , 2 , 5 , 3 };
    int J = 10;
    double C = 1.0e20;
    
    
    int * As;
    int AL;
    
    
    int * X;
    int * y;
    
    double * wT;
    double bT;
    
    double * w;
    double b;
    
    char logfn[1024];
    
    //
    
    As = ( int * ) calloc ( A , sizeof(int) );
    AL = L[0]; As[0] = 0; printf("As: %i , ",As[0]);
    for( a = 1 ; a < A ; a++ ) { AL += L[a]; As[a] = As[a-1] + L[a-1]; printf(" %i , ",As[a]); }
    printf("\n");
    
    // true parameters - random
    
    wT = ( double * ) calloc ( AL , sizeof(double) );
    
    bT = 0;
    for( a = 0 ; a < A ; a++ ) {
        W = 0;
        for( l = 0 ; l < L[a] ; l++ ) {
            wT[As[a]+l] = urand();
            W += wT[As[a]+l];
        }
        // bT += W;
        for( l = 0 ; l < L[a] ; l++ ) {
            wT[As[a]+l] -= W / L[a];
        }
    }
    
    // get profiles and classifications
    
    // this will do * all * profiles
    J = L[0]; for( a = 1 ; a < A ; a++ ) { J *= L[a]; }
    
    // restrict number of profiles here:
    // J = 10;
    
    X = ( int * ) calloc ( A * J , sizeof(int) );
    y = ( int * ) calloc (     J , sizeof(int) );
    
    for( a = 0 ; a < A ; a++ ) { X[a] = 1; }
    y[0] = classify( A , As , X , wT , bT );
    printclass( A , X , y[0] );
    
    for( j = 1 ; j < J ; j++ ) {
        for( a = 0 ; a < A ; a++ ) { X[A*j+a] = X[A*(j-1)+a]; }
        info = getnext( A , L , X + A*j );
        if( info > 0 ) {
            y[j] = classify( A , As , X + A*j , wT , bT );
            printclass( A , X + A*j , y[j] );
        }
    }
    
    // initialize parameter array for successful return
    
    w = ( double * ) calloc ( AL , sizeof(double) );
    
    // define log file
    
    strcpy( logfn , "log/svm_primal_snopt.log" );
    
    // run (primal) SVM
    
    info = svm_snopt( A , L , J , X , y , C , 'r' , NULL , 'p' , w , &b , NULL , logfn );
    
    // print resutls
    
    printf("true: \n");
    for( a = 0 ; a < A ; a++ ) {
        W = wT[As[a]];
        printf( "   %0.6f " , wT[As[a]] );
        for( l = 1 ; l < L[a] ; l++ ) {
            W += wT[As[a]+l];
            printf( ", %0.6f " , wT[As[a]+l] );
        }
        printf(" (%0.6f)\n",W);
    }
    printf( "   %0.6f\n" , bT );
    
    printf("estimated: \n");
    for( a = 0 ; a < A ; a++ ) {
        W = w[As[a]];
        printf( "   %0.6f " , w[As[a]] );
        for( l = 1 ; l < L[a] ; l++ ) {
            W += w[As[a]+l];
            printf( ", %0.6f " , w[As[a]+l] );
        }
        printf(" (%0.6f)\n",W);
    }
    printf( "   %0.6f\n" , b );
    
    // clean up
    
    free( As );
    free( wT );
    free( X );
    free( y );
    free( w );
    
    return 0;
    
}