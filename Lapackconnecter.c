/*
 * test.cpp
 *
 *  Created on: Jan 23, 2023
 *      Author: george
 */

#include <lapacke.h>
#include <stdio.h>
#include "utilities.h"


//double* testLA (double* A, int S);

void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_vector( char* desc, lapack_int n, lapack_int* vec );

double* testLA (double* A , int S)
{
	int N =  S;
	double* B = (double*)arr_alloc(N*N , double);

	for (int n = 0 ; n < N*N ; n++){
		B[n] = A[n];
	}

	int pivotArray[N]; //since our matrix has three rows
	int errorHandler;
	double lapackWorkspace[N*N];

    int NN = N*N;
	dgetrf_(&N, &N, B, &N, pivotArray, &errorHandler);
	if(errorHandler){
	    	printf ("dgetri eh, %d, should be zero\n", errorHandler);
	 }

	 dgetri_(&N, B, &N, pivotArray, lapackWorkspace, &NN, &errorHandler);
	 if(errorHandler){
	    	printf ("dgetri eh, %d, should be zero\n", errorHandler);
	  }

	 //free(lapackWorkspace);
	 return B;
} /* End of LAPACKE_dgels Example */
/* Auxiliary routine: printing a matrix */
void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );

        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i*ldm+j] );
                printf( "\n" );
        }
}


/* Auxiliary routine: printing a matrix */
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );

        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i+j*ldm] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_vector( char* desc, lapack_int n, lapack_int* vec ) {
        lapack_int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6" LAPACK_IFMT, vec[j] );
        printf( "\n" );
}


