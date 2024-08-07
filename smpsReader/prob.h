/*
 * prob.h
 *
 *  Created on: Sep 30, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef PROB_H_
#define PROB_H_

#include "smps.h"

#define DECOMPOSE_CHECK /*what is this*/

/* Structure that holds various dimensions of the stage problem */
typedef struct {
	int		rows;			/* number of rows */
	int		cols;			/* number of columns */
	int		intCols;		/* number of integer variables */
	int		binCols;		/* number of binary variables */
	int		prevRows;		/* number of rows in previous stage, is set to 0 for root-stage (master) */
	int		prevCols;		/* number of columns in previous stage, is set to 0 for root-stage (master) */
	int     cntCcols;  		/* number of columns in T which have at least one non-zero element */
	int     cntCrows;  		/* number of rows in T which have at least one non-zero element */
	int		numRV;			/* total number of random variables */
	int 	rvRowCnt;		/* number of rows effected by randomness */
	int 	rvColCnt;		/* number of columns effected by randomness */
    int 	rvaOmCnt;		/* number of RVs in dynamics noise dVector a */
	int 	rvbOmCnt;		/* number of RVs in right-hand side */
	int 	rvyuOmCnt;		/* number of RVs in right-hand side */
	int 	rvylOmCnt;		/* number of RVs in lower bound */
    int 	rvcOmCnt;		/* number of RVs in state-cost coefficients */
    int 	rvdOmCnt;		/* number of RVs in stage-cost coefficients */
    int 	rvAOmCnt;		/* number of RVs in dynamics state-matrix A */
    int 	rvBOmCnt;		/* number of RVs in dynamics decision-matrix B */
	int 	rvCOmCnt;		/* number of RVs in transfer matrix C */
    int 	rvDOmCnt;		/* number of RVs in recourse matrix D */
}numType;


/* structure to hold coordinate information at each stage */
typedef struct {
	iVector	allRVRows;		/* list of all random variable rows */
	iVector	allRVCols;		/* list of all random variable columns */
	iVector	CCols;			/* list of columns in transfer matrix with at least non-zero element in them */
	iVector	CRows;			/* list of rows in transfer matrix with at least non-zero element in them */
	iVector	rvCols;			/* list of all columns in transfer matrix with at least one random element */
	iVector	rvRows;			/* list of all rows with at least one random variable */
	iVector	rvbOmRows;		/* list of all right-hand sides with random variables */
	iVector	rvCOmCols;		/* list of all columns with coefficients with random variables */
	iVector	rvCOmRows;		/* list of all rows with coefficients with random variables */
	iVector	rvdOmCols;		/* list of all columns with cost coefficients with random variables */
	iVector	rvyuOmRows;		/* list of all upper bounds with random variables */
	iVector	rvylOmRows;		/* list of all lower bounds with random variables */
	iVector	rvOffset;		/* Index where the random variable begin - right-hand side, transfer matrix, and cost coefficients. */
}coordType;

typedef struct{
	oneProblem		*sp;			/* structure with complete problem information */
	numType			*num;			/* structure which holds the problem dimensions */
	coordType		*coord;			/* structure which holds the necessary coordinates of the problem */
	sparseVector	*aBar;			/* dynamics dVector a_{t+} */
	sparseVector	*bBar;			/* right-hand side b_t */
	sparseVector	*cBar;			/* state cost coefficients c_t */
	sparseVector	*dBar;			/* objective function coefficients d_t */
	sparseVector    *uBar;          /* the mean value of upper bound of variables */
	sparseVector    *lBar;          /* the mean value of lower bound of variables */
	sparseMatrix	*Abar;			/* dynamics state matrix A_{t+} */
	sparseMatrix	*Bbar;			/* dynamics decision matrix B_{t+} */
	sparseMatrix	*Cbar;			/* transfer matrix C_t */
	sparseMatrix	*Dbar;			/* recourse matrix D_t */
	dVector			mean;			/* Vector of mean values of random variables. */
	int				omegaBeg;		/* Beginning of omega dVector */
	dVector			meanX;			/* Mean value solution */
	double			lb;				/* lower bounds on cost-to-go function */
}probType;

/* subroutines in prob.c */
probType **newProbwSMPS(cString inputDir, cString probName, stocType **stoc, int *numStages);
dVector meanProblem(oneProblem *orig, stocType *stoc);
dVector calcLowerBound(oneProblem *orig, timeType *tim, stocType *stoc);
void freeProbType(probType **prob, int T);
void freeCoordType (coordType *coord);
void printDecomposeSummary(FILE *fptr, cString probName, timeType *tim, probType **prob);


#endif /* PROB_H_ */
