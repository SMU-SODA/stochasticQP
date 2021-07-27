/*
 * smps.h
 *
 *  Created on: Dec 30, 2020
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef SMPSREADER_SMPS_H_
#define SMPSREADER_SMPS_H_

#include "../solverUtilities/utilities.h"
#include "../solverUtilities/solver_gurobi.h"


#undef SMPS_CHECK

typedef struct{
	int		type;			/* type of problem: LP, QP, MIP or MIQP */
	void	*model;			/* problem pointer to be used by solver */
	int 	objSense;		/* Sense of the objective */
	cString	name;			/* name of the problem */

	int		mac;			/* number of columns */
	int 	mar;			/* number of rows */
	int		numBin;			/* number of binary variables in the problem */
	int		numInt;			/* number of integer variables  (includes both binary and general integer variables) in the problem */
	int		numnz;			/* number of non-zero elements in constraint matrix */



	dVector	objx;			                     /* objective function linear coefficients */
	sparseMatrix* objQ;                          /* objective function Hessian Matrix */




	dVector	bdl;			/* lower bound */
	dVector	bdu;			/* upper bound */
	cString	ctype;			/* type of decision variables: 'C' continuous, 'B' binary, 'I' general integer, 'S' semi-continuous, 'N' semi-integer */

	dVector	rhsx;			/* right-hand side */
	cString	senx;			/* constraint sense */

	iVector	matbeg;			/* sparse matrix representation: column beginning */
	iVector	matcnt;			/* sparse matrix representation: number of non-zero entries in a column */
	iVector	matind;			/* sparse matrix representation: rows with non-zero entries */
	dVector	matval;			/* sparse matrix representation: non-zero coefficients of the matrix */

	cString	objname;		/* objective function name */
	cString	*rname;			/* dVector of row names */
	cString	*cname;			/* dVector of column names */

	/* These don't seem necessary for Gurobi. */
	//	int		rstorsz;		/* memory size for storing row names */
	//	cString	rstore;			/* row names cString */
	//	int		cstorsz;		/* memory size for storing column names */
	//	cString	cstore;			/* column name cString */

	int		macsz;			/* extended column size ?///*/
	int		marsz;			/* extended row size???????? */
	int		matsz;			/* extended matrix size?????? */
}oneProblem;

typedef struct {
	int			type;			/* type of time file declaration, 0 for implicit and 1 for explicit */
	cString		probName;		/* name of the problem as read from time file */
	int			numStages;	    /* number of stages in the problem */
	cString		*stgNames;		/* unique cStrings to identify stages*/
	iVector		row;			/* a list of row names which mark the beginning of a new stage */
	iVector		col;			/* a list of column names which mark the beginning of a new stage */

	int			numRows;		/* used with explicit time file declaration only, set to numStages in implicit declaration */
	iVector		rowStg;			/* used with explicit time file declaration only */
	int			numCols;		/* used with explicit time file declaration only, set to numStages in implicit declaration */
	iVector		colStg;			/* used with explicit time file declaration only */
}timeType;

/* The statistical model structure is designed to hold information about processes which can be described as linear transformation of the past
 * values/observations, and past and current residual/error terms. A classical example of this type of model is the ARMA(p,q) model. The
 * description of the elements of this structure are written using ARMA as reference */
typedef struct {
	int				p;			/* autoregression order */
	int				q;			/* moving-average order */
	int				N;			/* dimension of time series */
	int				M; 			/* dimension of the residual process/noise */
	dVector			muEps;		/* Mean of the residual process/noise */
	sparseMatrix	*cvEps;		/* Covariance matrix of the residual process/noise */
	sparseMatrix	**AR;		/* autoregression coefficients */
	sparseMatrix	**MA;		/* moving-average coefficients */
	dVector			*eta;		/* trend time series */
	dVector			*sigma;		/* seasonal time series */
}statModel;

typedef struct {
	cString	type;				/* type of stocType being used */ 
	bool	sim;				/* set to TRUE if an external simulator is used */
	int		numOmega; 			/* number of stochastic elements stored in structure */
	int		numGroups;          /*number of groups of rv that need to be simulated seprately*/
	iVector	row; 				/* row number array in the original problem; -1 indicates objective function */
	iVector	col; 				/* column number array in the original problem; -1 indicates right-hand side */
	iVector	numVals;			/* number of realization for each random variable */
	dVector	*vals; 				/* indexed array of discrete realizations of random variable */
	dVector	*probs;				/* indexed array of probabilities associated with discrete realizations*/
	iVector	numPerGroup;        /*how many are in each group*/
	iVector	groupBeg;      
	dVector	mean;         		/* mean of each rv */
	statModel *mod;
}stocType;

/* subroutines in smps.c */
int readFiles(cString inputDir, cString probName, oneProblem **orig, timeType **tim, stocType **stoc);
oneProblem *readCore(cString inputDir, cString probName);
timeType *readTime(cString inputDir, cString probName, oneProblem *orig);
stocType *readStoc(cString inputDir, cString probName, oneProblem *orig, timeType *tim);
int readIndep(FILE *fptr, cString *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc, cString **rvRows, cString **rvCols);
int readIndepDiscrete(FILE *fptr, cString *fields, int maxOmegas, int maxVals, cString **rvRows, cString **rvCols, oneProblem *orig, stocType *stoc);
int readNormal(FILE *fptr, cString *fields, int maxOmegas, cString **rvRows, cString **rvCols, oneProblem *orig, stocType *stoc);
int readBlocks(FILE *fptr, cString *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc, cString **rvRows, cString **rvCols);
int readOneBlock(FILE *fptr, cString *fields, oneProblem *orig, int maxOmegas, int maxVals, bool origRV, stocType *stoc);
int readLinTrans(FILE *fptr, cString *fields, oneProblem *orig, stocType *stoc, int maxOmegas, cString **rvRows, cString **rvCols);
int readARMA(FILE *fptr, cString *fields, oneProblem *orig, stocType *stoc, int maxOmegas);
int readScenarios(FILE *fptr, cString *fields, oneProblem *orig, timeType *tim, int maxOmegas, int maxVals, stocType *stoc);

void freeOneProblem(oneProblem *p);
void freeTimeType(timeType *tim);
void freeStocType(stocType *stoc);
void freeStatModel(statModel *model);

/* subroutines in rvgen.c */
int generateOmegaIdx(stocType *stoc, long long *seed);
void generateOmega(stocType *stoc, dVector observ, double minVal, long long *seed, FILE **fptr);
void generateBlocks(stocType *stoc, dVector observ, int groupID, long long *seed);
void generateIndep(stocType *stoc, dVector observ, int groupID, long long *seed);
void generateLinTran(stocType *stoc, dVector observ, int groupID, double minVal, long long *seed);
int normal(dVector mu, dVector stdev, int numOmega, dVector observ, long long *seed);
int weibull(double scaleParam, double shapeParam, int numOmega, dVector observ, long long *seed);
double scalit(float lower, float upper, long long *RUN_SEED);
double randUniform();
int randInteger(int iMax);
int setupSAA(stocType *stoc, cString fname, long long *seed, dVector *simObservVals, dVector probs,
		int *numObs, int desiredSampleSize, double TOLERANCE);
int readSimData(cString fname, dVector *simObservVals, int numRV, int *numSamples);
int readSimLine(FILE **fid, dVector observ, int numRV, bool simulate);
void computeSampleMean(dVector *vals, iVector weights, int numRV, int sampleSize, int numObs, dVector sampleMean);



#endif /* SMPSREADER_SMPS_H_ */
