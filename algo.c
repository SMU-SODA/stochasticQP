/*
 * algo.c
 *
 *  Created on: Aug 30, 2021
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "stochasticQP.h"

extern configType config;

int runAlgo (probType **prob, stocType *stoc, cellType* cell) {

	
	double subset = 0.1 * cell->omega->cnt; /*initialize the number of samples you want to take from Omega in each iteration*/
	/*Assign memory to the basket of fixed sections of alpha and beta, if you are going to use dual-based algorithm*/
	double SigmaSize = subset * config.MAX_ITER;
	sigmaType* sigma; /* Sigma is a collection of fixed parts of alpha and beta which is independent of the observation */
	sigma = (sigmaType*)mem_malloc(sizeof(sigmaType));
	sigma->vals = (pixbCType**)arr_alloc(SigmaSize, pixbCType*);
	sigma->cnt = 0;
	for (int i = 0; i < SigmaSize; i++) {
		sigma->vals[i] = (pixbCType*)mem_malloc(sizeof(pixbCType));
	}
	for (int i = 0; i < SigmaSize; i++) {
		sigma->vals[i]->piCar = (double*)arr_alloc(prob[0]->num->cols + 1, double);
	}

	/* Assign memory to lambda structure which is the set of dual solutions we obtained so far. pi is related to equality constraints, mu2 
	corresponds to upper bounds and mu3 corresponds to lower bounds */

	lambdaType* lambda;
	lambda = (lambdaType*)mem_malloc(sizeof(lambdaType));
	lambda->pi = (double**)arr_alloc(SigmaSize, double*);
	lambda->mu2 = (double**)arr_alloc(SigmaSize, double*);
	lambda->mu3 = (double**)arr_alloc(SigmaSize, double*);
	for (int i = 0; i < SigmaSize; i++) {
		lambda->pi[i] = (double*)arr_alloc(prob[1]->num->rows + 1, double);
		lambda->mu2[i] = (double*)arr_alloc(prob[1]->num->cols + 1, double);
		lambda->mu3[i] = (double*)arr_alloc(prob[1]->num->cols + 1, double);
	}
	lambda->cnt = 0;

	/* assign memory to deta structure, this will record the deltaAlpha and deltaBetha associated with each observation*/
	deltaType* delta;
	delta = (deltaType*)mem_malloc(sizeof(deltaType));
	delta->vals = (lambdadeltaType***)arr_alloc(cell->omega->cnt, lambdadeltaType**);
	for (int i = 0; i < cell->omega->cnt; i++) {
		delta->vals[i] = (lambdadeltaType**)arr_alloc(SigmaSize, lambdadeltaType*);
	}
	for (int i = 0; i < cell->omega->cnt; i++) {
		for (int j = 0; j < SigmaSize; j++) {
			delta->vals[i][j] = (lambdadeltaType*)mem_malloc(sizeof(lambdadeltaType));
			delta->vals[i][j]->dbeta=(sparseVector*)mem_malloc(sizeof(sparseVector));
			delta->vals[i][j]->dbeta->val = (double*)arr_alloc(prob[0]->num->cols + 1, double); /*TO DO: Change it to a reasonable size*/
			delta->vals[i][j]->dbeta->col = (int*)arr_alloc(prob[0]->num->cols + 1,int);
		}
	}



	oneCut *cut = NULL;
	while ( cell->k < 25) {
		cell->k++;
		/* 1. Check optimality */

		/* 2. Switch between algorithms to add a new affine functions. */
		switch (config.ALGOTYPE) {
		case 0:
			cut = fullSolveCut(prob[1], cell, stoc, cell->candidX);
			if ( cut == NULL ) {
				errMsg("algorithm", "runAlgo", "failed to create the cut using full solve", 0);
				goto TERMINATE;
			}
			break;
		case 1:
			dualSolve(prob, cell, stoc,  sigma,  delta, lambda, cell->candidX, subset);
			break;
		case 2:
			partSolve();
			break;
		default:
			errMsg("ALGO", "main", "Unknown algorithm type", 0);
			goto TERMINATE;
		}

#if defined(ALGO_CHECK)
		double objEst = vXvSparse(cell->candidX, prob[0]->dBar) + cut->alpha - vXv(cut->beta, cell->candidX, NULL, prob[0]->num->cols);
		printf("\tCandidate estimate = %lf\n", objEst);
#endif

		/* 3. Add the affine function to the master problem. */
		/* 3a. Add cut to the set */
		cell->cuts->vals[cell->cuts->cnt] = cut;
		cell->cuts->cnt++;

		/* 3b. Update the master problem by adding the cut*/
		if ( addCut2Solver(cell->master->model, cut, prob[0]->num->cols) ) {
			errMsg("solver", "fullSolve", "failed to add cut to the master problem", 0);
			return 1;
		}

#if defined(WRITE_FILES)
		char str[BLOCKSIZE];
		sprintf(str, "cellMaster.lp");
		if (writeProblem(cell->master->model, str)) {
			sprintf(str, "failed to write the stage problem %s after adding the cut", cell->subprob->name);
			errMsg("solver", "fullSolve", str, 0);
			return 1;
		}
#endif

		/* 4. Solve the master problem. */
		if ( solveProblem(cell->master->model) ) {
			errMsg("solver", "fullSolve", "failed to solve the master problem", 0);
			return 1;
		}

		double objvalmaster;
		objvalmaster = getObjective(cell->master->model);

		if (getPrimal(cell->master->model, cell->candidX, 0, prob[0]->num->cols) ) {
			errMsg("solver", "fullSolve", "failed to obtain the candidate solution", 0);
			return 1;
		}


	}

	return 0;

	TERMINATE:
	return 1;
}//END runALgo()

int addCut2Solver(modelPtr *model, oneCut *cut, int lenX) {
	iVector rmatind;

	rmatind = (iVector) arr_alloc(lenX+1, int);

	rmatind[0] = lenX;
	for (int i = 1; i <= lenX;++i) {
		rmatind[i] = i-1;
	}

	static int cummCutNum = 0;

	/* Set up the cut name */
	sprintf(cut->name, "cut_%04d", cummCutNum++);

	/* Add a new linear constraint to a model. */
	if( addRow(model, lenX, cut->alpha, GE, rmatind, cut->beta, cut->name) ) {
		errMsg("solver", "addCut2Solver", "failed to addrow", 0);
		return 1;
	}

	mem_free(rmatind);

	return 0;
}//END addCuts2Solver()

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta dVector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->name = (cString) arr_alloc(NAMESIZE, char);
	cut->rowNum = -1;

	cut->beta = (dVector) arr_alloc(numX + 1, double);
	cut->beta[0] = 1.0;
	cut->alpha = 0.0;

	return cut;
}//END newCut

