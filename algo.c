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
    double opteta = 0;

    double status = 0;
    double quadcand = 0;
    double quadincum = 0;
    double fxk = 0;
    double fxk1 = 0;
    double fwk = 0;
    double fwk1 = 0;
    cell->incumbEst = -INFINITY;
	oneCut *cut = NULL;

	clock_t StartCut;
	clock_t EndCut;
	clock_t StartMas;
	clock_t EndMas;
	dVector temp1;
	double incumbObj;
	double subset = config.SAMPLE_FRACTION * cell->omega->cnt; /*initialize the number of samples you want to take from Omega in each iteration*/
	clock_t tStart = clock();
	double Tempobj = -100000;
	cell->obj = -50000;
	int* index;
	index = (int*) arr_alloc(prob[0]->num->cols + 1, int);

	double* meanx = (double*) arr_alloc(prob[0]->num->cols + 1, double);
	for (int i = 1; i <= prob[0]->num->cols; i++) {
		index[i] = i;
		meanx[i]=cell->candidX[i];
	}
	cell->k  = 0;

    /* update the master problem with the initial incumbent*/

	changeQPrhs(prob[0], cell, cell->incumbX);
	changeQPbds(cell->master->model , prob[0]->num->cols , prob[0]->lBar, prob[0]->uBar, cell->incumbX);
	changeQPcoef(cell->master->model , prob[0] ,  prob[0]->num->cols , prob[0]->dBar,  cell->incumbX);
    double* candidatesol = (double*) arr_alloc(prob[0]->num->cols + 1, double);

	if(prob[0]->sp->objQ != NULL){temp1 = vxMSparse(cell->incumbX, prob[0]->sp->objQ, prob[0]->num->cols);}
	else{temp1 = (double*)arr_alloc(prob[0]->num->cols + 1 , double);}

	incumbObj = vXv(temp1, cell->incumbX, index, prob[0]->num->cols) + vXvSparse( cell->incumbX , prob[0]->dBar);

    mem_free(temp1);
	while (cell->obj - Tempobj > 0.001 || cell->numit < 35){
		cell->numit++;

		/* 1. Check optimality */

		/* 2. Switch between algorithms to add a new affine functions. */
		switch (config.ALGOTYPE) {
		case 0:
			if(cell->numit > 1){
				for (int i = 1; i <= prob[0]->num->cols; i++) {
								candidatesol[i] = cell->candidX[i] + cell->incumbX[i] ;
							}
			}
			else{for (int i = 1; i <= prob[0]->num->cols; i++) {
				candidatesol[i] = cell->candidX[i]  ;
			}}
			StartCut = clock();
			cut = fullSolveCut(prob[1], cell, stoc, candidatesol);
			EndCut = clock();
			cell->Tcut = cell->Tcut + (EndCut - StartCut);
			//printf("\n %f", ((EndCut - StartCut) / CLOCKS_PER_SEC));
			//printf("\n %f", cell->Tcut / CLOCKS_PER_SEC);
			if ( cut == NULL ) {
				errMsg("algorithm", "runAlgo", "failed to create the cut using full solve", 0);
				goto TERMINATE;
			}
			break;

		case 1:
			StartCut = clock();
			cut = dualSolve(prob[1], cell, stoc, cell->candidX, subset);
			EndCut = clock();
			cell->Tcut = cell->Tcut + (EndCut - StartCut);
			if (cut == NULL) {
				errMsg("algorithm", "runAlgo", "failed to create the cut using dual solve", 0);
				goto TERMINATE;
			}
			break;
		case 2:
			/* calculate delta x*/

			if(cell->numit > 1){
				for (int i = 1; i <= prob[0]->num->cols; i++) {
								candidatesol[i] = cell->candidX[i] + cell->incumbX[i] ;
							}
			}
			else{for (int i = 1; i <= prob[0]->num->cols; i++) {
				candidatesol[i] = cell->candidX[i]  ;
			}}

			StartCut = clock();
			cut = partSolve(prob[1],  cell,  stoc, candidatesol, meanx, subset);
			EndCut = clock();
			cell->Tcut = cell->Tcut + (EndCut - StartCut);
			//printf("\n %f", ((EndCut - StartCut)/ CLOCKS_PER_SEC));
			//printf("\n %f", cell->Tcut  / CLOCKS_PER_SEC);
			if (cut == NULL) {
				errMsg("algorithm", "runAlgo", "failed to create the cut using partition-based solve", 0);
				goto TERMINATE;
			}
			break;
		default:
			errMsg("ALGO", "main", "Unknown algorithm type", 0);
			goto TERMINATE;
		}

#if defined(ALGO_CHECK)
		double objEst = vXvSparse(cell->candidX, prob[0]->dBar) + cut->alpha - vXv(cut->beta, cell->candidX, NULL, prob[0]->num->cols); // 
		printf("\tCandidate estimate = %lf\n", objEst);
#endif

		/* 3. Add the affine function to the master problem. */
		/* 3a. Add cut to the set */
        cut->alphaIncumb =  cut->alpha - vXv(cut->beta , cell->incumbX, index, prob[0]->num->cols) ;		/* right-hand side when using QP master, this is useful for quick updates */
		cell->cuts->vals[cell->cuts->cnt] = cut;

		cell->cuts->vals[cell->cuts->cnt]->rowNum = prob[0]->num->rows + cell->cuts->cnt;
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
		StartMas = clock();
		if ( solveProblem(cell->master->model) ) {
			errMsg("solver", "fullSolve", "failed to solve the master problem", 0);
			return 1;
		}
		EndMas = clock();
		cell->Tmas = cell->Tmas + (EndMas - StartMas);

#if defined(ALGO_CHECK)
		printf("\tObjective function value = %lf\n", getObjective(cell->master->model));
#endif
		Tempobj = cell->obj;
		printf("\t%d: Objective function value = %lf\n", cell->numit, getObjective(cell->master->model) + incumbObj);
		cell->obj = getObjective(cell->master->model) + incumbObj;
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
	if( addRow(model, lenX+1, cut->alphaIncumb, GE, rmatind, cut->beta, cut->name) ) {
		errMsg("solver", "addCut2Solver", "failed to addrow", 0);
		return 1;
	}


//	writeProblem(model, "addcut.lp");

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
