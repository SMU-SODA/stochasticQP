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
	oneCut *cut = NULL;
	clock_t StartCut;
	clock_t EndCut;
	clock_t StartMas;
	clock_t EndMas;
	double subset = config.SAMPLE_FRACTION * cell->omega->cnt; /*initialize the number of samples you want to take from Omega in each iteration*/
	clock_t tStart = clock();
	double Tempobj = -100000;
	cell->obj = -50000;

	cell->k  = 0;
	while (cell->obj - Tempobj > 0.001 ||  cell->numit < 50){
		cell->numit++;

		/* 1. Check optimality */

		/* 2. Switch between algorithms to add a new affine functions. */
		switch (config.ALGOTYPE) {
		case 0:
			StartCut = clock();
			cut = fullSolveCut(prob[1], cell, stoc, cell->candidX);
			EndCut = clock();
			cell->Tcut = cell->Tcut + (EndCut - StartCut);
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
			StartCut = clock();
			cut = partSolve(prob[1],  cell,  stoc, cell->candidX, subset);			
			EndCut = clock();
			cell->Tcut = cell->Tcut + (EndCut - StartCut);
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
		cell->cuts->vals[cell->cuts->cnt] = cut;
		cell->cuts->cnt++;




		/* 3b. Update the master problem by adding the cut*/
		if ( addCut2Solver(cell->master->model, cut, prob[0]->num->cols) ) {
			errMsg("solver", "fullSolve", "failed to add cut to the master problem", 0);
			return 1;
		}


		/*update the master for ncumbent*/

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
		printf("\t%d: Objective function value = %lf\n", cell->numit, getObjective(cell->master->model));
		cell->obj = getObjective(cell->master->model);
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
	if( addRow(model, lenX+1, cut->alpha, GE, rmatind, cut->beta, cut->name) ) {
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

