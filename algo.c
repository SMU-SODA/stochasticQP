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

int runAlgo (probType **prob, stocType *stoc, cellType *cell) {
	oneCut *cut;

	while ( cell->k < 10) {
		cell->k++;
		/* 1. Check optimality */

		/* 2. Switch between algorithms to add a new affine functions. */
		switch (config.ALGOTYPE) {
		case 0:
			cut = fullSolve(prob, cell, stoc, cell->candidX);
			break;
		case 1:
			partSolve();
			break;
		case 2:
			dualSolve();
			break;

		default:
			errMsg("ALGO", "main", "Unknown algorithm type", 0);

			goto TERMINATE;
		}

		/* 3. Add the affine function to the master problem. */
		/* 3a. Add cut to the set */
		cell->cuts->vals[cell->cuts->cnt] = cut;
		cell->cuts->cnt++;

		/* 3b. Update the master problem by adding the cut*/


		if ( addCut2Solver(cell->master, cut) ) {
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

	///*invoke the algorithm*/

	/* To be edited */
	//	double x = 0.1;
	//	for (int i = 0; i < omega->cnt; i++) {
	//		if (randUniform() < x) {
	//
	//		}
	//	}


	return 0;

	TERMINATE:
	return 1;
}//END runALgo()

int addCut2Solver(oneProblem *master, oneCut *cut) {
	dVector matvals;
	iVector rmatind;
	rmatind = (iVector) arr_alloc(master->mac, int);
	matvals = (dVector) arr_alloc(master->mac ,double);
	int j = 0;
	matvals[j] = 1;
	rmatind[j] = master->mac - 1;
	j++;
	for (int i = 0; i < master->mac-1 ;++i) {
		if (cut->beta[i] != 0) {
			matvals[j] = cut->beta[i];
			rmatind[j] = i;
			j++;
		}
	}

	static int cummCutNum = 0;
	/* Set up the cut name */
	sprintf(cut->name, "cut_%04d", cummCutNum++);

	/* Add a new linear constraint to a model. */
	if(addRow(master->model, j, cut->alpha, GE, rmatind, matvals, cut->name)) {
		errMsg("solver", "addCut2Solver", "failed to addrow", 0);
		return 1;
	}

	mem_free(rmatind);
	mem_free(matvals);

	return 0;
}//END addCuts2Cell()

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

