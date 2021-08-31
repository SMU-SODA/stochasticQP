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

int runAlgo (probType *prob, stocType *stoc, cellType *cell) {

	while ( cell->k < config.MAX_ITER ) {

		/* 1. Check optimality */

		/* 2. Switch between algorithms to add a new affine functions. */
		switch (config.ALGOTYPE) {
		case 0:
			fullSolve(cell, cell->omega, stoc, prob[1]->Cbar, cell->candidX, 0);
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

		/* 4. Solve the master problem. */


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

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta dVector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX, int numIstar, int numSamples) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->name = (cString) arr_alloc(NAMESIZE, char);
	cut->rowNum = -1;

	if ( numIstar > 0 ) {
		if (!(cut->iStar = arr_alloc(numIstar, int)))
			errMsg("allocation", "new_cut", "iStar", 0);
	}
	else
		cut->iStar = NULL;


	cut->beta = arr_alloc(numX + 1, double);
	cut->beta[0] = 1.0;
	cut->alpha = 0.0;
	cut->alphaIncumb = 0.0;

	return cut;
}//END newCut

