#include "stochasticQP.h"
oneCut* dualSolve(probType** prob, cellType* cell, stocType* stoch, double* x, double solveset) {

	/*initialization of the parameters*/
	sparseMatrix* COmega; /* Presenting the C matrix associated with an observation(I mean the difference from Cbar)*/
	sparseVector* bOmega;  /* Presenting the b vector associated with an observation(I mean the difference from bBar)*/
	sparseVector* ybar; /* Presenting the upperbound  vector associated with an observation(I mean the difference from yBar)*/
	sparseVector* yund; /* Presenting the lowerbound  vector associated with an observation(I mean the difference from mean of yunderscore)*/
	int* omegaP;
	int* omegaQ;
	bOmega = (sparseVector*)mem_malloc(sizeof(sparseVector));
	COmega = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));
	ybar = (sparseVector*)mem_malloc(sizeof(sparseVector));
	yund = (sparseVector*)mem_malloc(sizeof(sparseVector));

	DualType* dual = buildDual(prob[1]);
	buildbcOmega(bOmega, COmega, prob[1], yund, ybar);

	omegaP = (int*)arr_alloc(solveset, int);
	omegaQ = (int*)arr_alloc(2 * cell->omega->cnt - solveset, int);
	SampleGen(omegaP, omegaQ, cell->omega->cnt, solveset);


	double bestalpha;
	int* index = (int*)arr_alloc(prob[1]->num->prevCols + 1, int);
	double maxobj = -100;
	int bestindex = 0;
	for (int i = 1; i <= prob[1]->num->prevCols; i++) {
		index[i] = i;
	}

	/* 1. define a new cut */
	oneCut* cut = newCut(prob[1]->num->cols);

	double tempobj=0; /* holds the best dual objective function value when we are searching over lambda for a suboptimal solution*/
	double mubBar; /* sum over dual*bounds */


	int stat=0; /* if a dual exists in lambda or not*/
	double alpha;


	/* 2. loop through subset solveset and solve the subproblems */
	for (int i = 0; i < solveset; i++) {

		/* 2a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */

		if (solveSubprob(cell, prob[1], cell->subprob, cell->candidX, cell->omega->vals[omegaP[i]], dual->pi, &mubBar, dual->umu
			, dual->lmu)) {
			errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			goto TERMINATE;
		}

		/*2b. update sigma lambda delta*/
		double* fbeta = (double*)arr_alloc(prob[1]->num->prevCols + 1, double);
		stocUpdateQP(cell, prob[1], dual, &alpha, fbeta, mubBar, COmega, bOmega, ybar, yund, omegaQ[i]);

		/*2c. Aggregate the cut coefficients by weighting by observation probability. */
		cut->alpha += cell->omega->probs[omegaQ[i]] * alpha;
		for (int c = 1; c <= prob[1]->num->prevCols; c++) {
			cut->beta[c] += cell->omega->probs[omegaQ[i]] * fbeta[c];
		}
		mem_free(fbeta);
	}



	

	for (int i = 0; i < cell->omega->cnt - solveset; i++) {

		/* 3a. loop through remaining observations and complete delta */

		for (int j = 0; j < cell->lambda->cnt; j++) {
			AddtoDel(cell, prob[1], COmega, bOmega, ybar, yund, omegaQ[i], j);
		}

		double* bestbeta = (double*)arr_alloc(prob[1]->num->prevCols + 1, double);
		maxobj = -1000;
		bestindex = -1;
		
		for (int j = 0; j < cell->lambda->cnt; j++) {
			tempobj = 0;
			alpha = cell->sigma->vals[j]->interceptBar + cell->delta->vals[j][omegaQ[i]]->dalpha;
			tempobj = alpha + vXv(cell->sigma->vals[j]->piCar, cell->candidX, index, prob[1]->num->prevCols) + vXvSparse(cell->candidX, cell->delta->vals[j][i]->dbeta);

		/*3b. calculate estimated Obj value beta x+alpha*/

			if (tempobj > maxobj) {
				maxobj = tempobj;
				bestindex = j;
				bestalpha = alpha;
				
			}
		}
		VVsparse(bestbeta, cell->sigma->vals[bestindex]->piCar, cell->delta->vals[bestindex][omegaQ[i]]->dbeta, prob[1]->num->prevCols + 1);
		/* 3c. Aggregate the cut coefficients by weighting by observation probability. */

		cut->alpha += cell->omega->probs[omegaQ[i]] * bestalpha;
		for (int c = 1; c <= prob[1]->num->prevCols; c++) {
			cut->beta[c] += cell->omega->probs[omegaQ[i]] * cell->sigma->vals[bestindex]->piCar[c];
		}
		mem_free(bestbeta);
	}	
		mem_free(omegaP);
		mem_free(omegaQ);
		return cut;
	TERMINATE:
		return NULL;
}




