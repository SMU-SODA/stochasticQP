#include "stochasticQP.h"

oneCut* dualSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset) {
	bool *omegaP;
	/* Elements for argmax calculations */
	double 	bestalpha;
	int* 	index = (int*) arr_alloc(prob->num->prevCols + 1, int);
	double 	maxobj = -100;
	int 	bestindex = 0;
	for (int i = 1; i <= prob->num->prevCols; i++) {
		index[i] = i;
	}
	double tempobj=0; /* holds the best dual objective function value when we are searching over lambda for a suboptimal solution*/

	double mubBar; /* sum over dual*bounds */
	int stat=0; /* if a dual exists in lambda or not*/
	double alpha;

	/* initialization of the parameters */
	sparseVector* bOmega;  	/* Presenting the b vector associated with an observation(I mean the difference from bBar)*/
	sparseMatrix* COmega; 	/* Presenting the C matrix associated with an observation(I mean the difference from Cbar)*/
	sparseVector* dOmega;	/* Presenting the cost coefficient vector associated with an observation */
	sparseVector* uOmega;	/* Presenting the upperbound  vector associated with an observation(I mean the difference from yBar)*/
	sparseVector* lOmega;	/* Presenting the lowerbound  vector associated with an observation(I mean the difference from mean of yunderscore)*/
	bOmega = (sparseVector*)mem_malloc(sizeof(sparseVector));
	dOmega = (sparseVector*)mem_malloc(sizeof(sparseVector));
	uOmega = (sparseVector*)mem_malloc(sizeof(sparseVector));
	lOmega = (sparseVector*)mem_malloc(sizeof(sparseVector));
	COmega = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));

	bOmega->cnt = prob->num->rvbOmCnt;
	bOmega->col = prob->coord->rvbOmRows;

	COmega->cnt = prob->num->rvCOmCnt;
	COmega->col = prob->coord->rvCOmCols;
	COmega->row = prob->coord->rvCOmRows;

	dOmega->cnt = prob->num->rvdOmCnt;
	dOmega->col = prob->coord->rvdOmCols;

	uOmega->cnt = prob->num->rvyuOmCnt;
	uOmega->col = prob->coord->rvyuOmRows;

	lOmega->cnt = prob->num->rvylOmCnt;
	lOmega->col = prob->coord->rvylOmRows;

	/* Structure to hold dual solutions */
	solnType * dual = buildDual(prob->num);

	/* 1. define a new cut */
	oneCut* cut = newCut(prob->num->cols);

	/* 2. Generate a subset */
	omegaP = subsetGenerator(cell->omega->cnt);

	/* 3. loop through subset omegaP and solve the subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		if ( omegaP[obs] ) {
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0];
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1];
			dOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[2];
			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3];
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4];

			/* 3a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
			if (solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], bOmega, COmega, dOmega, lOmega, uOmega, dual)) {
				errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
				goto TERMINATE;
			}

			/*3b. update sigma lambda delta*/
			double* fbeta = (double*)arr_alloc(prob->num->prevCols + 1, double);
			stocUpdateQP(cell, prob, dual, &alpha, fbeta, COmega, bOmega, uOmega, lOmega, obs);

			/*3c. Aggregate the cut coefficients by weighting by observation probability. */
			cut->alpha += cell->omega->probs[obs] * alpha;
			for (int c = 1; c <= prob->num->prevCols; c++) {
				cut->beta[c] += cell->omega->probs[obs] * fbeta[c];
			}
			mem_free(fbeta);
		}
	}

	///* 4. loop through subset omegaP and use argmax on subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		if ( !omegaP[obs] ) {
			/* 4a. loop through remaining observations and complete delta */
			for (int j = 0; j < cell->lambda->cnt; j++) {
				AddtoDel(cell, prob, COmega, bOmega, uOmega, lOmega, obs, j);
			}

			double* bestbeta = (double*)arr_alloc(prob->num->prevCols + 1, double);
			maxobj = -10000000;
			bestindex = -1;

			for (int j = 0; j < cell->lambda->cnt; j++) {
				tempobj = 0;
				alpha = cell->sigma->vals[j]->interceptBar + cell->delta->vals[j][obs]->dalpha;
				tempobj = alpha + vXv(cell->sigma->vals[j]->piCar, cell->candidX, index, prob->num->prevCols) +
				vXvSparse(cell->candidX, cell->delta->vals[j][obs]->dbeta);

	//			/*4b. calculate estimated Obj value beta x+alpha*/
				if (tempobj > maxobj) {
					maxobj = tempobj;
					bestindex = j;
					bestalpha = alpha;
				}
			}
			VsumVsparse(bestbeta, cell->sigma->vals[bestindex]->piCar, cell->delta->vals[bestindex][obs]->dbeta, prob->num->prevCols + 1);

	//		/* 4c. Aggregate the cut coefficients by weighting by observation probability. */
			cut->alpha += cell->omega->probs[obs] * bestalpha;
			for (int c = 1; c <= prob->num->prevCols; c++) {
				cut->beta[c] += cell->omega->probs[obs] * cell->sigma->vals[bestindex]->piCar[c];
			}
			mem_free(bestbeta);
		}
	}

	mem_free(omegaP);

	return cut;
	TERMINATE:
	return NULL;
}//END dualSolve()




