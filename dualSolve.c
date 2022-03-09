#include "stochasticQP.h"

oneCut* dualSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset) {
	bool *omegaP;
	int  lambdaIdx;
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

			/*3b. update sigma, lambda, and delta structures */
			lambdaIdx = stocUpdateQP(cell, prob, dual, COmega, bOmega, uOmega, lOmega);

			/*3c. Calculate observations specific coefficients. */
			double* beta = (double*)arr_alloc(prob->num->prevCols + 1, double);
			alpha = cell->sigma->vals[lambdaIdx]->alpha + cell->delta->vals[lambdaIdx][obs]->alpha;
			for (int c = 1; c <= prob->num->cntCcols; c++)
				beta[prob->coord->CCols[c]] += cell->sigma->vals[lambdaIdx]->beta[c];
			for (int c = 1; c <= prob->num->rvCOmCnt; c++)
				beta[prob->coord->rvCols[c-1]] += cell->delta->vals[lambdaIdx][obs]->beta[c];

#if defined(STOCH_CHECK)
			printf("Reconstructed objective function    = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
#endif

			/*3d. Aggregate the cut coefficients by weighting by observation probability. */
			cut->alpha += cell->omega->probs[obs] * alpha;
			for (int c = 1; c <= prob->num->prevCols; c++) {
				cut->beta[c] += cell->omega->probs[obs] * beta[c];
			}
			mem_free(beta);
		}
	}

	/* 4. loop through subset omegaP and use argmax on subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		if ( !omegaP[obs] ) {
			/* 4a. Identify the best dual using the argmax operation */
			lambdaIdx = argmax(prob, cell->sigma, cell->delta, cell->candidX, obs);

			/* 4b. Calculate observations specific coefficients. */
			double* beta = (double*) arr_alloc(prob->num->prevCols + 1, double);
			alpha = cell->sigma->vals[lambdaIdx]->alpha + cell->delta->vals[lambdaIdx][obs]->alpha;
			for (int c = 1; c <= prob->num->cntCcols; c++)
				beta[prob->coord->CCols[c]] += cell->sigma->vals[lambdaIdx]->beta[c];
			for (int c = 1; c <= prob->num->rvCOmCnt; c++)
				beta[prob->coord->rvCols[c - 1]] += cell->delta->vals[lambdaIdx][obs]->beta[c];

#if defined(STOCH_CHECK)
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

			//printf("Reconstructed objective function apppro    = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
#endif

			/* 4c. Aggregate the cut coefficients by weighting by observation probability. */
			cut->alpha += cell->omega->probs[obs] * alpha;

			for (int c = 1; c <= prob->num->prevCols; c++) {
				cut->beta[c] += cell->omega->probs[obs] * beta[c];
			}

			mem_free(beta);
		}
	}

	mem_free(omegaP);

	return cut;
	TERMINATE:
	return NULL;
}//END dualSolve()

int argmax(probType *prob, sigmaType *sigma, deltaType *delta, dVector Xvect, int obs) {
	double maxobj = -INFINITY;
	int lambdaIdx = -1;

	for (int j = 0; j < sigma->cnt; j++) {
		double tempobj;

		tempobj = sigma->vals[j]->alpha + delta->vals[j][obs]->alpha
				- vXv(Xvect, sigma->vals[j]->beta, prob->coord->CCols, prob->num->cntCcols)
				- vXv(Xvect, delta->vals[j][obs]->beta, prob->coord->rvCOmCols, prob->num->rvCOmCnt);

		/*4b. calculate estimated Obj value beta x+alpha*/
		if (tempobj > maxobj) {
			maxobj = tempobj;
			lambdaIdx = j;
		}
	}

	return lambdaIdx;
}//END argmax()




