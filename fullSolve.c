#include "stochasticQP.h"

oneCut *fullSolveCut(probType *prob, cellType* cell, stocType* stoch, double* x) {
	sparseMatrix COmega;
	sparseVector bOmega;
	dVector 	 pi, piCBar, beta;
	double 		 alpha, mubBar;

	pi = (dVector) arr_alloc(prob->num->rows+1, double);
	bOmega.cnt = prob->num->rvbOmCnt; bOmega.col = prob->coord->rvbOmRows;
	COmega.cnt = prob->num->rvCOmCnt; COmega.col = prob->coord->rvCOmCols; COmega.row = prob->coord->rvCOmRows;

	/* 1. Create a new cut */
	oneCut *cut = newCut(prob->num->cols);

	/* 2. loop through observations and solve subproblem for all of them. */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {

		/* 2a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
		if ( solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], pi, &mubBar) ) {
			errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			goto TERMINATE;
		}
		cell->LPcnt++;

		/* 2b. Compute the cut coefficients for individual observation. */
		bOmega.val = cell->omega->vals[obs] + prob->coord->rvOffset[0];
		COmega.val = cell->omega->vals[obs] + prob->coord->rvOffset[1];

		/* Optimality cut calculations */
		alpha  = vXvSparse(pi, prob->bBar) + mubBar + vXvSparse(pi, &bOmega);
		beta   = vxMSparse(pi, prob->Cbar, prob->num->prevCols);
		piCBar = vxMSparse(pi, &COmega, prob->num->prevCols);
		for (int c = 1; c <= prob->num->rvCOmCnt; c++)
			beta[prob->coord->rvCOmCols[c]] += piCBar[c];
		mem_free(piCBar);

#if defined(STOCH_CHECK)
		printf("Objective estimate computed as cut height = %lf\n", alpha - vXv(beta, x, NULL, prob->num->prevCols));
#endif

		/* 2c. Aggregate the cut coefficients by weighting by observation probability. */
		cut->alpha += cell->omega->probs[obs]*alpha;
		for (int c = 1; c <= prob->num->prevCols; c++)
			cut->beta[c] += cell->omega->probs[obs]*beta[c];

		mem_free(beta);
	}

	return cut;

	TERMINATE:
	return NULL;
}//END fullSolve()
