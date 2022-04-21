#include "stochasticQP.h"

int partSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset) {
	double alpha;
	int lambdaIdx;

	/* initialization of the parameters */
	int up = 0, inact = 0, low = 0; /*Number of variables on their bounds*/
	bool newPartFlag = false;
	Mat W;
	Mat T;
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
	solnType* soln = buildSolnType(prob->num);

	bool* omegaP;


	/* 1. define a new cut */
	oneCut* cut = newCut(prob->num->cols);

	/* 2. Generate a subset */
	omegaP = subsetGenerator(cell->omega->cnt);

	/*3. Initialize the partition vector*/

	int partIndx = 0;
	/* 4. loop through subset omegaP and solve the subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		if (omegaP[obs]) {
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0] - 1;
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1] - 1;
			dOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[2] - 1;
			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3] - 1;
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4] - 1;

			/* 4a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */

			if (solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], bOmega, COmega, dOmega, lOmega, uOmega, soln)) {
				errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
				goto TERMINATE;
			}

			/* 4b. Calculate the partition */

			newPartFlag = false;
			partIndx = AddtoPart(prob, cell, uOmega, lOmega, soln, &newPartFlag, &up, &inact, &low);

			/* 4d. Store the fixed parts of current partition if needed*/

			if (newPartFlag) {

				/* 4d.1 Extract the WT matrices*/

				CalcWT(cell, prob, prob->sp->objQ, prob->Dbar, &W, &T);

				/* 4d.2 Add the obtained solution to the lambda structure*/

				addtoLambdaP(cell, soln, &W, &T, prob, uOmega, dOmega, inact, up);


				/* 4d.2 Add  to alpha and beta Bar*/



				/* 4d.3 add to  delta sol and complete a row*/

				addtoDeltaSol(cell, soln, &W, &T, prob, uOmega, dOmega, inact, up);

				/* 4d.4 add to alpha and beta delta and complete a row*/

				AddtoDeltaP(cell, soln, prob, bOmega, uOmega, lOmega);


				/*3c. Calculate observations specific coefficients. */
				double* beta = (double*)arr_alloc(prob->num->prevCols + 1, double);

				alpha = cell->sigma->vals[partIndx]->alpha + cell->delta->vals[partIndx][obs]->alpha;

				for (int c = 1; c <= prob->num->cntCcols; c++)
					beta[prob->coord->CCols[c]] += cell->sigma->vals[partIndx]->beta[c];
				for (int c = 1; c <= prob->num->rvCOmCnt; c++)
					beta[prob->coord->rvCOmCols[c]] += cell->delta->vals[lambdaIdx][obs]->beta[c];


				double obj;
				obj = getObjective(cell->subprob->model);
#if defined(STOCH_CHECK)
			printf("Reconstructed objective function (exact)  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
			printf("Dif  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols)-obj );
#endif

			/*3d. Aggregate the cut coefficients by weighting by observation probability. */
				cut->alpha += cell->omega->probs[obs] * alpha;
				for (int c = 1; c <= prob->num->prevCols; c++) {
					cut->beta[c] += cell->omega->probs[obs] * beta[c];
				}
				mem_free(beta);
			}
		}
	}
	return cut;

TERMINATE:
	return NULL;
}
	