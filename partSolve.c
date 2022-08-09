#include "stochasticQP.h"

oneCut *partSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset) {

	double alpha;
	int lambdaIdx;

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

	bool* omegaP;

	/* 1. define a new cut */
	oneCut* cut = newCut(prob->num->cols);

	/* 2. Generate a subset */
	omegaP = subsetGenerator(cell->omega->cnt);

	/*3. Initialize the partition vector*/
	int partIndx = 0;
	
	/* Build a basis */
	long long int* basis = (long long int*) arr_alloc(prob->num->cols + 1,long long int);
	Buildbase(basis, prob->num->cols, 3);

	/* 4. loop through subset omegaP and solve the subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {

		if (omegaP[obs]) {
			solnType* soln = buildSolnType(prob->num);
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
			_CrtDumpMemoryLeaks();
			StocUpdatePart( cell,  prob,  bOmega,  COmega, lOmega,  uOmega,  soln, basis , &partIndx);

			//			/*4b. Calculate observations specific coefficients. */
			//			double* beta = (double*)arr_alloc(prob->num->prevCols + 1, double);
			//
			//			alpha = cell->sigma->vals[partIndx]->alpha + cell->delta->vals[partIndx][obs]->alpha;
			//
			//			for (int c = 1; c <= prob->num->cntCcols; c++)
			//				beta[prob->coord->CCols[c]] += cell->sigma->vals[partIndx]->beta[c];
			//			for (int c = 1; c <= prob->num->rvCOmCnt; c++)
			//				beta[prob->coord->rvCOmCols[c]] += cell->delta->vals[lambdaIdx][obs]->beta[c];
			//
			//
			//			double obj;
			//			obj = getObjective(cell->subprob->model);
			//#if defined(STOCH_CHECK)
			//			printf("Reconstructed objective function (exact)  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
			//			printf("Dif  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols)-obj );
			//#endif
			//
			//			/* 4c. Aggregate the cut coefficients by weighting by observation probability. */
			//			cut->alpha += cell->omega->probs[obs] * alpha;
			//			for (int c = 1; c <= prob->num->prevCols; c++) {
			//				cut->beta[c] += cell->omega->probs[obs] * beta[c];
			//			}
			//			mem_free(beta);
			freeSolnType(soln);
		}
	}

	//	/* loop through the rest of the observations */
	//	/* 4. loop through subset omegaP and use argmax on subproblems */
	//	for (int obs = 0; obs < cell->omega->cnt; obs++) {
	//		if (!omegaP[obs]) {
	//			/* 4a. Identify the best dual using the argmax operation */
	//			lambdaIdx = argmax(prob, cell->sigma, cell->delta, cell->candidX, obs);
	//
	//			/* 4b. Calculate observations specific coefficients. */
	//			double* beta = (double*)arr_alloc(prob->num->prevCols + 1, double);
	//			alpha = cell->sigma->vals[lambdaIdx]->alpha + cell->delta->vals[lambdaIdx][obs]->alpha;
	//			for (int c = 1; c <= prob->num->cntCcols; c++)
	//				beta[prob->coord->CCols[c]] += cell->sigma->vals[lambdaIdx]->beta[c];
	//			for (int c = 1; c <= prob->num->rvCOmCnt; c++)
	//				beta[prob->coord->rvCOmCols[c]] += cell->delta->vals[lambdaIdx][obs]->beta[c];
	//
	//#if defined(STOCH_CHECK)
	//			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0] - 1;
	//			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1] - 1;
	//			dOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[2] - 1;
	//			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3] - 1;
	//			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4] - 1;
	//
	//			/* 3a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
	//			if (solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], bOmega, COmega, dOmega, lOmega, uOmega, soln)) {
	//				errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
	//				goto TERMINATE;
	//			}
	//
	//			printf("Reconstructed objective function (approx) = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
	//#endif
	//
	//			/* 4c. Aggregate the cut coefficients by weighting by observation probability. */
	//			cut->alpha += cell->omega->probs[obs] * alpha;
	//			for (int c = 1; c <= prob->num->prevCols; c++) {
	//				cut->beta[c] += cell->omega->probs[obs] * beta[c];
	//			}
	//			mem_free(beta);
	//		}
	//	}

	mem_free(omegaP);
	mem_free(bOmega);
	mem_free(dOmega);
	mem_free(uOmega);
	mem_free(lOmega);
	mem_free(COmega);
	mem_free(basis);
	
	return cut;

	TERMINATE:
	return NULL;
}//END partSolve()

int StocUpdatePart(cellType* cell, probType* prob, sparseVector* bOmega, sparseMatrix* COmega , sparseVector* lOmega,
		sparseVector* uOmega, solnType* soln , long long int* basis , int* partIndx) {
	int up = 0, inact = 0, low = 0; /*Number of variables on their bounds*/
	bool newPartFlag = false;
	Mat* W ;
	Mat* T ;

	/* 4b. Calculate the partition */
	(*partIndx) = addtoPartition(prob, cell, uOmega, lOmega, soln, &newPartFlag, &up, &inact, &low, basis);
	_CrtDumpMemoryLeaks();
	/* 4d. Store the fixed parts of current partition if needed*/
	if (newPartFlag) {

		/* 4d.1 Extract the WT matrices*/
		CalcWT(cell, prob, prob->sp->objQ, prob->Dbar, &W, &T, low, up, inact);

		Mat* WT = CombineWT(prob, W, T, low, up, inact);

		/* 4d.2 Add the obtained solution to the lambda structure*/
		addtoLambdaP(cell, soln, WT, prob, bOmega, uOmega, lOmega, low, up, inact);

		/* 4d.2 Add to alpha and beta Bar*/
		AddtoSigmaP(cell, soln, prob);

		/* 4d.3 add to  delta sol and complete a row*/
		for (int obs = 0; obs < cell->omega->cnt; obs++) {
			/* Add a new row to the delta structure for all observations and the latest lambda (lambdaIdx) */
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0] - 1;
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1] - 1;
			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3] - 1;
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4] - 1;
			addtoDeltaP(cell, soln, W, T, WT, prob, COmega, bOmega, uOmega, lOmega, obs, (*partIndx), inact, up, low);
		}
		freemat(WT);
		freemat(W);
		freemat(T);
	}

	return 0;
}; //EndStocUpdatePart
