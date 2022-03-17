#include "stochasticQP.h"

extern configType config;

bool *subsetGenerator(int numObs) {
	bool* omegaP = (bool*)arr_alloc(numObs, bool);
	int cnt = 0;

	Repeat:
	for (int obs = 0; obs < numObs; obs++) {
		omegaP[obs] = randUniform() <= config.SAMPLE_FRACTION;
		if(omegaP[obs] == true){
			cnt++;
		}
	}
	if (cnt == 0) {
		goto Repeat;
	}
	return omegaP;
}//END subsetGenerator()

void  buildOmegaCoordinates (probType *prob, sparseVector bOmega, sparseMatrix COmega, sparseVector dOmega, sparseVector uOmega, sparseVector lOmega) {

	bOmega.cnt = prob->num->rvbOmCnt;
	bOmega.col = prob->coord->rvbOmRows;

	COmega.cnt = prob->num->rvCOmCnt;
	COmega.col = prob->coord->rvCOmCols;
	COmega.row = prob->coord->rvCOmRows;

	dOmega.cnt = prob->num->rvdOmCnt;
	dOmega.col = prob->coord->rvdOmCols;

	uOmega.cnt = prob->num->rvyuOmCnt;
	uOmega.col = prob->coord->rvyuOmRows;

	lOmega.cnt = prob->num->rvylOmCnt;
	lOmega.col = prob->coord->rvylOmRows;
}//END buildOmegaCoordinates()

void addtoSigma(cellType* cell, probType* prob, solnType *soln) {
	int* index;
	int obs = cell->lambda->cnt-1;
	double x;
	index = (int*)arr_alloc(prob->num->cols + 1,int);
	for (int i = 1; i <= prob->num->cols; i++) {
		index[i] = i ;
	}

	cell->sigma->vals[obs] = (pixbCType *) mem_malloc(sizeof(pixbCType));

	/*Calculate fixed section of beta = Cbar*pi (pi is the dual vector associated with equality constraints)*/
	dVector fbeta = vxMSparse(soln->pi, prob->Cbar, prob->num->prevCols);
	cell->sigma->vals[obs]->beta = reduceVector(fbeta, prob->coord->CCols, prob->num->cntCcols);

	/*Calculate fixred section of alpha -.5 yQy  + bbar* pi + yunder nu - ybar mu)*/
	cell->sigma->vals[obs]->alpha = vXvSparse(soln->pi, prob->bBar) + vXvSparse(soln->lmu, prob->lBar) + vXvSparse(soln->umu, prob->uBar);
	if ( prob->sp->objQ != NULL) {
		dVector yTopQbar = vxMSparse(soln->y, prob->sp->objQ, prob->num->cols);
		x = 0.5 * vXv(yTopQbar, soln->y, index, prob->num->cols);
		cell->sigma->vals[obs]->alpha -= vXv(yTopQbar, soln->y, index, prob->num->cols);
		mem_free(yTopQbar);
	}
	cell->sigma->cnt++;

	mem_free(fbeta);
	mem_free(index);

}/* End addtoSigma() */

void addtoDelta(cellType* cell, probType* prob, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* uOmega, sparseVector* lOmega, int obs ,int numPi) {

	double* dbeta = vxMSparse(cell->lambda->pi[numPi], COmega, prob->num->prevCols);

	/* If this is a new lambda, we add a row to the delta structure */
	if ( obs ==  0 ) {
		cell->delta->vals[numPi] = (pixbCType **) arr_alloc(cell->omega->cnt, pixbCType *);
		cell->delta->cnt++;
	}

	cell->delta->vals[numPi][obs] = (pixbCType *) mem_malloc(sizeof(pixbCType));

	/* calculate alpha and beta*/
	cell->delta->vals[numPi][obs]->alpha = vXvSparse(cell->lambda->pi[numPi], bOmega)
														+ vXvSparse(cell->lambda->umu[numPi], uOmega) + vXvSparse(cell->lambda->lmu[numPi], lOmega); /*TO DO: ybar and yund vals start from index 0*/

	if ( prob->num->rvCOmCnt > 0 )
		cell->delta->vals[numPi][obs]->beta = reduceVector(dbeta, prob->coord->rvCOmCols, prob->num->rvCOmCnt);
	else
		cell->delta->vals[numPi][obs]->beta = NULL;

	free(dbeta);

}//END addtoDelta()

int addtoLambda(lambdaType* lambda, solnType *dual, int numRows, int numCols, bool *newLambdaFlag) {

	/* Loop to see if the current lambda already exists */
	int idx = 0;
	while (idx < lambda->cnt) {
		if (dual->pi[0] == lambda->pi[idx][0] && dual->umu[0] == lambda->umu[idx][0] && dual->lmu[0] == lambda->lmu[idx][0]) {
			break;
		}
		idx++;
	}

	if ( idx == lambda->cnt ) {
		/* TODO: New lambda discovered */
		(*newLambdaFlag) = true;
		lambda->pi[lambda->cnt]  = duplicVector(dual->pi, numRows+1);
		lambda->umu[lambda->cnt] = duplicVector(dual->umu, numCols+1);
		lambda->lmu[lambda->cnt] = duplicVector(dual->lmu, numCols+1);
		lambda->cnt++;
	}

	return idx;
}//END addtoLambda()



/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a dVector to hold an array of
 * observation and the probability associated with it. */
omegaType* newOmega(stocType* stoc) {
	omegaType* omega;
	int cnt, i, base, idx;

	omega = (omegaType*)mem_malloc(sizeof(omegaType));
	omega->probs = (dVector)arr_alloc(config.MAX_OBS, double);
	omega->weights = (iVector)arr_alloc(config.MAX_OBS, int);
	omega->vals = (dVector*)arr_alloc(config.MAX_OBS, dVector);
	omega->cnt = 0; omega->numRV = stoc->numOmega;

	if (config.SAA == 1 ) { //|| !stoc->isDiscrete
		config.SAA = 1;
		omega->cnt = config.MAX_OBS;
		// return omega;
	}

	if (strstr(stoc->type, "BLOCKS") != NULL) {
		if ((omega->cnt = stoc->numVals[0]) <= config.MAX_OBS) {
			omega->vals = (dVector*)mem_realloc(omega->vals, omega->cnt * sizeof(dVector));
			omega->probs = (dVector)mem_realloc(omega->probs, omega->cnt * sizeof(double));
			for (cnt = 0; cnt < omega->cnt; cnt++) {
				omega->probs[cnt] = stoc->probs[0][cnt];
				if (!(omega->vals[cnt] = (dVector)arr_alloc(omega->numRV + 1, double)))
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				for (i = 0; i < omega->numRV; i++)
					omega->vals[cnt][i + 1] = stoc->vals[i][cnt] - stoc->mean[i];
				omega->vals[cnt][0] = oneNorm(omega->vals[cnt] + 1, omega->numRV);
			}
		}
		else {
			omega->cnt = config.MAX_OBS;
			config.SAA = 1;
		}
	}
	else if (strstr(stoc->type, "INDEP") != NULL) {
		omega->cnt = 1; i = 0;
		while (i < stoc->numOmega) {
			omega->cnt *= stoc->numVals[i];
			if (omega->cnt > config.MAX_OBS) {
				omega->cnt = config.MAX_OBS;
				config.SAA = 1;
				break;
			}
			i++;
		}

		if (!config.SAA) {
			omega->vals = (dVector*)mem_realloc(omega->vals, omega->cnt * sizeof(dVector));
			omega->probs = (dVector)mem_realloc(omega->probs, omega->cnt * sizeof(double));
			for (cnt = 0; cnt < omega->cnt; cnt++) {
				if (!(omega->vals[cnt] = (dVector)arr_alloc(omega->numRV + 1, double)))
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				omega->probs[cnt] = 1; base = omega->cnt;
				for (i = 0; i < omega->numRV; i++) {
					base /= stoc->numVals[i];
					idx = (int)((double)cnt / (double)base) % stoc->numVals[i];
					omega->vals[cnt][i ] = stoc->vals[i][idx] - stoc->mean[i];
					omega->probs[cnt] *= stoc->probs[i][idx];
				}
			}
		}
	}
	else {
		omega->cnt = config.MAX_OBS;
		config.SAA = 1;
	}

	return omega;
}//END newOmega()

lambdaType* newLambda(double SigmaSize, probType** prob) {
	lambdaType* lambda = NULL;
	/* Assign memory to lambda structure which is the set of dual solutions we obtained so far. pi is related to equality constraints, mu2
corresponds to upper bounds and mu3 corresponds to lower bounds */

	lambda = (lambdaType*)mem_malloc(sizeof(lambdaType));
	lambda->pi = (double**)arr_alloc(SigmaSize, double*);
	lambda->umu = (double**)arr_alloc(SigmaSize, double*);
	lambda->lmu = (double**)arr_alloc(SigmaSize, double*);
	lambda->cnt = 0;

	return lambda;
}//END newLambda()

sigmaType* newSigma(double SigmaSize, probType** prob ) {
	sigmaType* sigma = NULL; /* Sigma is a collection of fixed parts of alpha and beta which is independent of the observation */

	sigma = (sigmaType*)mem_malloc(sizeof(sigmaType));
	sigma->vals = (pixbCType**) arr_alloc(SigmaSize, pixbCType*);
	sigma->cnt = 0;

	return sigma;
}//END newSigma()

deltaType* newDelta(double SigmaSize, probType** prob , cellType* cell) {
	deltaType* delta = NULL;
	/* assign memory to deta structure, this will record the deltaAlpha and deltaBetha associated with each observation*/

	delta = (deltaType*)mem_malloc(sizeof(deltaType));
	delta->vals = (pixbCType ***)arr_alloc(SigmaSize, pixbCType **);
	delta->cnt = 0;

	return delta;

}//END newDelta()

void freeSigma(sigmaType* sigma) {
	if (sigma) {
		if (sigma->vals) {
			for (int i = 0; i < sigma->cnt; i++) {

				if (sigma->vals[i])
				{
					if (sigma->vals[i]->beta) {
						mem_free(sigma->vals[i]->beta);
					}


				}
				mem_free(sigma->vals[i]);
			}
			mem_free(sigma->vals);
		}
		mem_free(sigma);
	}

}; /**EndFreeSigma**/

void freeLambda(lambdaType* lambda) {

	if (lambda)
	{
		if (lambda->pi)
		{
			for (int i = 0; i < lambda->cnt;i++) {
				mem_free(lambda->pi[i]);
			}
			mem_free(lambda->pi);
		}
		if (lambda->umu)
		{
			for (int i = 0; i < lambda->cnt; i++) {
				mem_free(lambda->umu[i]);
			}
			mem_free(lambda->umu);
		}
		if (lambda->lmu) {
			for (int i = 0; i < lambda->cnt; i++) {
				mem_free(lambda->lmu[i]);
			}
			mem_free(lambda->lmu);
		}
		mem_free(lambda);
	}

}; /*EndfreeLambda*/

void freeDelta(deltaType* delta, int numobs) {

	if (delta) {
		if (delta->vals) {
			for (int i = 0; i < delta->cnt; i++) {
				if (delta->vals[i]) {
					for (int j = 0; j < numobs; j++) {
						freeLambdaDelta(delta->vals[i][j]);
					}
					mem_free(delta->vals[i]);
				}
				mem_free(delta->vals);
			}
		}
		mem_free(delta);
	}

};/*End freeDelta*/

void freeLambdaDelta(pixbCType* lambdadelta) {
	if (lambdadelta) {
		mem_free(lambdadelta->beta);
	}
}

void freeOmegaType(omegaType* omega, bool partial) {
	int n;
	if (omega->vals) {
		for (n = 0; n < omega->cnt; n++)
			if (omega->vals[n])
				mem_free(omega->vals[n]);
		if (partial) {
			omega->cnt = 0;
			return;
		}
		mem_free(omega->vals);
	}
	if (omega->probs) mem_free(omega->probs);
	if (omega->weights) mem_free(omega->weights);
	mem_free(omega);
}//END freeOmegaType()

//
//int calcSigma(sigmaType* sigma, cellType* cell  ,probType** prob, dVector pi, dVector mu2, dVector mu3 , sparseVector* bOmega, sparseMatrix* COmega,
//		sparseVector* yuOmega , int obs) {
//	double fixedAlpha , alpha1 , lql , alpha, Bx ;
//	/*Check if a the sigma structure should be updated */
//
//	/* put pi.Cbar in sigma (fixed part of beta which does not change when C is deterministic) */
//	sigma->vals[obs]->piCar = vxMSparse(pi, prob[1]->Cbar, prob[1]->num->rows);
//
//
//	/* finding fixed part of alpha equal to -1/2 lambda Q lambda - xiBar .pi - yl.mu3 + yuBar.mu2 */
//
//	/* first  find -1/2 lambda Q lambda = obj -(beta x  - xi(obs) .pi - yl.mu3 + yu(obs).mu2)*/
//	/*- xi(obs) .pi - yl.mu3 + yu(obs).mu2*/
//
//	fixedAlpha = -vXv(prob[1]->bBar, pi + 1, NULL, prob[1]->num->rows)
//														- vXv(prob[1]->lBar, mu3 + 1, NULL, prob[1]->num->cols) + vXv(prob[1]->uBar, mu2 + 1, NULL, prob[1]->num->cols);
//
//
//	alpha1 = fixedAlpha - vXvSparse(pi + 1, bOmega) + vXvSparse(mu2 + 1, yuOmega);
//
//	Bx = vXv(sigma->vals[obs]->piCar, cell->candidX + 1, NULL, prob[0]->num->cols);
//
//	lql = getObjective(cell->subprob->model) - Bx - alpha1;
//
//
//	alpha = fixedAlpha + lql;
//	sigma->vals[obs]->fixed = alpha + Bx;
//	sigma->vals[obs]->interceptBar = alpha;
//	return 0;
//}

