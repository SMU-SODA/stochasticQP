#include "stochasticQP.h"

extern configType config;

bool *subsetGenerator(int numObs) {

	bool* omegaP = (bool *) arr_alloc(numObs, bool);
	for (int obs = 0; obs < numObs; obs++) {
		omegaP[obs] = randUniform() <= config.SAMPLE_FRACTION;
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

	index = (int*)arr_alloc(prob->num->cols + 1,int);
	for (int i = 1; i <= prob->num->cols; i++) {
		index[i] = i ;
	}

	/*Calculate fixed section of beta = Cbar*pi (pi is the dual vector associated with equality constraints)*/
	dVector fbeta = vxMSparse(soln->pi, prob->Cbar, prob->num->prevCols);
	copyVector(fbeta, cell->sigma->vals[obs]->piCar, prob->num->prevCols);

	/*Calculate fixred section of alpha -.5 yQy  + bbar* pi + yunder nu - ybar mu)*/
	cell->sigma->vals[obs]->interceptBar = -0.5 * vXv(vxMSparse(soln->y, prob->sp->objQ, prob->num->cols), soln->y, index, prob->num->cols) +
			vXvSparse(soln->pi, prob->bBar) + vXvSparse(soln->lmu, prob->lBar) - vXvSparse(soln->umu, prob->uBar);
			cell->sigma->cnt++;

			free(index);

}/* End addtoSigma() */

void AddtoDel(cellType* cell, probType* prob, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* ybar, sparseVector* yund, int obs ,int num) {

	/* 2b. Compute the cut coefficients for individual observation. */
	int* dbetaindex;
	dbetaindex = (int*)arr_alloc(prob->num->prevCols + 1, int);
	double* dbeta = vxMSparse(cell->lambda->pi[num], COmega, prob->num->prevCols);
	if (cell->delta->vals[num][obs]->state == 0) {
		cell->delta->vals[num][obs]->state = 1;

		/* calculate deltaalpha */
		cell->delta->vals[num][obs]->dalpha = vXvSparse(cell->lambda->pi[num], bOmega) - vXvSparse(cell->lambda->umu[num], ybar) + vXvSparse(cell->lambda->lmu[num], yund); /*TO DO: ybar and yund vals start from index 0*/
		cell->delta->vals[num][obs]->dbeta->cnt = 0;
		/*find out the number of nonzero elements and their location*/
		for (int j = 0; j <= prob->num->prevCols; j++) {
			if (dbeta[j] != 0) {
				cell->delta->vals[num][obs]->dbeta->cnt++;
			}
		}

			/*Assign memory to the index section of deltab */
			cell->delta->vals[num][obs]->dbeta->col = (int*) arr_alloc(cell->delta->vals[num][obs]->dbeta->cnt + 1, int);
			cell->delta->vals[num][obs]->dbeta->val = (double*) arr_alloc(cell->delta->vals[num][obs]->dbeta->cnt + 1, double);
			int cnt = 0;
			for (int j = 0; j < prob->num->prevCols; j++) {
				if (dbeta[j]) {
					cnt++;
					cell->delta->vals[num][obs]->dbeta->val[cnt] = dbeta[j];
					cell->delta->vals[num][obs]->dbeta->col[cnt-1] = j;
				}
			}
		}
	free(dbetaindex);
	free(dbeta);
}

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
		copyVector(dual->pi, lambda->pi[lambda->cnt], numRows);
		copyVector(dual->umu, lambda->umu[lambda->cnt], numCols);
		copyVector(dual->lmu, lambda->lmu[lambda->cnt], numCols);
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
	omega->probs = (dVector) arr_alloc(config.MAX_OBS, double);
	omega->weights = (iVector) arr_alloc(config.MAX_OBS, int);
	omega->vals = (dVector*) arr_alloc(config.MAX_OBS, dVector);
	omega->cnt = 0;
	omega->numRV = stoc->numOmega;

	if (config.SAA == 1 ) {
		config.SAA = 0;
		omega->cnt = config.MAX_OBS;
		return omega;
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
	else if (strstr(stoc->type, "INDEP_DISCRETE") != NULL) {
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

		if (config.SAA) {
			omega->vals = (dVector*)mem_realloc(omega->vals, omega->cnt * sizeof(dVector));
			omega->probs = (dVector)mem_realloc(omega->probs, omega->cnt * sizeof(double));
			for (cnt = 0; cnt < omega->cnt; cnt++) {
				if (!(omega->vals[cnt] = (dVector)arr_alloc(omega->numRV + 1, double)))
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				omega->probs[cnt] = 1; base = omega->cnt;
				for (i = 0; i < omega->numRV; i++) {
					base /= stoc->numVals[i];
					idx = (int)((double)cnt / (double)base) % stoc->numVals[i];
					omega->vals[cnt][i + 1] = stoc->vals[i][idx] - stoc->mean[i];
					omega->probs[cnt] *= stoc->probs[i][idx];
				}}}}

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

	for (int i = 0; i < SigmaSize; i++) {
		lambda->pi[i] = (double*)arr_alloc(prob[1]->num->rows + 1, double);
		lambda->umu[i] = (double*)arr_alloc(prob[1]->num->cols + 1, double);
		lambda->lmu[i] = (double*)arr_alloc(prob[1]->num->cols + 1, double);
	}

	lambda->cnt = 0;
	return lambda;

}//END newLambda()

sigmaType* newSigma(double SigmaSize, probType** prob ) {
	sigmaType* sigma = NULL; /* Sigma is a collection of fixed parts of alpha and beta which is independent of the observation */
	sigma = (sigmaType*)mem_malloc(sizeof(sigmaType));
	sigma->vals = (pixbCType**)arr_alloc(SigmaSize, pixbCType*);
	sigma->cnt = 0;
	for (int i = 0; i < SigmaSize; i++) {
		sigma->vals[i] = (pixbCType*)mem_malloc(sizeof(pixbCType));
	}
	for (int i = 0; i < SigmaSize; i++) {
		sigma->vals[i]->piCar = (double*)arr_alloc(prob[0]->num->cols + 1, double);
	}
	return sigma;

}//END newSigma()

deltaType* newDelta(double SigmaSize, probType** prob , cellType* cell) {
	deltaType* delta = NULL;
	/* assign memory to deta structure, this will record the deltaAlpha and deltaBetha associated with each observation*/

	delta = (deltaType*)mem_malloc(sizeof(deltaType));
	delta->vals = (lambdadeltaType***)arr_alloc(cell->omega->cnt, lambdadeltaType**);
	for (int i = 0; i < cell->omega->cnt; i++) {
		delta->vals[i] = (lambdadeltaType**)arr_alloc(SigmaSize, lambdadeltaType*);
	}
	for (int i = 0; i < cell->omega->cnt; i++) {
		for (int j = 0; j < SigmaSize; j++) {
			delta->vals[i][j] = (lambdadeltaType*)mem_malloc(sizeof(lambdadeltaType));
			delta->vals[i][j]->dbeta = (sparseVector*)mem_malloc(sizeof(sparseVector));
			delta->vals[i][j]->dbeta->val = (double*)arr_alloc(prob[0]->num->cols + 1, double); /*TO DO: Change it to a reasonable size*/
			delta->vals[i][j]->dbeta->col = (int*)arr_alloc(prob[0]->num->cols + 1, int);
			delta->vals[i][j]->state = 0;
		}
	}

	return delta;

}//END newDelta()

void freeSigma(sigmaType* sigma) {
	if (sigma) {
		if (sigma->vals) {
			for (int i = 0; i < sigma->cnt; i++) {

				if (sigma->vals[i])
				{
					if (sigma->vals[i]->piCar) {
						mem_free(sigma->vals[i]->piCar);
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
		if (lambda->lmu)
		{
			for (int i = 0; i < lambda->cnt; i++)
			{
				mem_free(lambda->umu[i]);
			}
			mem_free(lambda->umu);
		}
		mem_free(lambda);
	}
}; /*EndfreeLambda*/

void freeDelta(deltaType* delta, int numobs) {
	if (delta) {
		if (delta->vals) {
			for (int i = 0; i < delta->cnt; i++) {
				if (delta->vals[i]) {
					for (int j = 0; j < numobs; i++) {
						freeLambda(delta->vals[i][j]);
					}
					mem_free(delta->vals[i]);
				}
				mem_free(delta->vals);
			}
		}
		mem_free(delta);
	}
};/*End freeDelta*/

void freeLambdaDelta(lambdadeltaType* lambdadelta){

	if (lambdadelta) {

		freeSparseVector(lambdadelta->dbeta);
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


int calcSigma(sigmaType* sigma, cellType* cell  ,probType** prob, dVector pi, dVector mu2, dVector mu3 , sparseVector* bOmega, sparseMatrix* COmega,
		sparseVector* yuOmega , int obs) {
	double fixedAlpha , alpha1 , lql , alpha, Bx ;
	/*Check if a the sigma structure should be updated */

	/* put pi.Cbar in sigma (fixed part of beta which does not change when C is deterministic) */
	sigma->vals[obs]->piCar = vxMSparse(pi, prob[1]->Cbar, prob[1]->num->rows);


	/* finding fixed part of alpha equal to -1/2 lambda Q lambda - xiBar .pi - yl.mu3 + yuBar.mu2 */

	/* first  find -1/2 lambda Q lambda = obj -(beta x  - xi(obs) .pi - yl.mu3 + yu(obs).mu2)*/
	/*- xi(obs) .pi - yl.mu3 + yu(obs).mu2*/

	fixedAlpha = -vXv(prob[1]->bBar, pi + 1, NULL, prob[1]->num->rows)
										- vXv(prob[1]->lBar, mu3 + 1, NULL, prob[1]->num->cols) + vXv(prob[1]->uBar, mu2 + 1, NULL, prob[1]->num->cols);


	alpha1 = fixedAlpha - vXvSparse(pi + 1, bOmega) + vXvSparse(mu2 + 1, yuOmega);

	Bx = vXv(sigma->vals[obs]->piCar, cell->candidX + 1, NULL, prob[0]->num->cols);

	lql = getObjective(cell->subprob->model) - Bx - alpha1;


	alpha = fixedAlpha + lql;
	sigma->vals[obs]->fixed = alpha + Bx;
	sigma->vals[obs]->interceptBar = alpha;
	return 0;
}

