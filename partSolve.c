#include "stochasticQP.h"

oneCut* partSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset) {
	double tol = 0.0000001;
	double alpha;
	int lambdaIdx;
	int* index = (iVector)arr_alloc(prob->num->prevCols + 1, int);
	clock_t Start;
	clock_t End;


	for (int i = 1; i <= prob->num->prevCols; i++) {
		index[i] = i;
	}
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

	bool flag = 0;
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
	long long int* basis = (long long int*) arr_alloc(prob->num->cols + 1, long long int);
	Buildbase(basis, prob->num->cols, 3);
	/* calculate deltaalpha and delta beta for all existing partitions*/

	/* calculate delta x*/
	dVector dx = (double*)arr_alloc(prob->num->prevCols + 1, double);
	for (int i = 0; i <= prob->num->prevCols; i++) {
		dx[i] = cell->candidX[i] - cell->incumbX[i];
	}
	int size = 1000;
	double*** dnu;
	double*** dmu;
	double*** dy;
	double*** dld;
	pixbCType** deltax; /*size has to be fixed- FREE*/
	dnu = (double**)arr_alloc(size, double*); /*size has to be fixed- FREE*/;
	dmu = (double**)arr_alloc(size, double*); /*size has to be fixed- FREE*/;
	dy = (double**)arr_alloc(size, double*); /*size has to be fixed- FREE*/;
	dld = (double**)arr_alloc(size, double*); /*size has to be fixed- FREE*/;
	deltax = (pixbCType**)arr_alloc(size, pixbCType*); /*size has to be fixed- FREE*/

	for (int i = 0; i <= prob->num->prevCols; i++) {
		dx[i] = cell->candidX[i] - cell->incumbX[i];
	}
	for (int p = 0; p < cell->partition->cnt; p++) {
		AddtoDettaX(prob, cell, deltax, dx, cell->partition->low[p], cell->partition->up[p], cell->partition->inact[p], p, dy, dld, dnu, dmu);
	}
	/* 4. loop through subset omegaP and solve the subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		if (omegaP[obs]) {
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0];
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1];
			dOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[2];
			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3];
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4];
			solnType* soln = buildSolnType(prob->num);

			Start = clock();

			/* 4a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
			if (solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], bOmega, COmega, dOmega, lOmega, uOmega, soln)) {
				errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
				goto TERMINATE;
			}

			End = clock();
			cell->Tsub = cell->Tsub + (End - Start);


			Start = clock();

			StocUpdatePart(cell, prob, bOmega, COmega, lOmega, uOmega, soln, basis, &partIndx, dx, deltax, dy, dld, dnu, dmu);

			End = clock();

			cell->stochupdate = cell->stochupdate + (End - Start)/ CLOCKS_PER_SEC;
			//printf("\n %f", cell->stochupdate );

			/*4b. Calculate observations specific coefficients. */
			double* beta = (double*)arr_alloc(prob->num->prevCols + 1, double);

			alpha = cell->sigma->vals[partIndx]->alpha + cell->delta->vals[partIndx][obs]->alpha + deltax[partIndx]->alpha; /*to do*/

			for (int c = 1; c <= prob->num->prevCols; c++)
				beta[c] += cell->sigma->vals[partIndx]->beta[c];
			for (int c = 1; c <= prob->num->prevCols; c++)
				beta[c] += cell->delta->vals[partIndx][obs]->beta[c] + deltax[partIndx]->beta[c];


			double obj;
			obj = getObjective(cell->subprob->model);
#if defined(STOCH_CHECK)
			//printf("Reconstructed objective function (exact)  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
			//printf("Dif  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols)-obj );
#endif

			/* 4c. Aggregate the cut coefficients by weighting by observation probability. */
			cut->alpha += cell->omega->probs[obs] * alpha;
			for (int c = 1; c <= prob->num->prevCols; c++) {
				cut->beta[c] += cell->omega->probs[obs] * beta[c];
			}
			mem_free(beta);
			freeSolnType(soln);
		}
	}

	/* loop through the rest of the observations */
	/* 4. loop through subset omegaP and use argmax on subproblems */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		if (!omegaP[obs]) {
			/* 4a. Identify the best dual using the argmax operation */


			lambdaIdx = argmaxPart(cell, tol, prob, cell->sigma, cell->delta, cell->candidX, obs, prob->num->prevCols, index, &flag, deltax, dy, dld, dnu, dmu, bOmega, lOmega, uOmega);

			if (flag == 0) {
				cell->IterPart = cell->numit;
				printf("I have to solve\n");

				bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0];
				COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1];
				dOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[2];
				uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3];
				lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4];
				solnType* soln = buildSolnType(prob->num);
				Start = clock();
				if (solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], bOmega, COmega, dOmega, lOmega, uOmega, soln)) {
					errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
					goto TERMINATE;
				}
				End = clock();
				cell->Tsub = cell->Tsub + (End - Start);
				Start = clock();
				StocUpdatePart(cell, prob, bOmega, COmega, lOmega, uOmega, soln, basis, &lambdaIdx, dx, deltax, dy, dld, dnu, dmu);
				End = clock();
				cell->stochupdate = cell->stochupdate + (End - Start);
				freeSolnType(soln);
			}
			double* beta = (double*)arr_alloc(prob->num->prevCols + 1, double);

			alpha = cell->sigma->vals[lambdaIdx]->alpha + cell->delta->vals[lambdaIdx][obs]->alpha + deltax[lambdaIdx]->alpha; /*to do*/

			for (int c = 1; c <= prob->num->prevCols; c++)
				beta[c] += cell->sigma->vals[lambdaIdx]->beta[c];
			for (int c = 1; c <= prob->num->prevCols; c++)
				beta[c] += cell->delta->vals[lambdaIdx][obs]->beta[c] + deltax[lambdaIdx]->beta[c];



			/* 4b. Calculate observations specific coefficients. */

#if defined(STOCH_CHECK)
			//bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0] ;
			//COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1] ;
			//dOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[2] ;
			//uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3] ;
			//lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4] ;

			/* 3a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
			//solnType* soln = buildSolnType(prob->num);
			//if (solveSubprob(prob, cell->subprob, cell->candidX, cell->omega->vals[obs], bOmega, COmega, dOmega, lOmega, uOmega, soln)) {
				///errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
				//goto TERMINATE;
			//}
			//double obj;
			//obj = getObjective(cell->subprob->model);
			//printf("Reconstructed objective function (approx) = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols));
			//printf("Dif  = %lf\n", alpha - vXv(cell->candidX, beta, NULL, prob->num->prevCols) -obj);

#endif

			/* 4c. Aggregate the cut coefficients by weighting by observation probability. */
			cut->alpha += cell->omega->probs[obs] * alpha;
			for (int c = 1; c <= prob->num->prevCols; c++) {
				cut->beta[c] += cell->omega->probs[obs] * beta[c];
			}
			mem_free(beta);
			//freeSolnType(soln);
		}
	}

	mem_free(omegaP);
	mem_free(bOmega);
	mem_free(index);
	mem_free(dOmega);
	mem_free(uOmega);
	mem_free(lOmega);
	mem_free(COmega);
	mem_free(basis);
	mem_free(dx);

	for (int i = 0; i < cell->partition->cnt; i++) {
		mem_free(deltax[i]->beta);
		mem_free(dy[i]);
		mem_free(dmu[i]);
		mem_free(dld[i]);
		mem_free(dnu[i]);
		mem_free(deltax[i]);
	}
	mem_free(deltax);
	mem_free(dnu);
	mem_free(dmu);
	mem_free(dy);
	mem_free(dld);




	return cut;

TERMINATE:
	return NULL;
}//END partSolve()













int StocUpdatePart(cellType* cell, probType* prob, sparseVector* bOmega, sparseMatrix* COmega, sparseVector* lOmega,
	sparseVector* uOmega, solnType* soln, long long int* basis, int* partIndx, dVector dx, pixbCType** deltax, double** dy, double** dld, double** dnu, double** dmu) {
	int up = 0, inact = 0, low = 0; /*Number of variables on their bounds*/
	bool newPartFlag = false;
	Mat* W;
	Mat* T;
	dVector lStat, uStat;

	lStat = (dVector)arr_alloc(prob->num->cols + 1, double);
	if (getDoubleAttributeArray(cell->subprob->model, "LB", 0, prob->num->cols, lStat + 1)) {
		errMsg("solver", "AddtoPart", "failed to obtain variable lower bounds", 0);

	}

	uStat = (dVector)arr_alloc(prob->num->cols + 1, double);
	if (getDoubleAttributeArray(cell->subprob->model, "UB", 0, prob->num->cols, uStat + 1)) {
		errMsg("solver", "AddtoPart", "failed to obtain variable upper bounds", 0);

	}

	/* 4b. Calculate the partition */
	(*partIndx) = addtoPartition(prob, cell, uOmega, lOmega, soln, &newPartFlag, &up, &inact, &low, basis, lStat, uStat);

	if (inact < prob->num->rows) {
		printf("inact is less than M");
	}
	/* 4d. Store the fixed parts of current partition if needed*/

	if (newPartFlag) {

		/* 4d.1 Extract the WT matrices*/
		CalC(cell, prob, prob->sp->objQ, prob->Dbar, &W, &T, low, up, inact);

		Mat* WT = CombineWT(prob, W, T, low, up, inact);

		/* Add the matrix wt to a structure */

		cell->wtSet->wt[(*partIndx)] = WT;
		cell->wtSet->cnt++;
		//showmat(WT);
		/* 4d.2 Add the obtained solution to the lambda structure*/

		addtoLambdaP(cell, soln, WT, prob, bOmega, uOmega, lOmega, low, up, inact, dx);

		/* 4d.2 Add to alpha and beta Bar*/
		AddtoSigmaP(cell, soln, prob);
		AddtoDettaX(prob, cell, deltax, dx, cell->partition->low[(*partIndx)], cell->partition->up[partIndx[0]], cell->partition->inact[partIndx[0]], partIndx[0], dy, dld, dnu, dmu);

		cell->delta->dy[(*partIndx)] = (double**)arr_alloc(cell->omega->cnt, double*);
		cell->delta->dmu[(*partIndx)] = (double**)arr_alloc(cell->omega->cnt, double*);
		cell->delta->dnu[(*partIndx)] = (double**)arr_alloc(cell->omega->cnt, double*);
		/* 4d.3 add to  delta sol and complete a row*/
		for (int obs = 0; obs < cell->omega->cnt; obs++) {
			/* Add a new row to the delta structure for all observations and the latest lambda (lambdaIdx) */
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0];
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1];
			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3];
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4];
			addtoDeltaP(cell, soln, W, T, WT, prob, COmega, bOmega, uOmega, lOmega, obs, (*partIndx), inact, up, low, lStat, uStat, dy[(*partIndx)], dld[(*partIndx)], dnu[(*partIndx)], dmu[(*partIndx)]);

		}

		freemat(W);
		freemat(T);
	}
	mem_free(lStat); mem_free(uStat);
	return 0;
}; //EndStocUpdatePart









int argmaxPart(cellType* cell, double tol, probType* prob, sigmaType* sigma, deltaType* delta, dVector Xvect, int obs, int prev, int* index, bool* flag, pixbCType** deltaxp, double** dy, double** dld, double** dnu, double** dmu, sparseVector* bomega, sparseVector* lomega, sparseVector* uomega) {
	double maxobj = -INFINITY;
	flag[0] = false;
	int lambdaIdx = -1;
	double tempobj = 0;
	double dalphx = 0;
	int cnt = 0;
	int feas = 1;
	int rand = 0;


	for (int j = 0; j < sigma->cnt; j++) {
		feas = 1;
		//check feasibility for y
		for (int i = 1; i <= prob->num->cols; i++) {
			if (cell->partition->part[j][i] == 0) {
				rand = 0;
				for (int z = 1; z <= uomega->cnt; z++) {

					if (uomega->col[z] == i) {
						rand = 1;
						if (cell->lambda->y[j][i] + delta->dy[j][obs][i] + dy[j][i] > uomega->val[z] + prob->uBar->val[i] + tol) {
							feas = 0;
						}
					}
				}
				if (rand == 0) {
					if (cell->lambda->y[j][i] + delta->dy[j][obs][i] + dy[j][i] > prob->uBar->val[i] + tol) {
						feas = 0;
					}
				}

			}
		}




		//for mu

		for (int i = 1; i <= prob->num->cols; i++) {

			if (cell->partition->part[j][i] == 2) {

				if (cell->lambda->umu[j][i] + delta->dmu[j][obs][i] + dmu[j][i] < -tol) {
					feas = 0;
				}

			}
		}

		for (int i = 1; i <= prob->num->cols; i++) {

			if (cell->partition->part[j][i] == 1) {

				if (cell->lambda->lmu[j][i] + delta->dnu[j][obs][i] + dnu[j][i] < -tol) {
					feas = 0;
				}

			}
		}


		if (feas == 1) {
			flag[0] = true;
			//dVector temp1, temp2;
			//temp1 = vxMSparse(dy[j], prob->sp->objQ, prob->num->cols);
			//temp2 = vxMSparse(delta->dy[j][obs], prob->sp->objQ, prob->num->cols);

			//dalphx = -vXv(temp1, delta->dy[j][obs], index, prob->num->cols) - vXv(temp2, dy[j], index, prob->num->cols)+ vXvSparse( dld, bomega ) + vXvSparse(dnu, lomega)  - vXvSparse(dmu, uomega);
			//mem_free(temp1);
			//mem_free(temp2);


			tempobj = sigma->vals[j]->alpha + delta->vals[j][obs]->alpha + deltaxp[j]->alpha +
				-vXv(Xvect, sigma->vals[j]->beta, index, prev)
				- vXv(Xvect, delta->vals[j][obs]->beta, index, prev)
				- vXv(Xvect, deltaxp[j]->beta, index, prev);
			/*4b. calculate estimated Obj value beta x + alpha*/
			if (tempobj > maxobj) {
				maxobj = tempobj;
				lambdaIdx = j;
			}
		}
	}



	if (flag[0] == false) {
		printf("nofeas");
	}












	return lambdaIdx;
}//END argmax()















void AddtoDettaX(probType* prob, cellType* cell, pixbCType** delta, dVector deltaX, int low, int up, int inact, int partindex, double** deltay, double** deltaLd, double** dnu, double** dmu) {

	int* index = (iVector)arr_alloc(prob->num->prevCols + 1, int);
	int* index2 = (iVector)arr_alloc(prob->num->cols + 1, int);
	for (int i = 1; i <= prob->num->prevCols; i++) {
		index[i] = i;
	}
	for (int i = 1; i <= prob->num->cols; i++) {
		index2[i] = i;
	}

	int elm = 0;

	Mat* deltaPD;

	deltay[partindex] = (dVector)arr_alloc(prob->num->cols + 1, double);

	deltaLd[partindex] = (dVector)arr_alloc(prob->num->rows + 1, double);

	delta[partindex] = (pixbCType*)mem_malloc(sizeof(pixbCType));

	dVector dbeta;
	dnu[partindex] = (dVector)arr_alloc(prob->num->cols + 1, double);
	dmu[partindex] = (dVector)arr_alloc(prob->num->cols + 1, double);
	/* 1. Calculate  vector [delta ubarU , delta rho] */

	/* 1.a. caculate -cbardeltax*/

	dVector deltarho = vxMSparse(deltaX, prob->Cbar, prob->num->rows);
	Mat* rhyu = newmat(prob->num->rows + up, 1, 0);


	for (int i = 1; i <= prob->num->rows; i++) {
		rhyu->entries[up + i - 1] = -deltarho[i];
	}

	/* 2. Calculate the delta [yI,lambda,nuL,muU] *//* calculate wt -cbardeltax for partition i*/

	deltaPD = multiply(cell->wtSet->wt[partindex], rhyu);
	//showmat(deltaPD);
	/* 3. put delta lambda in a vector deltaLd */

	int cnt = 0;

	copyVector(deltaPD->entries + inact, deltaLd[partindex] + 1, prob->num->rows - 1);
	//showmat(cell->wtSet->wt[partindex]);
	//showmat(rhyu);



	/* Calculate delta beta */

	dbeta = vxMSparse(deltaLd[partindex], prob->Cbar, prob->num->prevCols);

	delta[partindex]->beta = dbeta;

	/* Calculate deltaY */
	cnt = 0;
	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[partindex][i] == 1) {
			deltay[partindex][i] = 0;
		}

		if (cell->partition->part[partindex][i] == 2) {
			deltay[partindex][i] = 0;
		}

		if (cell->partition->part[partindex][i] == 0) {
			deltay[partindex][i] = deltaPD->entries[cnt];
			cnt++;
		}
	}

	/* Calculate deltanu */
	cnt = 0;

	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[partindex][i] == 0 || cell->partition->part[partindex][i] == 2) {
			dnu[partindex][i] = 0;
		}
		if (cell->partition->part[partindex][i] == 1) {
			dnu[partindex][i] = deltaPD->entries[inact + prob->num->rows + cnt];
			cnt++;
		}
	}

	/* Calculate delta mu */
	cnt = 0;


	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[partindex][i] == 0 || cell->partition->part[partindex][i] == 1) {
			dmu[partindex][i] = 0;
		}
		if (cell->partition->part[partindex][i] == 2) {
			dmu[partindex][i] = deltaPD->entries[inact + low + prob->num->rows + cnt];
			cnt++;
		}
	}

	/* Calculate alpha */

	dVector temp4, temp5;
	temp4 = vxMSparse(deltay[partindex], prob->sp->objQ, prob->num->cols);
	temp5 = vxMSparse(cell->lambda->y[partindex], prob->sp->objQ, prob->num->cols);
	delta[partindex]->alpha =
		-vXv(temp4, deltay[partindex], index2, prob->num->cols)
		- vXv(temp5, deltay[partindex], index2, prob->num->cols)
		- vXv(temp4, cell->lambda->y[partindex], index2, prob->num->cols)
		+ vXvSparse(deltaLd[partindex], prob->bBar)
		+ vXvSparse(dnu[partindex], prob->lBar)
		- vXvSparse(dmu[partindex], prob->uBar);


	mem_free(temp4); mem_free(temp5);
	mem_free(index);
	mem_free(index2);
	freemat(rhyu);
	freemat(deltaPD);
	mem_free(deltarho);


}//END DelttaX()





