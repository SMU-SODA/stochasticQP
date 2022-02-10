#include "stochasticQP.h"










oneCut* dualSolve(probType** prob, cellType* cell, stocType* stoch ,double* x, double solveset) {

	/*initialization of the parameters*/
	sparseMatrix *COmega; /* Presenting the C matrix associated with an observation(I mean the difference from Cbar)*/
	sparseVector *bOmega;  /* Presenting the b vector associated with an observation(I mean the difference from bBar)*/
	bOmega = (sparseVector*)mem_malloc(sizeof(sparseVector));
	COmega = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));
	bOmega->cnt = prob[1]->num->rvbOmCnt;
	bOmega->col = prob[1]->coord->rvbOmRows;
	COmega->cnt = prob[1]->num->rvCOmCnt;
	COmega->col = prob[1]->coord->rvCOmCols;
	COmega->row = prob[1]->coord->rvCOmRows;
	sparseVector* ybar; /* Presenting the upperbound  vector associated with an observation(I mean the difference from yBar)*/
	sparseVector* yund; /* Presenting the lowerbound  vector associated with an observation(I mean the difference from mean of yunderscore)*/
	ybar = (sparseVector*)mem_malloc(sizeof(sparseVector));
	yund = (sparseVector*)mem_malloc(sizeof(sparseVector));
	ybar->cnt = prob[1]->num->rvyuOmCnt; ybar->col = prob[1]->coord->rvyuOmRows;
	yund->cnt = prob[1]->num->rvbOmCnt; yund->col = prob[1]->coord->rvylOmRows;
	double tempobj; /* holds the best dual objective function value when we are searching over lambda for a suboptimal solution*/
	double alpha; /*presenting an alpha value we obtained using a dual solution while searching through lambda*/
	double* beta; /*presenting an beta value we obtained using a dual solution while searching through lambda*/
	double* dbeta; /**presenting delta beta while searching through lambda*/
	int *omegaP; /*representing the set of opservations we want to solve exactly*/
	int* omegaQ; /*representing the set of opservations we want to solve approximatly*/
	omegaP = (int*)arr_alloc(solveset, int); 
	omegaQ = (int*)arr_alloc(cell->omega->cnt - solveset,int);
	double  bestalpha;
	double*  bestbetha;
	double * s;
	double dalpha = 0;
	double* Y;
	Y = (double*)arr_alloc(prob[1]->num->cols+1, double);
	int cnt;
	double 		 mubBar;
	dVector pi;
	pi = (dVector)arr_alloc(prob[1]->num->rows + 1, double);

	 /* 1. Initialize a new cut */
	 oneCut* cut = newCut(prob[1]->num->cols);
	
	 /*2. Creat a subset of observations you want to solve*/
	 sample(omegaP , solveset , cell->omega->cnt);

	/*Put the remaining observations in another set omegaQ*/
	 subtractSample(omegaP, omegaQ , cell->omega->cnt , solveset);

	/* 3. loop through subset solveset and solve the subproblems */

	 for (int i = 0; i < solveset; i++) {
		 /* 2a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
		 cnt = cell->lambda->cnt;
		 if (solveSubprob(prob[1], cell->subprob, cell->candidX, cell->omega->vals[i], cell->lambda->pi[i], &mubBar, cell->lambda->mu2[i]
			 , cell->lambda->mu3[i])) {
			 errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			 goto TERMINATE;
		 }
	 }
 //    /*Calculate fixred section of beta = Cbar*pi (pi is the dual vector associated with equality constraints)*/
	//sigma->vals[obs]->piCar = vxMSparse(pi, prob[1]->Cbar, prob[1]->num->prevCols);

	///*Calculate fixred section of alpha -.5 yQy  + bbar* pi + yunder nu - ybar mu)*/
	//sigma->vals[obs]->interceptBar= -0.5*vXv(vxMSparse(Y, prob[1]->sp->objQ, prob[1]->num->cols),Y, indexx , prob[1]->num->cols) + vXvSparse(pi, prob[1]->bBar) + mubBar ;
	//sigma->cnt++;
	///* Store the information of the dual solutions obtaned in matrix lambda*/
	//lambda->pi = pi  ;
	//}
	//for(int obs = 0; obs < cell->omega->cnt-subset; obs++) {
	//	maxobj = -1000;
	//	for (int i = 0; i < sigma->cnt ;i++) {

	//		/* 2b. Compute the cut coefficients for individual observation. */
	//		bOmega->val = cell->omega->vals[obs] + prob[1]->coord->rvOffset[0];
	//		COmega->val = cell->omega->vals[obs] + prob[1]->coord->rvOffset[1];
	//		ybar->val = cell->omega->vals[obs] + prob[1]->coord->rvOffset[2];
	//		yund->val = cell->omega->vals[obs] + prob[1]->coord->rvOffset[3];

	//	/* calculate alpha */
	//	dalpha = vXvSparse(lambda->pi[i], bOmega ,prob[1]->num->rows) - vXvSparse(lambda->mu2[i], yund, prob[1]->num->cols) + vXvSparse(lambda->mu3[i], ybar, prob[1]->num->cols);
	//	alpha = sigma->vals[i]->interceptBar + dalpha ;

 //        /* Calculate beta*/
	//	dbeta = vxMSparse(lambda->pi[i], COmega, prob[1]->num->prevCols);
	//	int* dbetaindex;
	//	dbetaindex = (int*)arr_alloc(prob[1]->num->prevCols, int);

	//	/*find out the number of nonzero elements and their location*/
	//	for(int j = 0; j < prob[1]->num->prevCols ;j++) {
	//		if (dbeta[j]) {
	//			delta->vals[j]->dbeta->cnt++;
	//			}
	//		
	//	/*Assign memory to the index section of deltab */

	//		delta->vals[i]->dbeta->col = (int*)mem_malloc(prob[1]->num->prevCols + 1, int);
	//	    delta->vals[i]->dbeta->val= (double*)mem_malloc(prob[1]->num->prevCols + 1, double);
	//		int cnt = 0;
	//		for (int j = 0; j < prob[1]->num->prevCols; j++) {
	//			if (dbeta[j]) {
	//				cnt++;
	//				delta->vals[i]->dbeta->val[cnt]= dbeta[j];
	//				delta->vals[i]->dbeta->col[cnt]=j;
	//			}
	//		}
	//		delta->vals[i]->dalpha = dalpha;

	//	

	//	vPlusv(dbeta, sigma->vals[i]->piCar, 1, prob[1]->num->prevCols);
	//	beta = dbeta;

	//	tempobj = alpha + vXv(beta, cell->candidX, indexx, prob[1]->num->prevCols);
 //       /*calculate estimated Obj value beta x+alpha*/
	//	if (tempobj> maxobj) {
	//	maxobj = alpha + vXv(beta, cell->candidX, indexx , prob[1]->num->prevCols);
	//	bestalpha = alpha;
	//	bestbetha = beta;
	//	}
	//	/* 2c. Aggregate the cut coefficients by weighting by observation probability. */
	//	cut->alpha += cell->omega->probs[obs] * bestalpha;
	//	for (int c = 1; c <= prob[1]->num->prevCols; c++)
	//		cut->beta[c] += cell->omega->probs[obs] * bestbetha[c];

	//	mem_free(bestbetha);
	//}
	//}
	 mem_free(omegaP);
	 mem_free(omegaQ);
    TERMINATE:
	return NULL;
}//END dualsolve()


