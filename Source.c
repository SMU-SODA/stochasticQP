/*
 * algorithm.c
 *  Created on: Jul , 2021
 *      Author: Niloofar Fadavi, Harsha Gangammanavar
 */
#include "./solverUtilities/utilities.h"
#include "./solverUtilities/solver_gurobi.h"
#include "./smpsReader/smps.h"
#include "./smpsReader/prob.h"

#include "stochasticQP.h"
extern configType config;


/*building master problem*/
cellType* buildCell(probType** prob , omegaType* omega) {
	
	oneProblem* stage0 = NULL;
    stage0 = (oneProblem*)mem_malloc(sizeof(oneProblem));
	stage0->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
	stage0->objQ->col = (iVector)arr_alloc(prob[0]->sp->mac * prob[0]->sp->mac, int);
	stage0->objQ->row = (iVector)arr_alloc(prob[0]->sp->mac * prob[0]->sp->mac, int);
	stage0->objQ->val = (dVector)arr_alloc(prob[0]->sp->mac * prob[0]->sp->mac, double);
	stage0->model = NULL;
	stage0->name = (cString)arr_alloc(NAMESIZE, char);
	stage0->objname = (cString)arr_alloc(NAMESIZE, char);
	stage0->objx = (dVector)arr_alloc(prob[0]->sp->macsz, double);
	stage0->bdl = (dVector)arr_alloc(prob[0]->sp->macsz, double);
	stage0->bdu = (dVector)arr_alloc(prob[0]->sp->macsz, double);
	stage0->ctype = (cString)arr_alloc(prob[0]->sp->macsz, char);
	stage0->rhsx = (dVector)arr_alloc(prob[0]->sp->marsz, double);
	stage0->senx = (cString)arr_alloc(prob[0]->sp->marsz, char);
	stage0->matbeg = (iVector)arr_alloc(prob[0]->sp->macsz, int);
	stage0->matcnt = (iVector)arr_alloc(prob[0]->sp->macsz, int);
	stage0->matval = (dVector)arr_alloc(prob[0]->sp->matsz, double);
	stage0->matind = (iVector)arr_alloc(prob[0]->sp->matsz, int);
	stage0->cname = (cString*)arr_alloc(prob[0]->sp->macsz, cString);
	stage0->rname = (cString*)arr_alloc(prob[0]->sp->marsz, cString);
	stage0->mac = prob[0]->sp->mac;
	stage0->type= prob[0]->sp->type;			/* type of problem: LP, QP, MIP or MIQP */

	for (int i = 0; i < prob[0]->sp->mac; i++) {
		stage0->objx[i] = prob[0]->sp->objx[i];
	}

	stage0->objQ->cnt = prob[0]->sp->objQ->cnt;

	for (int i = 0; i < prob[0]->sp->objQ->cnt; i++) {
		stage0->objQ->val[i] = prob[0]->sp->objQ->val[i];
		stage0->objQ->col[i] = prob[0]->sp->objQ->col[i];
	}

	stage0->objSense = prob[0]->sp->objSense;		
	stage0->mac = prob[0]->sp->mac;			/* number of columns */
	stage0->mar = prob[0]->sp->mar;			/* number of rows */
	stage0->numBin = prob[0]->sp->numBin;			/* number of binary variables in the problem */
	stage0->numInt = prob[0]->sp->numInt;			/* number of integer variables  (includes both binary and general integer variables) in the problem */
	stage0->numnz = prob[0]->sp->numnz;			/* number of non-zero elements in constraint matrix */
	stage0->macsz = prob[0]->sp->macsz;			/* number of columns */
	stage0->marsz = prob[0]->sp->marsz;			/* number of rows */
	stage0->matsz = prob[0]->sp->matsz;			/* number of rows */
	

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->bdl[i] = prob[0]->sp->bdl[i];
	}

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->bdu[i] = prob[0]->sp->bdu[i];
	}

	for (int i = 0; i < prob[0]->sp->marsz; i++) {
		stage0->rhsx[i] = prob[0]->sp->rhsx[i];
	}

	for (int i = 0; i < prob[0]->sp->marsz; i++) {
		stage0->senx[i] = prob[0]->sp->senx[i];
	}

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->matbeg[i] = prob[0]->sp->matbeg[i];
	}

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->matcnt[i] = prob[0]->sp->matcnt[i];
	}

	for (int i = 0; i < prob[0]->sp->matsz; i++) {
		stage0->matind[i] = prob[0]->sp->matind[i];
	}


	for (int i = 0; i < prob[0]->sp->matsz; i++) {
		stage0->matval[i] = prob[0]->sp->matval[i];
	}


	for (int i = 0; i < NAMESIZE; i++) {
		stage0->objname[i] = prob[0]->sp->objname[i];
	}

	for (int i = 0; i < NAMESIZE; i++) {
		stage0->name[i] = prob[0]->sp->name[i];
	}

	for (int i = 0; i < prob[0]->sp->marsz; i++) {
		stage0->rname[i] = prob[0]->sp->rname[i];
	}
	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->cname[i] = prob[0]->sp->cname[i];
	}



	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->ctype[i] = prob[0]->sp->ctype[i];
	}
	cString	ctype;			/* type of decision variables: 'C' continuous, 'B' binary, 'I' general integer, 'S' semi-continuous, 'N' semi-integer */
	
	/* These don't seem necessary for Gurobi. */
	//	int		rstorsz;		/* memory size for storing row names */
	//	cString	rstore;			/* row names cString */
	//	int		cstorsz;		/* memory size for storing column names */
	//	cString	cstore;			/* column name cString */

	cellType* prb;
	prb = (cellType*)mem_malloc(sizeof(cellType));

	prb->master = stage0;

	prb->master = newSubprob(prb->master);

	
	oneProblem* stage1 = NULL;
	stage1 = (oneProblem*)mem_malloc(sizeof(oneProblem));
	stage1->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
	stage1->objQ->col = (iVector)arr_alloc(prob[1]->sp->mac * prob[1]->sp->mac, int);
	stage1->objQ->row = (iVector)arr_alloc(prob[1]->sp->mac * prob[1]->sp->mac, int);
	stage1->objQ->val = (dVector)arr_alloc(prob[1]->sp->mac * prob[1]->sp->mac, double);
	stage1->model = NULL;
	stage1->name = (cString)arr_alloc(NAMESIZE, char);
	stage1->objname = (cString)arr_alloc(NAMESIZE, char);
	stage1->objx = (dVector)arr_alloc(prob[1]->sp->macsz, double);
	stage1->bdl = (dVector)arr_alloc(prob[1]->sp->macsz, double);
	stage1->bdu = (dVector)arr_alloc(prob[1]->sp->macsz, double);
	stage1->ctype = (cString)arr_alloc(prob[1]->sp->macsz, char);
	stage1->rhsx = (dVector)arr_alloc(prob[1]->sp->marsz, double);
	stage1->senx = (cString)arr_alloc(prob[1]->sp->marsz, char);
	stage1->matbeg = (iVector)arr_alloc(prob[1]->sp->macsz, int);
	stage1->matcnt = (iVector)arr_alloc(prob[1]->sp->macsz, int);
	stage1->matval = (dVector)arr_alloc(prob[1]->sp->matsz, double);
	stage1->matind = (iVector)arr_alloc(prob[1]->sp->matsz, int);
	stage1->cname = (cString*)arr_alloc(prob[1]->sp->macsz, cString);
	stage1->rname = (cString*)arr_alloc(prob[1]->sp->marsz, cString);
	stage1->mac = prob[1]->sp->mac;
	stage1->type = prob[1]->sp->type;			/* type of problem: LP, QP, MIP or MIQP */
	for (int i = 0; i < prob[1]->sp->mac; i++) {
		stage1->objx[i] = prob[1]->sp->objx[i];
	}

	stage1->objQ->cnt = prob[1]->sp->objQ->cnt;

	for (int i = 0; i < prob[1]->sp->objQ->cnt; i++) {
		stage1->objQ->val[i] = prob[1]->sp->objQ->val[i];
		stage1->objQ->col[i] = prob[1]->sp->objQ->col[i];
	}

	stage1->objSense = prob[1]->sp->objSense;
	stage1->mac = prob[1]->sp->mac;			/* number of columns */
	stage1->mar = prob[1]->sp->mar;			/* number of rows */
	stage1->numBin = prob[1]->sp->numBin;			/* number of binary variables in the problem */
	stage1->numInt = prob[1]->sp->numInt;			/* number of integer variables  (includes both binary and general integer variables) in the problem */
	stage1->numnz = prob[1]->sp->numnz;			/* number of non-zero elements in constraint matrix */
	stage1->macsz = prob[1]->sp->macsz;			/* number of columns */
	stage1->marsz = prob[1]->sp->marsz;			/* number of rows */
	stage1->matsz = prob[1]->sp->matsz;			/* number of rows */


	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->bdl[i] = prob[1]->sp->bdl[i];
	}

	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->bdu[i] = prob[1]->sp->bdu[i];
	}

	for (int i = 0; i < prob[1]->sp->marsz; i++) {
		stage1->rhsx[i] = prob[1]->sp->rhsx[i];
	}

	for (int i = 0; i < prob[1]->sp->marsz; i++) {
		stage1->senx[i] = prob[1]->sp->senx[i];
	}

	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->matbeg[i] = prob[1]->sp->matbeg[i];
	}

	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->matcnt[i] = prob[1]->sp->matcnt[i];
	}

	for (int i = 0; i < prob[1]->sp->matsz; i++) {
		stage1->matind[i] = prob[1]->sp->matind[i];
	}


	for (int i = 0; i < prob[1]->sp->matsz; i++) {
		stage1->matval[i] = prob[1]->sp->matval[i];
	}


	for (int i = 0; i < NAMESIZE; i++) {
		stage1->objname[i] = prob[1]->sp->objname[i];
	}

	for (int i = 0; i < NAMESIZE; i++) {
		stage1->name[i] = prob[1]->sp->name[i];
	}

	for (int i = 0; i < prob[1]->sp->marsz; i++) {
		stage1->rname[i] = prob[1]->sp->rname[i];
	}
	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->cname[i] = prob[1]->sp->cname[i];
	}



	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->ctype[i] = prob[1]->sp->ctype[i];
	}

	/* These don't seem necessary for Gurobi. */
	//	int		rstorsz;		/* memory size for storing row names */
	//	cString	rstore;			/* row names cString */
	//	int		cstorsz;		/* memory size for storing column names */
	//	cString	cstore;			/* column name cString */
	prb->subprob = stage1;

	prb->subprob = newSubprob(prb->subprob);

	prb->omega = omega;

	return prb;
	}


/*This function gets the stoc file and creats a matrix of observations*/

int numObs(stocType* stoch) {
	int numObs=1; 
	int coef = 1;
	int precoef = 1;

	for (int i = 0; i < stoch->numOmega; ++i) {
		numObs = numObs * stoch->numVals[i];
	} /*number of observations*/



	//dVector probs;
	//probs = (dVector)arr_alloc(numObs, double);
	//for (int i = 0; i < numObs; ++i) {
	//	probs[i] = 1;
	//}
	//
	//for (int i = stoch->numOmega; i > 0; --i) {
	//	precoef = coef;
	//	coef = coef * stoch->numVals[i-1];
	//	int f = numObs / coef; /*number in groups*/

	//	for (int j = 0; j < f; ++j) {
	//		for (int k = 0; k < stoch->numVals[i-1]; ++k) {

	//			for (int t = 0; t < precoef; ++t) {
	//				int ind1 = j * stoch->numVals[i-1] + k * precoef + t;
	//				
	//				probs[ ind1] = probs[ind1] * stoch->probs[i-1][k];
	//			}
	//		}
	//	}
	//}
	return numObs;
}
/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a dVector to hold an array of
 * observation and the probability associated with it. */

omegaType* newOmega(stocType* stoc, int numObs) {
	omegaType* omega;
	int cnt, i, base, idx;

	omega = (omegaType*)mem_malloc(sizeof(omegaType));
	omega->probs = (dVector)arr_alloc(numObs, double);
	omega->weights = (iVector)arr_alloc(numObs, int);
	omega->vals = (dVector*)arr_alloc(numObs, dVector);
	omega->cnt = 0; omega->numRV = stoc->numOmega;

	

	if (config.SAA == 1 ) {
		config.SAA = 1;
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
					omega->vals[cnt][i + 1] = stoc->vals[i][idx] - stoc->mean[i];
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


 // This function gets a oneproblem type and changes the rhs 

oneProblem* setRhs(oneProblem* subProb, dVector rhs) {
	for (int i = 0; i < subProb->mac; i++)
	{
		subProb->rhsx[i] = rhs[i];
	};
	return subProb;
}


///* This function will solve a new subproblem. This involves replacing the right-hand side of the subproblem with new values, based upon some
// * observation of omega, and some X dVector of primal variables from the master problem.  Generally, the latest observation is used.  When
// * forming a normal cut, the candidate x should be used, while the incumbent x should be used for updating the incumbent cut. */
//int solveSubprob(probType* prob, oneProblem* subproblem, dVector Xvect, basisType* basis, lambdaType* lambda, sigmaType* sigma, deltaType* delta, int deltaRowLength,
//	omegaType* omega, int omegaIdx, bool* newOmegaFlag, int currentIter, double TOLERANCE, bool* subFeasFlag, bool* newBasisFlag,
//	double* subprobTime, double* argmaxTime) {
//	int  	status;
//	clock_t tic;
//
//	/* (a) compute and change the right-hand side using current observation and first-stage solution */
//	if (computeRHS(subproblem->lp, prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, omega->vals[omegaIdx])) {
//		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
//		return 1;
//	}
//
//	if (prob->num->rvdOmCnt > 0) {
//		/* (b) Compute and change the cost coefficients using current observation */
//		if (computeCostCoeff(subproblem->lp, prob->num, prob->coord, prob->dBar, omega->vals[omegaIdx])) {
//			errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
//			return 1;
//		}
//	}
//
//#if defined(ALGO_CHECK)
//	writeProblem(subproblem->lp, "subproblem.lp");
//#endif
//
//	/* (c) Solve the subproblem to obtain the optimal dual solution. */
//	tic = clock();
//	setIntParam(PARAM_PREIND, OFF);
//	changeLPSolverType(ALG_PRIMAL);
//	if (solveProblem(subproblem->lp, subproblem->name, subproblem->type, &status)) {
//		if (status == STAT_INFEASIBLE) {
//			/* Set the subproblem feasibility flag to false and proceed to complete stochastic updates. These updates are
//			 * used to generate the feasibility cuts later. */
//			printf("Subproblem is infeasible for current first-stage decision and observation.\n");
//			writeProblem(subproblem->lp, "infeasibleSP.lp");
//			(*subFeasFlag) = false;
//		}
//		else {
//			errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
//			return 1;
//		}
//	}
//	setIntParam(PARAM_PREIND, ON);
//	(*subprobTime) += ((double)(clock() - tic)) / CLOCKS_PER_SEC;
//
//#ifdef STOCH_CHECK
//	double obj;
//	obj = getObjective(subproblem->lp, PROB_LP);
//	printf("Objective value of Subproblem  = %lf\n", obj);
//#endif
//
//	if (newBasisFlag != NULL) {
//		tic = clock();
//		/* (d) update the stochastic elements in the problem */
//		status = stochasticUpdates(prob, subproblem->lp, basis, lambda, sigma, delta, deltaRowLength,
//			omega, omegaIdx, (*newOmegaFlag), currentIter, TOLERANCE, newBasisFlag, (*subFeasFlag));
//		(*newOmegaFlag) = false;
//		(*argmaxTime) += ((double)(clock() - tic)) / CLOCKS_PER_SEC;
//
//#ifdef STOCH_CHECK
//		obj = sigma->vals[status].pib - vXv(sigma->vals[status].piC, Xvect, prob->coord->CCols, prob->num->cntCcols);
//		obj += delta->vals[sigma->lambdaIdx[status]][omegaIdx].pib - vXv(delta->vals[sigma->lambdaIdx[status]][omegaIdx].piC,
//			omega->vals[omegaIdx], prob->coord->rvCOmCols, prob->num->rvCOmCnt);
//		printf("Objective function estimate    = %lf\n", obj);
//#endif
//	}
//
//	return 0;
//}// END solveSubprob()




/* This function computes the right hand side of the subproblem, based on a given X dVector and a given observation of omega.
 * It is defined as:
 * 			rhs = R(omega) - T(omega) x X
 * and is calculated as:
 * 			rhs = (Rbar - Tbar x X) + (Romega - Tomega x X)
 *
 * where the "bar" denotes the fixed or mean value, and the "omega" denotes a random variation from this mean. The function allocates an array
 * for the dVector, which must be freed by the customer.  Also, the zeroth position of this rhs dVector is reserved, and the actual values begin at rhs[1].
 * R is b, and T is C
 \***********************************************************************/

//int computeRHS(modelPtr* lp, numType* num, coordType* coord, sparseVector* bBar, sparseMatrix* Cbar, dVector X, dVector obs) {
//	sparseMatrix Comega;
//	sparseVector bomega;
//	dVector rhs;
//	int cnt, * indices;
//
//
//	if (!(indices = (iVector)arr_alloc(num->rows, int)))
//		errMsg("allocation", "solveSubprob", "indices", 0);
//
//
//	for (cnt = 0; cnt < num->rows; cnt++)
//		indices[cnt] = cnt;
//
//	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows; bomega.val = obs + coord->rvOffset[0];
//
//	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols + num->rvbOmCnt;
//	Comega.row = coord->rvCOmRows + num->rvbOmCnt; Comega.val = obs + coord->rvOffset[1];
//
//	/* Start with the values of b(omega) -- both fixed and varying */
//	rhs = expandVector(bBar->val, bBar->col, bBar->cnt, num->rows);
//	for (cnt = 1; cnt <= bomega.cnt; cnt++)
//		rhs[bomega.col[cnt]] += bomega.val[cnt];
//
//	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
//	rhs = MSparsexvSub(Cbar, X, rhs);
//	rhs = MSparsexvSub(&Comega, X, rhs);
//
//	/* change the right-hand side in the solver */
//	if (changeRHS(lp, num->rows, indices, rhs + 1)) {
//		errMsg("solver", "solveSubprob", "failed to change the right-hand side in the solver", 0);
//		return 1;
//	}
//
//	mem_free(indices); mem_free(rhs);
//	return 0;
//}//END computeRHS()
//
//int computeCostCoeff(modelPtr* lp, numType* num, coordType* coord, sparseVector* dBar, dVector observ) {
//	sparseVector dOmega;
//	dVector cost;
//	int	cnt, * indices;
//
//	if (!(indices = (iVector)arr_alloc(num->cols, int)))
//		errMsg("allocation", "solveSubprob", "indices", 0);
//	for (cnt = 0; cnt < num->cols; cnt++)
//		indices[cnt] = cnt;
//
//	dOmega.cnt = num->rvdOmCnt; dOmega.col = coord->rvdOmCols; dOmega.val = coord->rvOffset[2] + observ;
//
//	/* Extract the cost coefficients */
//	cost = expandVector(dBar->val, dBar->col, dBar->cnt, num->cols);
//	for (cnt = 1; cnt <= dOmega.cnt; cnt++)
//		cost[dOmega.col[cnt]] += dOmega.val[cnt];
//
//	/* change cost coefficients in the solver */
//	if (changeObjx(lp, num->cols, indices, cost + 1)) {
//		errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver", 0);
//		return -1;
//	}
//
//	mem_free(indices); mem_free(cost);
//	return 0;
//}//END computeCostCoeff()
//
//  void chgRHSwSoln(sparseVector* bBar, sparseMatrix* Cbar, dVector rhs, dVector X) {
//	int cnt;
//
//	/* copy the original right-hand side */
//	for (cnt = 1; cnt <= bBar->cnt; cnt++)
//		rhs[bBar->col[cnt]] = bBar->val[cnt];
//
//  change the right-hand side with first stage solution */
//	rhs = MSparsexvSub(Cbar, X, rhs);
//
//}//END chgRHSwMean()
//



//int chgRHSwObserv(modelPtr* lp, numType* num, coordType* coord, dVector observ, dVector spRHS, dVector X) {
//
//	sparseVector bomega;
//	sparseMatrix Comega;
//	dVector 	rhs;
//	iVector	indices;
//	int		cnt, stat1;
//
//	bomega.cnt = num->rvbOmCnt;
//	bomega.col = coord->rvbOmRows;
//	bomega.val = observ;
//
//	//Comega.cnt = num->rvCOmCnt; 
//	//Comega.col = coord->rvCOmCols + num->rvbOmCnt;
//	//Comega.row = coord->rvCOmRows + num->rvbOmCnt; 
//	//Comega.val = observ + num->rvbOmCnt;
//
//	if (!(indices = (iVector)arr_alloc(num->rows, int)))
//		errMsg("allocation", "chgRHSwObserv", "indices", 0);
//
//	if (!(rhs = (dVector)arr_alloc(num->rows + 1, double)))
//		errMsg("allocation", "chgRHSwObserv", "rhs", 0);
//
//
//	/* copy right-hand side modified with mean information */
//	for (cnt = 1; cnt <= num->rows; cnt++) {
//		rhs[cnt] = spRHS[cnt];
//		indices[cnt - 1] = cnt - 1;
//	}
//
//
//	/* change right-hand side with randomness in b */
//	for (cnt = 1; cnt <= bomega.cnt; cnt++)
//		rhs[bomega.col[cnt]] += bomega.val[cnt];
//
//
//
//
//	///* change right-hand side with randomness in transfer matrix */
//	//rhs = MSparsexvSub(&Comega, X, rhs);
//
//	/* change the right-hand side in the solver */
//	stat1 = changeRHS(lp, num->rows, indices, rhs + 1);
//	if (stat1) {
//		errMsg("solver", "chgRHSwObserv", "failed to change the right-hand side in the solver", 0);
//		return 1;
//	}
//	mem_free(rhs); mem_free(indices);
//	return 0;
//}//END chgRHSwRand()
//
//
//
//
////int chgObjxwObserv(modelPtr* lp, numType* num, coordType* coord, dVector cost, iVector indices, dVector observ) {
////	dVector vals;
////	int n;
////
////	if (!(vals = (dVector)arr_alloc(num->rvdOmCnt + 1, double)))
////		errMsg("allocation", "chgObjwObserv", "vals", 0);
////
////	for (n = 1; n <= num->rvdOmCnt; n++)
////		vals[n] = cost[n] + observ[coord->rvOffset[2] + n];
////
////	if (changeObjx(lp, num->rvdOmCnt, indices + 1, vals + 1)) {
////		errMsg("solver", "chgObjswObserv", "failed to change the cost coefficients in the solver", 0);
////		return 1;
////	}
////
////	mem_free(vals);
////	return 0;
////}//END chgObjwObserv()


//This function gets a oneproblem and builds a model to be solved by gurobi

oneProblem* newSubprob(oneProblem* sp) {

	/* since the basic structure of subproblem is not modified during the course of the algorithm, we just load it onto the solver */

	sp->model = setupProblem(sp->name, sp->mac, sp->mar, sp->objSense,0.0, sp->objx, sp->senx, sp->rhsx, sp->matbeg, sp->matcnt,
		sp->matind, sp->matval, sp->bdl, sp->bdu, sp->ctype, sp->cname, sp->rname);

	



	if (sp->model == NULL)
	{
		errMsg("Problem Setup", "newSubprob", "sp", 0);
		return NULL;
	}

#if 0
	int     status;
	char probName[NAMESIZE];
	sprintf(probName, "newSubprob%d.lp", agent);
	status = writeProblem(scell->sp->lp, probName);
	if (status) {
		errMsg("write problem", "new_subprob", "failed to write subproblems problem to file", 0);
		return NULL;
	}
#endif

	return sp;
}//END new_subprob


