/*
 * algorithm.c
 *  Created on: Jul , 2021
 *      Author: Niloofar Fadavi, Harsha Gangammanavar
 */

#include "stochasticQP.h"
extern configType config;

/* Building the cell for the 2-SQP algorithms */
cellType* buildCell(probType** prob , stocType* stoc) {

	cellType* cell;
	cell = (cellType*)mem_malloc(sizeof(cellType));
	cell->k = 0;
	cell->LPcnt = 0;
	cell->optFlag = false;

	/* 1. construct the master problem */
	cell->master = newMaster(prob[0]->sp);
	if ( cell->master == NULL ) {
		errMsg("setup", "buildCell", "failed to setup master", 0);
		return NULL;
	}

	/* 2. construct the subproblem */
	cell->subprob = newSubproblem(prob[1]->sp);
	if ( cell->subprob == NULL ) {
		errMsg("setup", "buildCell", "failed to setup subproblem", 0);
		return NULL;
	}

	/* 3. construct the omega structure */
	cell->omega = newOmega(stoc);

	/* 4. construct the cuts structure */
	cell->maxCuts = config.MAX_ITER;

	cell->cuts = (cutsType*)mem_malloc(sizeof(cutsType));
	cell->cuts->vals = (oneCut**) arr_alloc(cell->maxCuts, oneCut*);
	cell->cuts->cnt = 0;

	cell->fCuts = (cutsType*)mem_malloc(sizeof(cutsType));
	cell->fCuts->vals = (oneCut**) arr_alloc(cell->maxCuts, oneCut*);
	cell->fCuts->cnt = 0;

	/* 5. Allocate memory to candidate and incumbent solution */
	cell->incumbX = (dVector) duplicVector(prob[0]->meanX, prob[0]->num->cols);
	cell->candidX = (dVector) duplicVector(prob[0]->meanX, prob[0]->num->cols);
	cell->incumbEst = 0.0;
	cell->candidEst = 0.0;

	return cell;
}//END buildCell()

int solveSubprob(probType *prob, oneProblem *subproblem, dVector Xvect, dVector obsVals, dVector piS, double *mubBar) {
	dVector rhs = NULL, cost = NULL;
	iVector	indices;

	indices = (iVector) arr_alloc(maximum(prob->num->rows, prob->num->cols), int);
	for ( int n = 0; n < maximum(prob->num->rows,prob->num->cols); n++ )
		indices[n] = n;

	/* (a) compute the right-hand side using current observation and first-stage solution */
	rhs = computeRHS(prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, obsVals);
	if ( rhs == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
		return 1;
	}

	/* (b) change the right-hand side in the solver */
	if ( changeRHSArray(subproblem->model, prob->num->rows, indices, rhs + 1) ) {
		errMsg("solver", "solve_subprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	if ( prob->num->rvdOmCnt > 0 ) {
		/* (c) compute the cost coefficients using current observation */
		cost = computeCostCoeff(prob->num, prob->coord, prob->dBar, obsVals);
		if ( cost == NULL ) {
			errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
			return 1;
		}

		/* (d) change cost coefficients in the solver */
		if ( changeObjCoeffArray(subproblem->model, prob->num->cols, indices, cost+1) ) {
			errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver",0);
			return 1;
		}
	}

	/* (e) Solve the subproblem to obtain the optimal dual solution. */
	if ( solveProblem(subproblem->model) ) {
		errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
		return 1;
	}

#if defined(WRITE_FILES)
	writeProblem(subproblem->model, "cellSubprob.lp");
#endif

#if defined(STOCH_CHECK)
	double obj;
	obj = getObjective(subproblem->model);
	printf("\t\tObjective value of Subproblem  = %lf;\t", obj);
#endif

	/* Record the dual and reduced cost on bounds. */
	if ( getDual(subproblem->model, piS, 0, prob->num->rows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return 1;
	}

	if ( computeMU(subproblem->model, prob->num->cols, mubBar) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to compute mubBar for subproblem", 0);
		return 1;
	}

	if (rhs) mem_free(rhs);
	if (cost) mem_free(cost);
	if (indices) mem_free(indices);
	return 0;
}// END solveSubprob()

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
dVector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, dVector X, dVector observ) {
	int cnt;
	dVector rhs;
	sparseVector bomega;
	sparseMatrix Comega;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows; bomega.val = coord->rvOffset[0] + observ;
	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols; Comega.row = coord->rvCOmRows; Comega.val = coord->rvOffset[1] + observ;

	/* Start with the values of b(omega) -- both fixed and varying */
	rhs = expandVector(bBar->val, bBar->col, bBar->cnt, num->rows);
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
	rhs = MSparsexvSub(Cbar, X, rhs);
	rhs = MSparsexvSub(&Comega, X, rhs);

	return rhs;
}//END computeRHS()

dVector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, dVector observ) {
	dVector cost;
	sparseVector cOmega;
	int	cnt;

	cOmega.cnt = num->rvdOmCnt; cOmega.col = coord->rvdOmCols; cOmega.val = coord->rvOffset[2] + observ;
	cost = expandVector(dBar->val, dBar->col, dBar->cnt, num->cols);
	for (cnt = 1; cnt <= cOmega.cnt; cnt++)
		cost[cOmega.col[cnt]] += cOmega.val[cnt];

	return cost;
}//END computeCostCoeff()


/* This function compute the reduced cost of every second stage variables. They will be used to calculate the \mu x b and then added to the \pi x b. */
int computeMU(modelPtr *model, int numCols, double *mubBar) {
	dVector	dj, u;
	int		n;

	(*mubBar) = 0.0;

	dj = (dVector) arr_alloc(numCols+1, double);
	u = (dVector) arr_alloc(numCols+1, double);

	if ( getPrimal(model, u, 0, numCols) ) {
		errMsg("solver", "computeMU", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlack(model, dj, 0, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		(*mubBar) += dj[n]*u[n];
	}

	mem_free(u); mem_free(dj);
	return 0;
}//END compute_mu()

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a dVector to hold an array of
 * observation and the probability associated with it. */
omegaType* newOmega(stocType* stoc) {
	omegaType* omega;
	int cnt, i, base, idx;

	omega = (omegaType*)mem_malloc(sizeof(omegaType));
	omega->probs = (dVector) arr_alloc(config.MAX_OBS, double);
	omega->weights = (iVector) arr_alloc(config.MAX_OBS, int);
	omega->vals = (dVector*) arr_alloc(config.MAX_OBS, dVector);
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
				}}}}

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

oneProblem* newSubproblem(oneProblem* probSP) {

	oneProblem *stage1;
	stage1 = (oneProblem*)mem_malloc(sizeof(oneProblem));
	stage1->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
	stage1->objQ->col = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
	stage1->objQ->row = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
	stage1->objQ->val = (dVector)arr_alloc(probSP->mac * probSP->mac, double);
	stage1->model = NULL;
	stage1->name = (cString)arr_alloc(NAMESIZE, char);
	stage1->objname = (cString)arr_alloc(NAMESIZE, char);
	stage1->objx = (dVector)arr_alloc(probSP->macsz, double);
	stage1->bdl = (dVector)arr_alloc(probSP->macsz, double);
	stage1->bdu = (dVector)arr_alloc(probSP->macsz, double);
	stage1->ctype = (cString)arr_alloc(probSP->macsz, char);
	stage1->rhsx = (dVector)arr_alloc(probSP->marsz, double);
	stage1->senx = (cString)arr_alloc(probSP->marsz, char);
	stage1->matbeg = (iVector)arr_alloc(probSP->macsz, int);
	stage1->matcnt = (iVector)arr_alloc(probSP->macsz, int);
	stage1->matval = (dVector)arr_alloc(probSP->matsz, double);
	stage1->matind = (iVector)arr_alloc(probSP->matsz, int);
	stage1->cname = (cString*)arr_alloc(probSP->macsz, cString);
	stage1->rname = (cString*)arr_alloc(probSP->marsz, cString);
	stage1->mac = probSP->mac;
	stage1->type= probSP->type;			/* type of problem: LP, QP, MIP or MIQP */

	/* Copy the quadratic part of the objective */
	stage1->objQ->cnt = probSP->objQ->cnt;
	for (int i = 0; i < probSP->objQ->cnt; i++) {
		stage1->objQ->val[i] = probSP->objQ->val[i];
		stage1->objQ->col[i] = probSP->objQ->col[i];
	}

	stage1->objSense = probSP->objSense;
	stage1->mac = probSP->mac;				/* number of columns */
	stage1->mar = probSP->mar;				/* number of rows */
	stage1->numBin = probSP->numBin;		/* number of binary variables in the problem */
	stage1->numInt = probSP->numInt;		/* number of integer variables  (includes both binary and general integer variables) in the problem */
	stage1->numnz = probSP->numnz;			/* number of non-zero elements in constraint matrix */
	stage1->macsz = probSP->macsz;			/* number of columns */
	stage1->marsz = probSP->marsz;			/* number of rows */
	stage1->matsz = probSP->matsz;			/* number of rows */


	strcpy(stage1->objname, probSP->objname);
	strcpy(stage1->name, probSP->name);

	/* Loop over the variables to copy relevant sections */
	for (int i = 0; i < probSP->macsz; i++) {
		stage1->objx[i] = probSP->objx[i];
		stage1->bdl[i] = probSP->bdl[i];
		stage1->bdu[i] = probSP->bdu[i];
		stage1->matbeg[i] = probSP->matbeg[i];
		stage1->matcnt[i] = probSP->matcnt[i];
		stage1->ctype[i] = probSP->ctype[i];
		stage1->cname[i] = (cString)arr_alloc(NAMESIZE, char);
		strcpy(stage1->cname[i], probSP->cname[i]);
	}

	/* Loop over the constraints to copy relevant sections */
	for (int i = 0; i < probSP->marsz; i++) {
		stage1->rhsx[i] = probSP->rhsx[i];
		stage1->senx[i] = probSP->senx[i];
		stage1->rname[i] = (cString)arr_alloc(NAMESIZE, char);
		strcpy(stage1->rname[i], probSP->rname[i]);
	}

	/* Loop over the non-zero elements of constraint matrix to copy relevant sections */
	for (int i = 0; i < probSP->matsz; i++) {
		stage1->matind[i] = probSP->matind[i];
		stage1->matval[i] = probSP->matval[i];
	}

	/* Load the master problem into the solver */
	stage1->model = setupProblem(stage1->name, stage1->mac, stage1->mar, stage1->objSense, 0.0, stage1->objx, stage1->objQ, stage1->senx, stage1->rhsx, stage1->matbeg,
			stage1->matcnt, stage1->matind, stage1->matval, stage1->bdl, stage1->bdu, stage1->ctype, stage1->cname, stage1->rname);
	if ( stage1->model == NULL ) {
		errMsg("setup", "newSubproblem", "failed to load the master problem onto the solver", 0);
		return NULL;
	}

#if 0
	int     status;
	char probName[NAMESIZE];
	sprintf(probName, "newSubprob.lp");
	status = writeProblem(stage1->model, probName);
	if (status) {
		errMsg("write problem", "newSubproblem", "failed to write subproblems problem to file", 0);
		return NULL;
	}
#endif

	return stage1;
}//END newSubproblem

oneProblem *newMaster(oneProblem *probSP) {
	oneProblem* stage0 = NULL;

	stage0 = (oneProblem*)mem_malloc(sizeof(oneProblem));
	stage0->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
	stage0->objQ->col = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
	stage0->objQ->row = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
	stage0->objQ->val = (dVector)arr_alloc(probSP->mac * probSP->mac, double);
	stage0->model = NULL;
	stage0->name = (cString)arr_alloc(NAMESIZE, char);
	stage0->objname = (cString)arr_alloc(NAMESIZE, char);
	stage0->objx = (dVector)arr_alloc(probSP->macsz + 1, double);
	stage0->bdl = (dVector)arr_alloc(probSP->macsz +1, double);
	stage0->bdu = (dVector)arr_alloc(probSP->macsz +1, double);
	stage0->ctype = (cString)arr_alloc(probSP->macsz +1, char);
	stage0->rhsx = (dVector)arr_alloc(probSP->marsz, double);
	stage0->senx = (cString)arr_alloc(probSP->marsz, char);
	stage0->matbeg = (iVector)arr_alloc(probSP->macsz+1, int);
	stage0->matcnt = (iVector)arr_alloc(probSP->macsz+1, int);
	stage0->matval = (dVector)arr_alloc(probSP->matsz, double);
	stage0->matind = (iVector)arr_alloc(probSP->matsz, int);
	stage0->cname = (cString*)arr_alloc(probSP->macsz+1, cString);
	stage0->rname = (cString*)arr_alloc(probSP->marsz, cString);
	stage0->mac = probSP->mac;
	stage0->type= probSP->type;			/* type of problem: LP, QP, MIP or MIQP */

	/* Copy the quadratic part of the objective */
	stage0->objQ->cnt = probSP->objQ->cnt;
	for (int i = 0; i < probSP->objQ->cnt; i++) {
		stage0->objQ->val[i] = probSP->objQ->val[i];
		stage0->objQ->col[i] = probSP->objQ->col[i];
	}

	stage0->objSense = probSP->objSense;
	stage0->mac = probSP->mac;				/* number of columns */
	stage0->mar = probSP->mar;				/* number of rows */
	stage0->numBin = probSP->numBin;		/* number of binary variables in the problem */
	stage0->numInt = probSP->numInt;		/* number of integer variables  (includes both binary and general integer variables) in the problem */
	stage0->numnz = probSP->numnz;			/* number of non-zero elements in constraint matrix */
	stage0->macsz = probSP->macsz;			/* number of columns */
	stage0->marsz = probSP->marsz;			/* number of rows */

	strcpy(stage0->objname, probSP->objname);
	strcpy(stage0->name, probSP->name);

	/* Loop over the variables to copy relevant sections */
	for (int i = 0; i < probSP->macsz; i++) {
		stage0->objx[i] = probSP->objx[i];
		stage0->bdl[i] = probSP->bdl[i];
		stage0->bdu[i] = probSP->bdu[i];
		stage0->matbeg[i] = probSP->matbeg[i];
		stage0->matcnt[i] = probSP->matcnt[i];
		stage0->ctype[i] = probSP->ctype[i];
		stage0->cname[i] = (cString)arr_alloc(NAMESIZE, char);
		strcpy(stage0->cname[i], probSP->cname[i]);
	}

	/* Loop over the constraints to copy relevant sections */
	for (int i = 0; i < probSP->marsz; i++) {
		stage0->rhsx[i] = probSP->rhsx[i];
		stage0->senx[i] = probSP->senx[i];
		stage0->rname[i] = (cString)arr_alloc(NAMESIZE, char);
		strcpy(stage0->rname[i], probSP->rname[i]);
	}

	/* Loop over the non-zero elements of constraint matrix to copy relevant sections */
	for (int i = 0; i < probSP->matsz; i++) {
		stage0->matind[i] = probSP->matind[i];
		stage0->matval[i] = probSP->matval[i];
	}

	/* TODO: add an additional column for the approximation */
	stage0->macsz = probSP->macsz + 1;
	stage0->mac = probSP->mac + 1;
	stage0->objx[stage0->macsz -1] = 1;
	stage0->bdl[stage0->macsz - 1] = -GRB_INFINITY;
	stage0->bdu[stage0->macsz - 1] = GRB_INFINITY;
	stage0->matbeg[stage0->macsz - 1] = probSP->matsz;			
	stage0->matcnt[stage0->macsz - 1] =0;
	stage0->ctype[stage0->macsz -1] = 'C';
	stage0->cname[stage0->macsz -1] = (cString)arr_alloc(NAMESIZE, char);
	strcpy(stage0->cname[stage0->macsz -1], "eta");

	/* Load the master problem into the solver */
	stage0->model = setupProblem(stage0->name, stage0->mac, stage0->mar, stage0->objSense, 0.0, stage0->objx, stage0->objQ, stage0->senx, stage0->rhsx, stage0->matbeg,
			stage0->matcnt, stage0->matind, stage0->matval, stage0->bdl, stage0->bdu, stage0->ctype, stage0->cname, stage0->rname);
	if ( stage0->model == NULL ) {
		errMsg("setup", "newMaster", "failed to load the master problem onto the solver", 0);
		return NULL;
	}

#if defined(WRITE_FILES)
	int     status;
	char probName[NAMESIZE];
	sprintf(probName, "newMaster.lp");
	status = writeProblem(stage0->model, probName);
	if (status) {
		errMsg("write problem", "newMaster", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return stage0;
}//END newMaster()
