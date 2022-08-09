/*
 * algorithm.c
 *  Created on: Jul , 2021
 *      Author: Niloofar Fadavi, Harsha Gangammanavar
 */

#include "stochasticQP.h"
#define _CRTDBG_MAP_ALLOC
#include "crtdbg.h"
extern configType config;

/* Building the cell for the 2-SQP algorithms */
cellType* buildCell(probType** prob , stocType* stoc) {
	cellType* cell;
	cell = (cellType*)mem_malloc(sizeof(cellType));
	cell->k = 0;
	cell->LPcnt = 0;
	cell->optFlag = false;
	cell->piM = NULL;
	cell->time = NULL;

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

	/*6. Allocate memory to SigmaType, LambdaType, DeltaType*/
	cell->lambda 	= NULL;
	cell->sigma 	= NULL;
	cell->delta 	= NULL;
	cell->partition = NULL;

	int structSize = config.SAMPLE_FRACTION * cell->omega->cnt * config.MAX_ITER;
	if (config.ALGOTYPE == DUALLBASED) {
		cell->lambda = newLambda(structSize, prob);
		cell->sigma  = newSigma(structSize, prob);
		cell->delta  = newDelta(structSize, prob,cell);
	}
	else if (config.ALGOTYPE == PARTITIONBASED) {
		cell->partition = newPartition(structSize);
		cell->sigma  = newSigma(structSize, prob);
		cell->lambda = newLambda(structSize, prob);
		cell->delta  = newDelta(structSize, prob,cell);
	}

	return cell;
}//END buildCell()


void newDeltaSol(cellType* cell , int partSize , int obsnum ) {

	cell->deltaSol = (deltaSolType*)mem_malloc(sizeof(deltaSolType)) ;
	cell->deltaSol->cnt = 0;
	cell->deltaSol->sol = (solnType***)arr_alloc(partSize, solnType**);
}; //EndnewDeltaSol



void freecut(cutsType * cut){
	if (cut) {

		if (cut->vals)	{
			for (int i = 0; i < cut->cnt; i++) {
				freeonecut(cut->vals[i]);
			}
			mem_free(cut->vals);
		}
		mem_free(cut);
	}
}

void freeonecut(oneCut* cut) {
	if (cut) {
		if (cut->beta)mem_free(cut->beta);
		if (cut->name)mem_free(cut->name);
		mem_free(cut);
	}

};



int solveSubprob(probType* prob, oneProblem* subproblem, dVector Xvect, dVector obsVals,
		sparseVector* bOmega, sparseMatrix* COmega, sparseVector* dOmega, sparseVector* lOmega, sparseVector* uOmega, solnType *dual) {
	dVector rhs = NULL, cost = NULL, bds = NULL;
	iVector	indices;

	indices = (iVector) arr_alloc(maximum(prob->num->rows, prob->num->cols), int);
	for ( int n = 0; n < maximum(prob->num->rows,prob->num->cols); n++ )
		indices[n] = n;

	/* (a1) compute the right-hand side using current observation and first-stage solution */
	rhs = computeRHS(prob->bBar, prob->Cbar, bOmega, COmega, Xvect, prob->num->rows);
	if ( rhs == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
		return 1;
	}

	/* (a2) change the right-hand side in the solver */
	if ( changeRHSArray(subproblem->model, prob->num->rows, indices, rhs + 1) ) {
		errMsg("solver", "solve_subprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	if (prob->num->rvdOmCnt > 0) {
		/* (b1) compute the cost coefficients using current observation */
		cost = computeCostCoeff(prob->dBar, dOmega, prob->num->cols);
		if (cost == NULL) {
			errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
			return 1;
		}

		for ( int n = 0; n < dOmega->cnt; n++ )
			indices[n] = dOmega->col[n+1]-1;

		/* (b2) change cost coefficients in the solver */
		if (changeObjCoeffArray(subproblem->model, dOmega->cnt, indices, cost + 1)) {
			errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver", 0);
			return 1;
		}
	}

	if ( prob->num->rvyuOmCnt > 0 ) {
		/* (c1) compute the upper bounds using current observation */
		bds = computeBDS(prob->uBar, uOmega, prob->num->cols);
		if (bds == NULL) {
			errMsg("algorithm", "solveSubprob", "failed to compute subproblem upper bounds", 0);
			return 1;
		}

		for ( int n = 0; n < uOmega->cnt; n++ )
			indices[n] =uOmega->col[n+1]-1;

		/* (c2) change upper bounds in the solver */
		if (changeBDSArray(subproblem->model, "UB", uOmega->cnt, indices, bds+1)) {
			errMsg("solver", "solve_subprob", "failed to change the upper bounds in the solver", 0);
			return 1;
		}
		mem_free(bds);
	}

	if ( prob->num->rvylOmCnt > 0 ) {
		/* (c1) compute the lower bounds using current observation */
		bds = computeBDS(prob->uBar, lOmega, prob->num->cols);
		if (cost == NULL) {
			errMsg("algorithm", "solveSubprob", "failed to compute subproblem lower bounds", 0);
			return 1;
		}

		for ( int n = 0; n < lOmega->cnt; n++ )
			indices[n] = lOmega->col[n+1]-1;

		/* (c2) change lower bounds in the solver */
		if (changeBDSArray(subproblem->model, "LB", lOmega->cnt, indices, bds+1)) {
			errMsg("solver", "solve_subprob", "failed to change the lower bounds in the solver", 0);
			return 1;
		}
		mem_free(bds);
	}

#if defined(WRITE_FILES)
	writeProblem(subproblem->model, "cellSubprob.lp");
#endif

	/* (e) Solve the subproblem to obtain the optimal dual solution. */
	if ( solveProblem(subproblem->model) ) {
		errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
		return 1;
	}

#if defined(STOCH_CHECK)
	double obj;
	obj = getObjective(subproblem->model);
	//printf("\t\tObjective value of Subproblem  = %lf;\t", obj);
#endif

	/* Record the primal value */
	if ( getPrimal(subproblem->model , dual->y, 0, prob->num->cols) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the primal", 0);
		return 1;
	}

	/* Record the dual and reduced cost on bounds. */
	if ( getDual(subproblem->model, dual->pi, 0, prob->num->rows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return 1;
	}

	if ( getBoundDual(subproblem->model, prob->num->cols, dual->umu, dual->lmu) ) {
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
dVector computeRHS(sparseVector *bBar, sparseMatrix *Cbar, sparseVector* bOmega, sparseMatrix* COmega, dVector X, int numRows) {
	int cnt;
	dVector rhs;

	/* Start with the values of b(omega) -- both fixed and varying */
	rhs = expandVector(bBar->val , bBar->col, bBar->cnt, numRows);
	for (cnt = 1; cnt <= bOmega->cnt; cnt++)
		rhs[bOmega->col[cnt]] += bOmega->val[cnt];


	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
	rhs = MSparsexvSub(Cbar, X, rhs);
	rhs = MSparsexvSub(COmega, X, rhs);

	return rhs;
}//END computeRHS()

dVector computeCostCoeff(sparseVector *dBar, sparseVector* dOmega, int numCols) {
	dVector costFull, cost;

	costFull = expandVector(dBar->val, dBar->col, dBar->cnt, numCols);
	for (int n = 1; n < dOmega->cnt; n++)
		costFull[dOmega->col[n]] += dOmega->val[n];

	cost = reduceVector(costFull, dOmega->col, dOmega->cnt);

	return cost;
}//END computeCostCoeff()

dVector computeBDS(sparseVector* bdsBar, sparseVector* bdsOmega, int numCols) {
	dVector bdsFull, bds;

	bdsFull = expandVector(bdsBar->val, bdsBar->col, bdsBar->cnt, numCols);
	for ( int n = 1; n <= bdsOmega->cnt; n++ ) {
		bdsFull[bdsOmega->col[n]] += bdsOmega->val[n];
	}

	bds = reduceVector(bdsFull, bdsOmega->col, bdsOmega->cnt);
	mem_free(bdsFull);

	return bds;
}//END computeCostCoeff()

/* This function compute the reduced cost of every second stage variables. They will be used to calculate the \mu x b and then added to the \pi x b. */
int getBoundDual(modelPtr *model, int numCols, double* mu_up, double* mu_low) {
	dVector u, dj;

	u = (dVector) arr_alloc(numCols+1, double);
	dj = (dVector) arr_alloc(numCols+1, double);

	if ( getPrimal(model, u, 0, numCols) ) {
		errMsg("solver", "computeMU", "failed to obtain primal solution", 0);
		return 1;
	}

	if (getDualSlack(model, dj, 0, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	for (int i = 1; i < numCols;i++) {
		if (dj[i] <= 0) {
			mu_up[i] = dj[i];
			mu_low[i] = 0.0;
		}
		else
		{
			mu_up[i] = 0.0;
			mu_low[i] = dj[i];
		}
	}

	mem_free(u);
	mem_free(dj);
	return 0;
}//END compute_mu()

oneProblem* newSubproblem(oneProblem* probSP) {

	oneProblem *stage1;
	stage1 = (oneProblem*)mem_malloc(sizeof(oneProblem));

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

	if (probSP->objQ != NULL ) {
		stage1->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
		stage1->objQ->col = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
		stage1->objQ->row = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
		stage1->objQ->val = (dVector)arr_alloc(probSP->mac * probSP->mac, double);

		/* Copy the quadratic part of the objective */
		stage1->objQ->cnt = probSP->objQ->cnt;
		for (int i = 0; i <= probSP->objQ->cnt; i++) {
			stage1->objQ->val[i] = probSP->objQ->val[i];
			stage1->objQ->col[i] = probSP->objQ->col[i];
			stage1->objQ->row[i] = probSP->objQ->row[i];
		}
	}
	else { stage1->objQ = NULL; }
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

#if defined(WRITE_FILES)
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
	if (probSP->objQ != NULL )  {
		stage0->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
		stage0->objQ->col = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
		stage0->objQ->row = (iVector)arr_alloc(probSP->mac * probSP->mac, int);
		stage0->objQ->val = (dVector)arr_alloc(probSP->mac * probSP->mac, double);


		stage0->objQ->cnt = probSP->objQ->cnt;
		for (int i = 0; i < probSP->objQ->cnt; i++) {
			stage0->objQ->val[i] = probSP->objQ->val[i];
			stage0->objQ->col[i] = probSP->objQ->col[i];
			stage0->objQ->row[i] = probSP->objQ->row[i];
		}
	}
	else { stage0->objQ = NULL; }
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


solnType* buildSolnType (numType *num) {
	solnType* dual;

	dual      = (solnType*)mem_malloc(sizeof(solnType));
	dual->y   = (dVector)arr_alloc(num->cols + 1, double);
	dual->pi  = (dVector)arr_alloc(num->rows + 1, double);
	dual->lmu = (dVector)arr_alloc(num->cols + 1, double);
	dual->umu = (dVector)arr_alloc(num->cols + 1, double);

	return dual;
}//END buildSolnType()


void freeSolnType(solnType *soln) {

	if ( soln ) {
		if ( soln->y) mem_free(soln->y);
		if ( soln->pi) mem_free(soln->pi);
		if ( soln->lmu) mem_free(soln->lmu);
		if ( soln->umu) mem_free(soln->umu);
		mem_free(soln);
	}

}//END freeSolnType()


void VsumVsparse(dVector result, dVector v, sparseVector* vs, int len) {

	copyVector(v, result, len);

	for (int i = 0; i < vs->cnt ; i++) {
		result[vs->col[i]] =  vs->val[i];
	}
}

int stocUpdateQP(cellType* cell, probType* prob, solnType* dual, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* uOmega, sparseVector* lOmega) {
	int lambdaIdx = 0;
	bool newLambdaFlag = false;

	/* Update the lambda structure */
	lambdaIdx = addtoLambda(cell->lambda, dual, prob->num->rows, prob->num->cols, &newLambdaFlag);

	/* add to sigma the alpha and beta and to delta the dalpha and dbeta*/
	if ( newLambdaFlag ) {

		addtoSigma(cell, prob, dual);

		for (int obs = 0; obs < cell->omega->cnt; obs++) {
			/* Add a new row to the delta structure for all observations and the latest lambda (lambdaIdx) */
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0]-1;
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1]-1;

			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3]-1;
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4]-1;

			addtoDelta(cell, prob, COmega, bOmega, uOmega, lOmega, obs, lambdaIdx);
		}

	}
	//	else {
	//		AddtoDel(cell, prob, COmega, bOmega, ybar, yund,obs, lambdaIdx - 1); /* num position in lambda*/
	//	}

	return lambdaIdx;
}//END stocUpdateQP()

void freeCellType (cellType *cell) {

	if (cell->master) {
		mem_free(cell->master->name);
		for (int n = 0; n < cell->master->mac; n++) {
			mem_free(cell->master->cname[n]);
		}
		for (int n = 0; n < cell->master->mar; n++) {
			mem_free(cell->master->rname[n]);
		}
		freeOneProblem(cell->master);
	}
	if (cell->subprob) {
		mem_free(cell->subprob->name);
		for (int n = 0; n < cell->subprob->mac; n++) {
			mem_free(cell->subprob->cname[n]);
		}
		for (int n = 0; n < cell->subprob->mar; n++) {
			mem_free(cell->subprob->rname[n]);
		}
		freeOneProblem(cell->subprob);
	}

	if (cell->candidX) mem_free(cell->candidX);
	if (cell->incumbX) mem_free(cell->incumbX);
	if (cell->piM) mem_free(cell->piM);

	freecut(cell->cuts);
	freecut(cell->fCuts);

	freeSigma(cell->sigma);
	freeLambda(cell->lambda);
	freeDelta(cell->delta, cell->omega->cnt);
	freePartition(cell->partition);
	freeOmegaType(cell->omega,false);

	if (cell->time)mem_free(cell->time);

	mem_free(cell);

} /** EndfreeCellType() **/

PartitionType *newPartition(int Partsize) {
	PartitionType *partition;

	partition = (PartitionType*) mem_malloc(sizeof(PartitionType));
	partition->cnt    = 0;
	partition->part   = (int**) arr_alloc(Partsize, int*);
	partition->basnum = (long long int*) arr_alloc(Partsize, long long int);
	// TODO: Allocate memory to individual elements of partition->part as and when you discover new partitions.

	return partition;
}; /*End newPartition*/

void freePartition( PartitionType* partition) {

	for (int i = 0; i <= partition->cnt ; i++) {
		if ((partition->part[i]) = ! NULL) {
			mem_free(partition->part[i]);
		}
	}
	mem_free(partition->part);
	mem_free(partition->basnum);
	mem_free(partition);

}; //EndfreePartition

int addtoPartition(probType* prob, cellType* cell, sparseVector* uOmega, sparseVector* lOmega, solnType* soln,
		bool* flag, int* up , int* inact , int* low , long long int* base ) {
	long long int index = 0;

	dVector lStat, uStat;

	lStat = (dVector) arr_alloc(prob->num->cols+1, double);
	if ( getDoubleAttributeArray (cell->subprob->model, "LB", 0, prob->num->cols, lStat+1) ) {
		errMsg("solver", "AddtoPart", "failed to obtain variable lower bounds", 0);
		goto TERMINATE;
	}

	uStat = (dVector) arr_alloc(prob->num->cols+1, double);
	if ( getDoubleAttributeArray (cell->subprob->model, "UB", 0, prob->num->cols, uStat+1) ) {
		errMsg("solver", "AddtoPart", "failed to obtain variable upper bounds", 0);
		goto TERMINATE;
	}

	iVector part;
		
	part =	(iVector) arr_alloc(prob->num->cols + 1, int);
	/* for current solution check and save the status of bound constraints, if set to lower bound 
	put 1 in vector part, if on upper bound put 2 in vector part, if inactive put 0 in part */
	for (int i = 1; i <= prob->num->cols; i++) {

		if ( fabs(soln->y[i] - lStat[i]) < config.TOLERANCE ) {
			part[i] = 1;
			(*low)++;
		}
		else if(fabs(uStat[i] - soln->y[i]) < config.TOLERANCE) {
			part[i] = 2;
			(*up)++;
		}
		else {
			part[i] = 0;
			(*inact)++;
		}
	}

	/* produce a number using 3  as a basis from the obtained partition*/
	for (int i = 1; i <= prob->num->cols; i++) {
		index = base[i] * part[i]  + index;
	}

	/* Check to see if this partition is obtained before */
	int partIdx = 0;
	while ( partIdx < cell->partition->cnt ) {
		if (cell->partition->basnum[partIdx] == index) {
			break;
		}
		partIdx++;
	}

	if ( partIdx == cell->partition->cnt ) {
		/* New partition discovered */
		(*flag) = true;

		cell->partition->part[partIdx] = (int*) arr_alloc(prob->num->cols + 1, int);
		copyIntvec(part, cell->partition->part[partIdx] , prob->num->cols);
		cell->partition->basnum[partIdx] = index;
		cell->partition->cnt = cell->partition->cnt + 1;
	}

	mem_free(lStat); mem_free(uStat);
	mem_free(part);
	return partIdx;
	
	TERMINATE:
	mem_free(lStat); mem_free(uStat);
	mem_free(part);
	return -1;
}; /*END addtoPartition()*/

void freeSolSet(solutionSetType* SolSet) {
	for (int i = 0; i < SolSet->cnt; i++) {
		if (SolSet->vals[i]) {
			freeSolnType(SolSet->vals[i]);
		}
	}
	if (SolSet->vals) {
		mem_free(SolSet->vals);
	}
	if (SolSet) {
		mem_free(SolSet);
	}
}; /*EndfreeSolset*/

void addtoDeltaP(cellType* cell, solnType* soln, Mat* W, Mat* T, Mat* WT, probType* prob, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* uOmega, sparseVector* lOmega, int obs, int lambdaIdx, int inact, int up , int low) {

	int elm = 0;

	Mat* rhyu = newmat(prob->num->rows + up, 1, 0);
	Mat* deltaPD;
	double* deltaLd;
	dVector deltay = (dVector)arr_alloc(prob->num->cols + 1, double);
	deltaLd = (dVector)arr_alloc(prob->num->rows + 1, double);



	/* 1. Calculate  vector [delta ubarU , delta rho] */

	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[lambdaIdx][i] == 2) {
			for (int j = 1; j <= uOmega->cnt; j++) {

				{
					if (uOmega->col[j] == i) {
						rhyu->entries[elm] = uOmega->val[i];
					}
				}
			}
			elm++;
		}
	}


	for (int i = 1; i <= bOmega->cnt; i++) {

		rhyu->entries[up + bOmega->col[i] - 1 ] = bOmega->val[i];
	}


	/* 2. Calculate the delta [yI,lambda,nuL,muU] */
	deltaPD =  multiply(WT, rhyu);

	/* 3. put pibar in lambda */
	int cnt = 0;
	copyVector(deltaPD->entries + inact, deltaLd + 1, prob->num->rows-1);

	int* index = (iVector)arr_alloc(prob->num->prevCols + 1,int);
	int* index2 = (iVector)arr_alloc(prob->num->cols + 1, int);
	for (int i = 1; i <= prob->num->prevCols; i++) {
		index[i] = i;
	}

	dVector dbeta = (dVector) arr_alloc(prob->num->prevCols + 1,double);

	dVector temp1, temp2, temp3;

	temp1 = vxMSparse(cell->lambda->pi[lambdaIdx], COmega, prob->num->prevCols);
	temp2 = vxMSparse(deltaLd , prob->Cbar, prob->num->rows);
	temp3 = vxMSparse(deltaLd , COmega, prob->num->rows);

	addVectors(dbeta, temp1, index, prob->num->prevCols );
	addVectors(dbeta, temp2, index , prob->num->prevCols );
	addVectors(dbeta, temp3, index, prob->num->prevCols );

	mem_free(temp1); mem_free(temp2); mem_free(temp3);

	/* If this is a new lambda, and we are in the first obsevation we add a new row first */
	if (obs == 0) {
		cell->delta->vals[lambdaIdx] = (pixbCType**) arr_alloc(cell->omega->cnt, pixbCType*);
		cell->delta->cnt++;
	}
	cell->delta->vals[lambdaIdx][obs] = (pixbCType*) mem_malloc(sizeof(pixbCType));

	/* Calculate delta beta in sparce struct */
	/*To do: should we add any conditions on null beta like dual subrutine? */

	cell->delta->vals[lambdaIdx][obs]->beta = reduceVector(dbeta, prob->coord->rvCOmCols, prob->num->rvCOmCnt);

	/* Calculate deltaY */

	dVector deltaU = (dVector)arr_alloc(prob->num->cols + 1, double);
	for (int i = 1; i <= uOmega->cnt; i++) {
		deltaU[uOmega->col[i]] = uOmega->val[i];
	}

	dVector deltaL = (dVector)arr_alloc(prob->num->cols + 1, double);
	for (int i = 1; i <= lOmega->cnt; i++) {
		deltaL[lOmega->col[i]] = lOmega->val[i];
	}

	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[lambdaIdx][i] == 1) {

			deltay[i] = deltaL[i];
		}
		if (cell->partition->part[lambdaIdx][i] == 2 ) {

			deltay[i] = deltaU[i];
		}
		if (cell->partition->part[lambdaIdx][i] == 0) {
			deltay[i] = deltaPD->entries[cnt];
			cnt++;
		}
	}

	/* Calculate deltanu */

	cnt = 0;
	dVector dnu = (dVector)arr_alloc(prob->num->cols + 1, double);
	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[lambdaIdx][i] == 0 || cell->partition->part[lambdaIdx][i] == 2) {

			dnu[i] = 0;

		}

		if (cell->partition->part[lambdaIdx][i] == 1) {
			dnu[i] = cell->lambda->pd[lambdaIdx]->entries[inact + prob->num->rows + cnt];
			cnt++;
		}
	}

	/* 10. Calculate delta mu */

	cnt = 0;
	dVector dmu = (dVector)arr_alloc(prob->num->cols + 1, double);

	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[lambdaIdx][i] == 0 || cell->partition->part[lambdaIdx][i] == 1) {

			dmu[i] = 0;
		}

		if (cell->partition->part[lambdaIdx][i] == 2) {
			dmu[i] = cell->lambda->pd[lambdaIdx]->entries[inact + low + prob->num->rows + cnt];
			cnt++;
		}
	}

	for (int i = 1; i <= prob->num->cols ; i++) {
		index2[i] = i;
	}


	/* calculate alpha */
	dVector temp4, temp5;
	temp4 = vxMSparse(deltay, prob->sp->objQ, prob->num->cols);
	temp5 = vxMSparse(cell->lambda->y[lambdaIdx], prob->sp->objQ, prob->num->cols);

	cell->delta->vals[lambdaIdx][obs]->alpha = -vXv(temp4, deltay, index2, prob->num->cols)
								- 2 * vXv(temp5, deltay, index2, prob->num->cols) + vXvSparse(deltaLd , prob->bBar)
								+ vXvSparse(cell->lambda->pi[lambdaIdx], bOmega)
								+ vXvSparse(deltaLd, bOmega) + vXvSparse(dnu, prob->lBar)
								+ vXvSparse(cell->lambda->lmu[lambdaIdx], lOmega)
								+ vXvSparse(dnu, lOmega) + vXvSparse(dmu, prob->uBar)
								+ vXvSparse(dmu, uOmega)
								+ vXvSparse(cell->lambda->umu[lambdaIdx], uOmega);

	mem_free(temp4); mem_free(temp5);
	mem_free(dbeta);
	mem_free(index);
	mem_free(index2);
	freemat(rhyu);
	freemat(deltaPD);
	mem_free(deltaLd);
	mem_free(deltay);
	mem_free(deltaL);
	mem_free(deltaU);
	mem_free(dmu);
	mem_free(dnu);
}; //EndAddtoDeltaP

void AddtoSigmaP(cellType* cell ,solnType* sol , probType* prob ) {
	int* index;
	int cnt = cell->partition->cnt - 1;
	index = (int*) arr_alloc(prob->num->cols + 1, int);

	for (int i = 1; i <= prob->num->cols; i++) {
		index[i] = i;
	}
	cnt = cell->sigma->cnt;
	cell->sigma->cnt++;
	cell->sigma->vals[cnt] = (pixbCType*) mem_malloc(sizeof(pixbCType));

	dVector temp;
	temp = vxMSparse(cell->lambda->y[cnt], prob->sp->objQ, prob->num->cols);

	cell->sigma->vals[cnt]->alpha = -vXv(temp, cell->lambda->y[cnt], index, prob->num->cols) +
			vXvSparse(cell->lambda->pi[cnt], prob->bBar) + vXvSparse(cell->lambda->lmu[cnt], prob->lBar) + vXvSparse(cell->lambda->umu[cnt], prob->uBar);
	cell->sigma->vals[cnt]->beta = vxMSparse(cell->lambda->pi[cnt], prob->Cbar, prob->num->prevCols);

	mem_free(temp);
	mem_free(index);
}; //EndAddtoSigmaP

void addtoLambdaP(cellType* cell, solnType* soln, Mat* WT ,  probType* prob, sparseVector*  bOmega, sparseVector* uOmega,
		sparseVector* lOmega, int low, int up, int inact) {

	int idx = cell->lambda->cnt;
	cell->lambda->cnt++;
	cell->lambda->y[idx]   = (dVector) arr_alloc(prob->num->cols + 1, double);
	cell->lambda->umu[idx] = (dVector) arr_alloc(prob->num->cols + 1, double);
	cell->lambda->lmu[idx] = (dVector) arr_alloc(prob->num->cols + 1, double);
	cell->lambda->pi[idx]  = (dVector) arr_alloc(prob->num->rows + 1, double);

//	cell->lambda->pd[idx] = (Mat*) mem_malloc(sizeof(Mat));
//	cell->lambda->pd[idx]->entries = (dVector)arr_alloc(prob->num->cols + prob->num->rows, double);  /* stores the vector [yI,PI,nuL,muU] bar */

	Mat* rhyu = newmat(prob->num->rows + up, 1, 0);

	//*** Calculate yBar = y - [w,T]delta ***//

	//** Build a vector [yi,lambda,nul,muu]
	Mat* pdsol = newmat(prob->num->cols + prob->num->rows, 1 , 0); /* The current primal-dual vector [yI,PI,nuL,muU] in this partition */

	/* 1. Copy the inactive variables in yI segment of pdsol*/

	for (int i = 1; i  <= prob->num->cols; i++) {
		int num = 0;
		if (cell->partition->part[idx][i] == 0)
		{
			pdsol->entries[num] = soln->y[i];
			//cell->lambda->y[idx][i]= soln->y[i];
			num++;
		}
	}

	/* 2. Put the duals of equality constraints after yI in pdsol */
	copyVector(soln->pi+1, pdsol->entries + inact, prob->num->rows-1);

	/* 3. copy the nonzero nu and mu values*/
	int num1 = 0;
	int num2 = 0;

	for (int i = prob->num->cols; i >= 1; i--) {

		if (cell->partition->part[idx][i] == 1 )
		{
			pdsol->entries[inact + prob->num->rows + num1] = soln->lmu[i];
			num1++;
		}
		if ( cell->partition->part[idx][i] == 2)
		{
			pdsol->entries[inact + prob->num->rows + low + num2] = soln->umu[i];
			num2++;
		}
	}

	/* 5.Build a vector of [delta yubar , delta rho ] */

	int elm = 0;

	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[idx][i] == 2){
			for (int j = 1; j <= uOmega->cnt; j++) {

				{/*check i-1*/
					if (uOmega->col[j] == i) {  /* col is from 0 or 1?*/
						rhyu->entries[elm] = uOmega->val[i];
					}
				}
			}
			elm++;
		}
	}


	for (int i = 1; i <= bOmega->cnt; i++) {

		/*check i-1*/
		rhyu->entries[up + bOmega->col[i]] = bOmega->val[i];
	}


	/* 6.Calculate ybar = y - [WT]delta */
	Mat *tempMat1, *tempMat2;
	tempMat1 = scalermultiply(WT,-1);
	tempMat2 = multiply(tempMat1, rhyu);

	cell->lambda->pd[idx] = sum(pdsol, tempMat2); /* to do : add minus */
	freemat(tempMat1); freemat(tempMat2);

	/* 7. Complete ybar in lambda */
	dVector Ubar = (dVector)arr_alloc(prob->num->cols + 1 ,double);
	for (int i = 1; i <= prob->uBar->cnt ; i++) {
		Ubar[prob->uBar->col[i]] = prob->uBar->val[i];
	}
	dVector Lbar = (dVector)arr_alloc(prob->num->cols + 1, double);
	for (int i = 1; i <= prob->lBar->cnt; i++) {
		Lbar[prob->lBar->col[i]] = prob->lBar->val[i];
	}
	int cnt = 0;
	for (int i = 1; i <= prob->num->cols; i++) {

		if (cell->partition->part[idx][i] == 2 ) {

			cell->lambda->y[idx][i] = Ubar[i];
		}

		if ( cell->partition->part[idx][i] == 1) {

			cell->lambda->y[idx][i] = Lbar[i];
		}

		if (cell->partition->part[idx][i] == 0) {

			cell->lambda->y[idx][i] = cell->lambda->pd[idx]->entries[cnt];
			cnt++;
		}
	}

	/* 8. Complete pibar in lambda*/
	copyVector(cell->lambda->pd[idx]->entries + inact, cell->lambda->pi[idx] + 1, prob->num->rows-1);

	/* 9. Complete nu in lambda*/

	cnt = 0;
	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[idx][i] == 0 || cell->partition->part[idx][i] == 2) {

			cell->lambda->lmu[idx][i] = 0;

		}
		if (cell->partition->part[idx][i] == 1) {
			cell->lambda->lmu[idx][i] = cell->lambda->pd[idx]->entries[ inact + prob->num->rows +cnt];
			cnt++;
		}
	}

	/* 10. Complete mu in lambda*/

	cnt = 0;
	for (int i = 1; i <= prob->num->cols; i++) {
		if (cell->partition->part[idx][i] == 0 || cell->partition->part[idx][i] == 1) {

			cell->lambda->umu[idx][i] = 0;

		}

		if (cell->partition->part[idx][i] == 2) {
			cell->lambda->umu[idx][i] = cell->lambda->pd[idx]->entries[inact + low + prob->num->rows + cnt];
			cnt++;
		}
	}

	freemat(pdsol);
	freemat(rhyu);
	mem_free(Ubar);
	mem_free(Lbar);

}; //EndaddtoLambdaP


void Buildbase(long long int* basis, int cols, int base) {

	long long int pow = 1;

	for (int i = 1; i <= cols; i++) {

		basis[i] = pow;
		pow = (long long int) pow * base;

	}
};//EndBuildbase

Mat* CombineWT(probType* prob,Mat* W, Mat* T , int low , int up  , int inact) {
	Mat* WT = newmat(prob->num->cols + prob->num->rows, up + prob->num->rows, 0);

	WT->col = up + prob->num->rows;
	WT->row = prob->num->cols + prob->num->rows;
	copyVector(W->entries, WT->entries, (inact + prob->num->rows) * (up + prob->num->rows));
	copyVector(T->entries, WT->entries + (inact + prob->num->rows) * (up + prob->num->rows), (low + up) * (up + prob->num->rows));

	return WT;
};

//* Matrix W calculation  *//
void CalcWT(cellType* cell, probType* prob, sparseMatrix* Q, sparseMatrix* D , Mat** W, Mat** T, int low, int up, int inact) {
	_CrtDumpMemoryLeaks();
	int elm = 0;
	Mat* QII = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QIU = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QLU = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QUU = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QUI = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QLI = transSparsM(Q, prob->num->cols, prob->num->cols);
	showmat(QLI);
	_CrtDumpMemoryLeaks();
	
	Mat* DMU = transSparsM(D, prob->num->cols, prob->num->rows);
	_CrtDumpMemoryLeaks();
	Mat* DML = transSparsM(D, prob->num->cols, prob->num->rows);
	Mat* DMI = transSparsM(D, prob->num->cols, prob->num->rows);
	_CrtDumpMemoryLeaks();
	
	Mat* M1;
	Mat* M2;
	Mat* minvM1;
	Mat* DMLT;
	Mat* DMUT;
	Mat* w2;
	Mat* W2W;
	Mat* DMIT;
	Mat* invM1;
	Mat* w1;
	/*Build Q(II) in a mat strcture*/
	int cnt = cell->partition->cnt - 1;
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2) {
			removeCol(QII, i);
			removeRow(QII, i);
		}
	}

	showmat(DMI);
	/*Build D(MI)*/
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2) {
			removeCol(DMI, i);
			printf("ASCII value = %d, Character = %c\n", i);
			showmat(DMI);
		}
	}
	showmat(DMI);

	/*Build D(MI) transpose*/

	DMIT = transpose(DMI);

	M1 = newmat(prob->num->rows + inact, prob->num->rows + inact, 0);

	//* Build Matrix M1 which is equal to [ QII , DMIT ; DMI , 0 ] *//
	// 1. place the first I rows

	for (int i = 0; i < inact; i++) {
		for (int j = 0; j < inact; j++) {
			M1->entries[elm] = QII->entries[i * inact + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			M1->entries[elm] = DMIT->entries[i * prob->num->rows + j];
			elm++;
		}
	}

	// 2. place the next M rows
	for (int i = 0; i < prob->num->rows; i++) {
		for (int j = 0; j < inact; j++) {
			M1->entries[elm] = DMI->entries[(i)*inact + j];
			elm++;
		}
		for (int j = 1; j <= prob->num->rows; j++) {
			M1->entries[elm] = 0;
			elm++;
		}
	}


	showmat(M1);
	invM1 = inverse(M1);

	/* Build Matrix M2 Which is equal to [QIU , 0 ; DMU , -I] */
	M2 = newmat(prob->num->rows + inact, prob->num->rows + up, 0);

	/*Build QIU in a mat structure*/
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2) {
			removeRow(QIU, i);
		}
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1) {
			removeCol(QIU, i);
		}
	}

	/*Build DMU in a mat strcture*/
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1) {
			removeCol(DMU, i);
		}
	}

	/* first I rows of M2*/

	elm = 0;
	for (int i = 0; i < inact; i++) {

		for (int j = 0; j < up; j++) {
			M2->entries[elm] = QIU->entries[i * up + j];
			elm++;
		}

		for (int j = 0; j < prob->num->rows; j++) {
			M2->entries[elm] = 0;
			elm++;
		}

	}

	/* next M rows of M2*/
	for (int i = 0; i < prob->num->rows; i++) {
		for (int j = 0; j < up; j++) {
			M2->entries[elm] = DMU->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			if (i == j) {
				M2->entries[elm] = -1;
				elm++;
			}
			else {
				M2->entries[elm] = 0;
				elm++;
			}
		}
	}

	/*Calculate 4 components of the W*/

	minvM1 = scalermultiply(invM1, -1);
	(*W) = multiply(minvM1, M2);

	//* Matrix T caculation  *//

	/* QLU */
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 2)
		{
			removeRow(QLU, i);
		}
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removeCol(QLU, i);
		}
	}

	/* QUU */
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removeRow(QUU, i);
			removeCol(QUU, i);
		}
	}

	/* QUI */
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removeRow(QUI, i);
		}
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			removeCol(QUI, i);
		}
	}

	/* QLI */
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 2)
		{
			removeRow(QLI, i);
		}
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			removeCol(QLI, i);
		}
	}

	/*DML*/
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 2) {
			removeCol(DML, i);
		}
	}

	/* Build T=  w1 - w2 * w*/
	w1 = newmat(low + up, prob->num->rows + up, 0);

	/*Build w1*/
	/* first L rows of W1*/

	elm = 0;
	for (int i = 0; i < low; i++) {
		for (int j = 0; j < up; j++) {
			w1->entries[elm] = -QLU->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			w1->entries[elm] = 0;
			elm++;
		}
	}

	/* next U rows of W1*/

	for (int i = 0; i < up; i++) {
		for (int j = 0; j < up; j++) {
			w1->entries[elm] = -QUU->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {

			w1->entries[elm] = 0;
			elm++;

		}
	}

	/*Build w2*/

	w2 = newmat(low + up, prob->num->rows + inact, 0);
	elm = 0;
	DMLT = transpose(DML);
	DMUT = transpose(DMU);

	/* first L rows of W2*/
	for (int i = 0; i < low; i++) {
		for (int j = 0; j < inact; j++) {
			w2->entries[elm] = -QLU->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			w2->entries[elm] = -DMLT->entries[i * prob->num->rows + j];
			elm++;
		}
	}

	/* next U rows of W2*/

	for (int i = 0; i < up; i++) {
		for (int j = 0; j < inact; j++) {
			w2->entries[elm] = -QUI->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			w2->entries[elm] = -DMUT->entries[i * prob->num->rows + j];
			elm++;
		}
	}

	/*Build T*/

	W2W = multiply(w2, (*W));
	(*T) = sum(w1, W2W);

	if (DMU != NULL) {
		freemat(DMU);
	}
	if (DML != NULL) {
		freemat(DML);
	}
	if (DMI != NULL) {
		freemat(DMI);
	}
	if (M1 != NULL) {
		freemat(M1);
	}
	if (M2 != NULL) {
		freemat(M2);
	}
	if (invM1 != NULL) {
		freemat(invM1);
	}
	if (DMLT != NULL) {
		freemat(DMLT);
	}
	if (DMUT != NULL) {
		freemat(DMUT);
	}
	if (w1 != NULL) {
		freemat(w1);
	}
	if (w2 != NULL) {
		freemat(w2);
	}
	if (W2W != NULL) {
		freemat(W2W);
	}
	if (minvM1 != NULL) {
		freemat(minvM1);
	}
	if (DMIT != NULL) {
		freemat(DMIT);
	}
	if (QII != NULL) {
		freemat(QII);
	}
	if (QIU != NULL) {
		freemat(QIU);
	}
	if (QLU != NULL) {
		freemat(QLU);
	}
	if (QUU != NULL) {
		freemat(QUU);
	}
	if (QUI != NULL) {
		freemat(QUI);
	}
	if (QLI != NULL) {
		freemat(QLI);
	}
}
