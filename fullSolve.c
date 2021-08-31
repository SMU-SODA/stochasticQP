#include "stochasticQP.h"
extern configType config;

int fullSolve(cellType* cell, omegaType* omega, stocType* stoch, sparseMatrix* C, double* x, int iternum) {

	double* CX; /*The rhs of the second stage*/
	double** dualMatrix; /*A matrix of dual vectors accosiated with observations*/
	double* pi; /*The dual solution of a given o*/
	double* beta; /*Intercet*/
	double** alpha; /*Slope*/
	double* newalpha; /*Intercet*/
	double  newbeta; /*Slope*/
	int* vind = (int*)arr_alloc(omega->numRV, int);
	double* val = (double*)arr_alloc(omega->numRV +1, double);
	int k = 0;

	/*Memory assignment*/
	CX = (double*) arr_alloc(cell->subprob->mar, double); /*Number of rows from the second-stage*/
	dualMatrix = (dVector*)arr_alloc(omega->cnt, double*); /*A natrix of dual vectors for each observation*/
	beta = (dVector)arr_alloc(omega->cnt, double); /*A vector of intersepts of the cuts associated with observations*/
	alpha = (double**)arr_alloc(omega->cnt, double*);  /*A matrix of cut coefficints associated with observations*/
	newalpha = (double*)arr_alloc(cell->master->mac, double); /*The coefficient of the produced cut*/
	pi = (dVector)arr_alloc(cell->subprob->mar, double); /*A vector dual*/


	/* 1, Create a new cut */


	/* 2. Compute the right-hand side to solve subproblems */
	/* Calculate the multiplication of C and X*/
	for (int i = 0; i < C->cnt; i++) {
		CX[C->row[i]] = CX[C->row[i]] + C->val[i + 1] * x[C->col[i]];
	}

	/* loop through observations and solve subporblem for all of them. */
	for (int i = 0; i < omega->cnt; i++) {

		/* change the rhs of subproblem */
		for (int j = 0; j < omega->numRV; j++) {
			vind[j] = stoch->row[j] - cell->master->mar; 					/* save the index of stochastic constraints in subproblem*/
			val[j] = omega->vals[i][j + 1] + stoch->mean[j] - CX[vind[j]]; 	/* The values associated with the stochastic rhs*/ }

		if (changeRHSArray(cell->subprob->model, omega->numRV, vind, val)) {
			errMsg("solver", "fullsolve", "failed to change the right-hand side in the solver", 0);
			return 1;}

#if defined(WRITE_FILES)
		char str[BLOCKSIZE];
		sprintf(str, "subproblem.lp");
		if (writeProblem(cell->subprob->model, str)) {
			sprintf(str, "failed to write the stage problem %s after modification", cell->subprob->name);
			errMsg("solver", "fullsolve", str, 0);
			return 1;
		}
#endif

		if (solveProblem(cell->subprob->model)) {
			errMsg("solver", "fullsolve", "failed to solve the problem", 0);
			return 1;
		}

		if (getDual(cell->subprob->model, pi, 0, cell->subprob->mar )) {
			errMsg("solver", "fullsolve", "failed to getdual", 0);
			return 1;
		}

		dualMatrix[i] = pi;

		beta[i] = getObjective(cell->subprob->model);
	}

	/*Calculate alpha and beta for each observation*/
	for (int t = 0; t < omega->cnt; t++) {
		alpha[t] = (dVector)arr_alloc(cell->master->mac, double);
		for (int l = 0; l < cell->master->mac; ++l) {
			alpha[t][l] = 0;
		}
	}

	for (int t = 0; t < omega->cnt; t++) {

		for (int l = 0; l < C->cnt; ++l) {
			alpha[t][C->col[l]] =  alpha[t][C->col[l]] + dualMatrix[t][C->row[l]] * C->val[l+1]; /** dual*C **/
		}

		for (int p = 0; p < cell->master->mac; p++) {
			beta[t] = beta[t] - alpha[t][p] * x[p]; /*extract the alpha.x from obj to get beta*/
		}
	}

	/*construct the cuts and add it to the master*/

	for (int m = 0; m < cell->master->mac; m++) {
		for (int n = 0; n < omega->cnt; n++) {
			newalpha[m] = newalpha[m] + omega->probs[n] * alpha[n][m];
		}
	}
	newalpha[cell->master->mac] = -1;
	newbeta = 0;
	for (int n = 0; n < omega->cnt; n++) {
		newbeta = newbeta + omega->probs[n] * beta[n];
	}


	/* create a cut*/

	oneCut* newcut;
	newcut = (oneCut*)mem_malloc(sizeof(oneCut));
	newcut->ck = iternum;

	newcut->alpha = newbeta;
	newcut->beta = newalpha;
	newcut->rowNum = cell->cuts->cnt + 1;				/* row number for master problem in solver */
	cell->cuts->cnt = cell->cuts->cnt + 1;
	cell->cuts->vals[cell->cuts->cnt] = newcut;

	/*update the master problem*/
	dVector matvals;
	int* rmatind;
	rmatind = (int*)arr_alloc(cell->master->mac, double);
	matvals = (double*)arr_alloc(cell->master->mac ,double);
	int j=0;

	for (int i = 0; i< cell->master->mac ;++i) {
		if (newalpha[i] != 0) {
			matvals[j] = newalpha[i];
			rmatind[j] = i;
			j = j + 1;
		}
	}


	/* Add a new linear constraint to a model. */
	if(addRow(cell->master->model, j, -newcut->alpha, GRB_LESS_EQUAL, rmatind, matvals, "cons"))
	{
		errMsg("solver", "fullsolve", "failed to addrow", 0);
		return 1;
	}

#if 1
	char str[BLOCKSIZE];
	sprintf(str, "newmaster.lp");
	if (writeProblem(cell->master->model, str)) {
		sprintf(str, "failed to write the stage problem %s after adding the cut", cell->subprob->name);
		errMsg("solver", "fullsolve", str, 0);
		return 1;
	}
#endif

	solveProblem(cell->master);
	getPrimal(cell->master,x,0, cell->master->mac-1);
	return 0;
}
