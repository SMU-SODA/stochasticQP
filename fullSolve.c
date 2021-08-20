#include "stochasticQP.h"
extern configType config;

int fullSolve(cellType* cell, omegaType* omega, stocType* stoch, sparseMatrix* C, double* x, int iternum) {

	modelPtr* subModel;
	double* CX;
	CX = (double*)mem_malloc(cell->subprob->mar, sizeof(double));
	double** dualMatrix;
	double* pi;
	dualMatrix = (dVector*)mem_malloc(omega->cnt, sizeof(double*));
	
	/*Calculate the multiplication of C and X*/

	for (int i = 0; i < cell->subprob->mar; i++) {

		CX[i] = 0;
	}

	for (int i = 0; i < C->cnt; i++) {

		CX[C->row[i]] = CX[C->row[i]] + C->val[i + 1] * x[C->col[i]];

	}

	/*loop through observations and solve them all*/
	int* vind = (int*)mem_malloc(omega->numRV, sizeof(int));
	double* val = (double*)mem_malloc(omega->numRV, sizeof(double));

	for (int i = 0; i < omega->cnt; i++) {

		/*change the rhs of subproblem*/
		int k = 0;

		for (int j = 0; j < omega->numRV; j++) {

			vind[j] = stoch->row[j] - cell->master->mar;
			k = stoch->row[j] - cell->master->mar;
			val[j] = omega->vals[i][j + 1] + stoch->mean[j] ;
		}
		if (changeRHSArray(cell->subprob->model, omega->numRV, vind, val)) {
			errMsg("solver", "fullsolve", "failed to change the right-hand side in the solver", 0);
			return 1;
		}
		

#if 1
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


		pi = (dVector)mem_malloc(cell->subprob->marsz, sizeof(double));
		
		
		
		if (getDual(cell->subprob->model, pi, 0, cell->subprob->mar)) {
			errMsg("solver", "fullsolve", "failed to getdual", 0);
			return 1;
		}
		dualMatrix[i] = pi;
	}

	/*Calculate alpha and beta for each observation*/

	double* beta;
	double** alpha;
	for (int t = 0; t < omega->cnt; t++) {
		for (int p = 0; p < cell->master->mac; p++) {
			for (int l = 0; l < cell->subprob->mar; ++l) {
				alpha[t][p] = alpha[t][p] + dualMatrix[t][l]; /** C[l][p]*/
			}
		}
		beta[t] = getObjective(subModel);
		for (int p = 0; p < cell->master->mac; p++) {
			beta[t] = beta[t] - alpha[t][p] * x[p];
		}
	}

	/*construct the cuts and add it to the master*/

	double* newalpha;
	double  newbeta;
	for (int m = 0; m < cell->master->mac; m++) {
		for (int n = 0; n < omega->cnt; n++) {
			newalpha[m] = newalpha[m] + omega->probs[n] * alpha[n][m];
		}
	}
	for (int n = 0; n < omega->cnt; n++) {
		newbeta = newbeta + omega->probs[n] * beta[n];
	}


	/* creat a cut*/

	oneCut* newcut;
	newcut->ck = iternum;

	newcut->alpha = newbeta;
	newcut->beta = newalpha;

	newcut->rowNum = cell->cuts->cnt + 1;				/* row number for master problem in solver */

	cell->cuts->cnt = cell->cuts->cnt + 1;
	cell->cuts->vals[cell->cuts->cnt] = newcut;


	return 0;
}