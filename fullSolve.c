#include "stochasticQP.h"
extern configType config;






oneCut *fullSolve(probType **prob, cellType* cell, stocType* stoch, double* x) {

	double* CX; /*The rhs of the second stage*/
	double** dualMatrix; /*A matrix of dual vectors accosiated with observations*/
	double* pi; /*The dual solution of a given o*/
	double* beta; /*Intercet*/
	double** alpha; /*Slope*/
	double* newalpha; /*Intercet*/
	double  newbeta; /*Slope*/
	int* vind = (int*)arr_alloc(cell->omega->numRV, int);
	double* val = (double*)arr_alloc(cell->omega->numRV +1, double);
	int k = 0;

	/*Memory assignment*/
	CX = (double*) arr_alloc(prob[1]->num->rows+1, double); /*Number of rows from the second-stage*/
	dualMatrix = (dVector*)arr_alloc(cell->omega->cnt, double*); /*A natrix of dual vectors for each observation*/
	beta = (dVector)arr_alloc(cell->omega->cnt, double); /*A vector of intersepts of the cuts associated with observations*/
	alpha = (double**)arr_alloc(cell->omega->cnt, double*);  /*A matrix of cut coefficints associated with observations*/
	newalpha = (double*)arr_alloc(cell->master->mac, double); /*The coefficient of the produced cut*/
	pi = (dVector)arr_alloc(cell->subprob->mar, double); /*A vector dual*/


	/* 1. Create a new cut */
	oneCut *cut = newCut(prob[0]->num->cols);

	/* 2. Compute the right-hand side to solve subproblems */
	/* Calculate the multiplication of C and X*/
	CX= MSparsexvAdd(prob[1]->Cbar,x, CX);


	/* loop through observations and solve subporblem for all of them. */

	for (int i = 0; i < cell->omega->cnt; i++) {

		/* change the rhs of subproblem */
		for (int j = 0; j < cell->omega->numRV; j++) {
			vind[j] = stoch->row[j] - cell->master->mar; 				         	/* save the index of stochastic constraints in subproblem*/
			val[j] = cell->omega->vals[i][j + 1] + stoch->mean[j] - CX[vind[j]+1]; 	/* The values associated with the stochastic rhs*/ }

		if (changeRHSArray(cell->subprob->model, cell->omega->numRV, vind, val)) {
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

		double obj = getObjective(cell->subprob->model);

		beta = (dVector) arr_alloc(prob[0]->num->cols, double);
		for (int l = 0; l < prob[1]->Cbar->cnt; ++l) {
			beta[prob[1]->Cbar->col[l]] +=  pi[prob[1]->Cbar->row[l]] * prob[1]->Cbar->val[l+1]; /** dual*C **/
			cut->beta[prob[1]->Cbar->col[l]] += cell->omega->probs[i]* pi[prob[1]->Cbar->row[l]] * prob[1]->Cbar->val[l + 1];
		}

		for (int p = 0; p < cell->master->mac; p++) {
			cut->alpha += cell->omega->probs[i]*(obj + beta[p]*x[p]); /*extract the alpha.x from obj to get beta*/
		}

		mem_free(beta);
	}

	return cut;
}//END fullSolve()
