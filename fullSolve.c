#include "stochasticQP.h"
extern configType config;






oneCut *fullSolve(probType **prob, cellType* cell, stocType* stoch, double* x) {

	double* rhs=NULL; /*The rhs of the second stage*/
	double* pi; /*The dual solution of a given o*/
	double* beta; /*Intercet*/
	int k = 0;

	/*Memory assignment*/
	pi = (dVector)arr_alloc(cell->subprob->mar, double); /*A vector dual*/
	int* vind = (int*)arr_alloc(cell->omega->numRV, int);

	/* 1. Create a new cut */
	oneCut *cut = newCut(prob[0]->num->cols);

	/* 2. Compute the right-hand side to solve subproblems */
	/* Calculate the multiplication of C and X*/
	/*CX= MSparsexvAdd(prob[1]->Cbar,x, CX);*/


	/* loop through observations and solve subporblem for all of them. */
	for (int j = 0; j < prob[1]->num->rows; j++) {
		vind[j] = j; 				         	/* save the index of stochastic constraints in subproblem*/
	}
	for (int i = 0; i < cell->omega->cnt; i++) {
		
		/* change the rhs of subproblem */

		updateRHSwState(prob[1]->num, prob[1]->coord, prob[1]->bBar, prob[1]->Cbar, cell->candidX,
			cell->omega->vals[i], &rhs);

		
		
			
		if (changeRHSArray(cell->subprob->model, prob[1]->num->rows, vind, rhs +1)) {
			errMsg("solver", "fullsolve", "failed to change the right-hand side in the solver", 0);
			return NULL;
		}

		mem_free(rhs);

#if defined(WRITE_FILES)
		char str[BLOCKSIZE];
		sprintf(str, "subproblem.lp");
		if (writeProblem(cell->subprob->model, str)) {
			sprintf(str, "failed to write the stage problem %s after modification", cell->subprob->name);
			errMsg("solver", "fullsolve", str, 0);
			return NULL;
		}
#endif

		if (solveProblem(cell->subprob->model)) {
			errMsg("solver", "fullsolve", "failed to solve the problem", 0);
			return NULL;
		}

		if (getDual(cell->subprob->model, pi, 0, cell->subprob->mar )) {
			errMsg("solver", "fullsolve", "failed to getdual", 0);
			return NULL;
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
	mem_free(vind);
	mem_free(pi);
	return cut;
}//END fullSolve()


/* This function computes the right hand side of the subproblem, based on a given X dVector and a given observation of omega.
 * It is defined as:
 *             rhs = r(omega) - C(omega) x X
 * and is calculated as:
 *             rhs = (rbar - Cbar x X) + (rOmega - Comega x X)
 *
 * where the "bar" denotes the fixed or mean value, and the "omega" denotes a random variation from this mean. The function
 * allocates an array for the dVector, which must be freed by the customer.  Also, the zeroth position of this rhs dVector is
 * reserved, and the actual values begin at rhs[1].
 \***********************************************************************/
int updateRHSwState(numType* num, coordType* coord, sparseVector* bBar, sparseMatrix* Cbar, dVector X,
	dVector obs, dVector *rhs) {
	(*rhs) = expandVector(bBar->val, bBar->col, bBar->cnt, num->rows);
	sparseMatrix Comega;
	sparseVector bomega;

	bomega.cnt = num->rvbOmCnt;    bomega.col = coord->rvbOmRows; bomega.val = obs + coord->rvOffset[0];

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols + num->rvbOmCnt;
	Comega.row = coord->rvCOmRows + num->rvbOmCnt; Comega.val = obs + coord->rvOffset[1];

	/* Start with the values of b(omega) -- both fixed and varying */

	for (int cnt = 1; cnt <= bomega.cnt; cnt++)
		(*rhs)[bomega.col[cnt]] += bomega.val[cnt];

	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
	(*rhs) = MSparsexvSub(Cbar, X, (*rhs));
	(*rhs) = MSparsexvSub(&Comega, X, (*rhs));

	return 0;
}//END updateRHSwState()