#include "stochasticQP.h"
int stochasticUpdates(probType** prob, cellType* cell, stocType* stoch, lambdaType* lambda,sigmaType* sigma  ,double* x ,int rand) {

	double* pi, * mu2, * mu3;
	double mubBar;
	sparseVector* bOmega;
	sparseMatrix* COmega;
	sparseVector* yuOmega;

	bOmega = (sparseVector*)mem_malloc(sizeof(sparseVector)); /* The random section of rhs */
	yuOmega = (sparseVector*)mem_malloc(sizeof(sparseVector)); /* The random section of y upperbounds */
	COmega = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix)); /* The random section of C */
	mu2 = (double*)arr_alloc(prob[1]->num->cols + 1, double);
	mu3 = (double*)arr_alloc(prob[1]->num->cols + 1, double);
	pi = (double*)arr_alloc(prob[1]->num->rows + 1, double);


	/* 2a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
	if (solveSubprobdual(prob[1], cell->subprob, cell->candidX, cell->omega->vals[rand], pi, &mubBar, mu2, mu3))
	{
		errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
		return 1;
			}
	/*Update lambda*/



	/*obtain the random sections associated with obs*/
	bOmega->cnt = prob[1]->num->rvbOmCnt;
	bOmega->col = prob[1]->coord->rvbOmRows;
	bOmega->val = cell->omega->vals[rand] + prob[1]->coord->rvOffset[0];

	COmega->cnt = prob[1]->num->rvCOmCnt; COmega->col = prob[1]->coord->rvCOmCols;
	COmega->row = prob[1]->coord->rvCOmRows;
	COmega->val = cell->omega->vals[rand] + prob[1]->coord->rvOffset[1];

	yuOmega->cnt = prob[1]->num->rvyuOmCnt;
	yuOmega->col = prob[1]->coord->rvyuOmRows;
	yuOmega->val = cell->omega->vals[rand] + prob[1]->coord->rvOffset[3];

	calcSigma(sigma, cell, prob, pi, mu2, mu3, bOmega, COmega, yuOmega, rand);

	mem_free(mu2);
	mem_free(mu3);
	mem_free(pi);
}

int calcSigma(sigmaType* sigma, cellType* cell  ,probType** prob, dVector pi, dVector mu2, dVector mu3 , sparseVector* bOmega, sparseMatrix* COmega, sparseVector* yuOmega , int obs) {
	double fixedAlpha , alpha1 , lql , alpha, Bx ;
	/*Check if a the sigma structure should be updated */
	
	/* put pi.Cbar in sigma (fixed part of beta which does not change when C is deterministic) */
	sigma->vals[obs]->piCar = vxMSparse(pi, prob[1]->Cbar, prob[1]->num->rows);


	/* finding fixed part of alpha equal to -1/2 lambda Q lambda - xiBar .pi - yl.mu3 + yuBar.mu2 */

	/* first  find -1/2 lambda Q lambda = obj -(beta x  - xi(obs) .pi - yl.mu3 + yu(obs).mu2)*/
	/*- xi(obs) .pi - yl.mu3 + yu(obs).mu2*/

	fixedAlpha = -vXv(prob[1]->bBar, pi + 1, NULL, prob[1]->num->rows)
		- vXv(prob[1]->ylbar, mu3 + 1, NULL, prob[1]->num->cols) + vXv(prob[1]->yubar, mu2 + 1, NULL, prob[1]->num->cols);


	alpha1 = fixedAlpha - vXvSparse(pi + 1, bOmega) + vXvSparse(mu2 + 1, yuOmega);

	Bx = vXv(sigma->vals[obs]->piCar, cell->candidX + 1, NULL, prob[0]->num->cols);

	lql = getObjective(cell->subprob->model) - Bx - alpha1;


	alpha = fixedAlpha + lql;
	sigma->vals[obs]->fixed = alpha + Bx;
	sigma->vals[obs]->interceptBar = alpha;
	return 0;
}