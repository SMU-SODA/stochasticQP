#include "./solverUtilities/utilities.h"
#include "solverUtilities/solver_gurobi.h"
#include "./smpsReader/smps.h"
#include "./smpsReader/prob.h"
#include "stochasticQP.h"
//
//
//
///* In the regularized QP method, we need to change the rhs of x to d. The
// * 		 A * x 			= b
// * 		 eta + beta * x >= alpha
// * Since x = xbar + d, the corresponding changes will therefore be:
// * 		 A * d = b - A * xbar
// * 		 eta + beta * d >= alpha - beta * xbar
// * But as long as the incumbent sulotion does not change, b - A * xbar and alpha - beta * xbar (for the existing cuts) won't change. So we only need
// * to change it when the incumbent changes.
// *
// * On the other hand, in each iteration, a new cut will be added (and/or some cuts may be dropped) and therefore we need to shift the rhs of the
// * added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of in the routine addCut() in cuts.c. We do not need to worry about the shift
// * of rhs for the dropped cuts.
// * This function performs the change of rhs when the incumbent changes, as described above. */
int changeQPrhs(probType *prob, cellType *cell, dVector xk) {
	int 	status = 0, cnt;
	dVector 	rhs;
	iVector 	indices;

	if (!(rhs =(dVector) arr_alloc(prob->num->rows + cell->cuts->cnt+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(iVector) arr_alloc(prob->num->rows + cell->cuts->cnt, int)))
		errMsg("Allocation", "changeRhs", "indices",0);
	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned dVector, the T sparse
	 dVector, and the x dVector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, xk, rhs);

	/*** new rhs = alpha - beta * xbar (benders cuts)***/
	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[prob->num->rows + cnt +1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, xk, NULL, prob->sp->mac);
		indices[prob->num->rows + cnt] = cell->cuts->vals[cnt]->rowNum;

	}

	/* Now we change the right-hand of the master problem. */
	status = changeRHSArray (cell->master->model, prob->num->rows + cell->cuts->cnt, indices, rhs + 1 );

	if (status)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrhs()

///* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
// * (incumbent solution). */
int changeQPbds(modelPtr* modle, int numCols, sparseVector* lbar, sparseVector* ubar, dVector xk) {
	int 	status = 0, cnt;
	dVector	lbounds, ubounds;
	iVector	lindices, uindices;
	char 	*llu, *ulu;
	dVector bdl = expandVector(lbar->val, lbar->col, lbar->cnt, numCols) ;
	dVector bdu = expandVector(ubar->val, ubar->col, ubar->cnt, numCols) ;
	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);


	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);


	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt + 1] - xk[cnt + 1];
		uindices[cnt] = cnt;
	}

	if (changeBDSArray(modle, "UB", numCols , uindices, ubounds )) {
		errMsg("solver", "solve_subprob", "failed to change the upper bounds in the solver", 0);
		return 1;
	}


	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt+1] - xk[cnt + 1];
	lindices[cnt] = cnt;
	}

	if (changeBDSArray(modle, "LB", numCols , lindices, lbounds )) {
		errMsg("solver", "solve_subprob", "failed to change the upper bounds in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices);
	mem_free(ubounds); mem_free(uindices);
	mem_free(bdl); mem_free(bdu);
	return 0;
}//END changeQPbds()


int changeQPcoef(modelPtr* model, probType* prob ,int numCols, sparseVector* dBar, dVector xk) {

	dVector costFull, cost;

	costFull = expandVector(dBar->val, dBar->col, dBar->cnt, numCols);
    dVector temp1;
    dVector temp2;
    int* indices;

    if(prob->sp->objQ != NULL){
	temp1 = vxMSparse(xk, prob->sp->objQ, prob->num->cols);
	temp2 = MSparsexv(xk, prob->sp->objQ, prob->num->cols);}
    else{
    	temp1 = (dVector) arr_alloc(prob->num->cols + 1 , double);
		temp2 =  (dVector)arr_alloc(prob->num->cols + 1 , double  );
    }

	cost = (double*)arr_alloc(prob->num->cols + 1,double);
	indices = (int*)arr_alloc(prob->num->cols ,int);

	for(int i = 1; i <= prob->num->cols; i++){
		cost[i] = temp1[i] + temp2[i] + costFull[i];
		indices[i - 1] = i - 1;
	}

	/* (b2) change cost coefficients in the solver */
	if (changeObjCoeffArray(model, prob->num->cols, indices, cost + 1)) {
		errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver", 0);
		return 1;
	}
	mem_free(temp1); mem_free(temp2);
	mem_free(cost); mem_free(indices); mem_free(costFull);
}


