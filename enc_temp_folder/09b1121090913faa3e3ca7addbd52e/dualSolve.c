#include "stochasticQP.h"

oneCut* dualSolve(probType** prob, cellType* cell, stocType* stoch, double* x) {
	double mubBar;
	int subset = 20;

	/* Allocate memory to sigma, sigma structure contains a vector of pixbCType.  pixbCType is the fixed section we need to store for each solution
	, in our quadratic case, it is beta=mu1.Cbar (dvector) and the fixed section of -1/2 lambda Q lambda - xiBar mu1 + [-yunderscore, ybar] [mu2 , mu3] (double) */


	sigmaType* sigma;
	sigma = (sigmaType*)mem_malloc(sizeof(sigmaType));
	sigma->vals = (pixbCType**)arr_alloc(subset, pixbCType*);


	for (int i = 0; i < subset; i++) {
		sigma->vals[i] = (pixbCType*)mem_malloc(sizeof(pixbCType));}

	for (int i = 0; i < subset; i++) {
		sigma->vals[i]->piCar  = (double*)arr_alloc(prob[0]->num->cols + 1, double);
	}


	/* Allocate memory to delta, delta structure contains a vector of pixbCType.  pixbCType is the fixed section we need to store for each solution
	, in our quadratic case, it is beta=mu1.C (dvector)and the fixed section of -1/2 lambda Q lambda - xi mu1 + [-yunderscoreBar, ybarBar] [mu2 , mu3] (double) */
	
    deltaType* delta;
	delta = (deltaType*)mem_malloc(sizeof(deltaType));
	delta->vals = (pixbCType**)arr_alloc( subset, pixbCType*);
	for (int i = 0; i < subset;i++) {
	   delta->vals[i] = (pixbCType*)mem_malloc(sizeof(pixbCType));
    }

   for (int i = 0; i < subset; i++) {
	   delta->vals[i]->piCar = (double*)arr_alloc(  prob[0]->num->cols +1 , double);
   }

   /* Create lambdatype stracture to store the reduced costs and  pi */
   lambdaType* lambda;
   lambda = (lambdaType*)mem_malloc(sizeof(lambdaType));
   /* alocate memory to lambda structure*/
   lambda->pi = (double**)arr_alloc(subset, double*);
   lambda->mu2 = (double**)arr_alloc(subset , double*);
   lambda->mu3 = (double**)arr_alloc(subset, double*);
   for (int i = 0; i < subset; i++) {
    lambda->pi[i] = (double*)arr_alloc(prob[0]->num->rows + 1, double);
	lambda->mu2[i] = (double*)arr_alloc(prob[0]->num->cols + 1, double);
	lambda->mu3[i] = (double*)arr_alloc(prob[0]->num->cols + 1, double);
	}



	 /* 1. Create a new cut */
	 oneCut* cut = newCut(prob[1]->num->cols);

	/*2. Create the subset of observations for which you want to solve and put the others in another set*/
	int* solveSet;
	int* notsolve;
	int j = 0;
	solveSet = (int*)arr_alloc(subset,int);
	notsolve = (int*)arr_alloc(cell->omega->cnt - subset, int);


	for (int rand = 0; rand < subset; rand++) {
		solveSet[rand] = randInteger(cell->omega->cnt); /*the set we want to solve*/
	}

	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		for (int i = 0; i < subset; i++) {
			if (solveSet[i] == obs) {
				notsolve[j] = obs;  /*the set we do not want to solve*/
				j++;
			}
		}
	}




	sparseVector* bOmega;
	sparseMatrix* COmega;
	sparseVector* yuOmega;

	bOmega = (sparseVector*)mem_malloc(sizeof(sparseVector)); /* The random section of rhs */
	yuOmega = (sparseVector*)mem_malloc(sizeof(sparseVector)); /* The random section of y upperbounds */
	COmega = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix)); /* The random section of C */

	double fixedAlpha;
	double alpha,alpha1,Bx;
	double lql;


	/* 2. loop through subset solveset and solve the subpeoblem */

	for (int obs = 0; obs < subset; obs++) {
		double muBar;

		/* 2a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
		if (solveSubprobdual(prob, cell->subprob, cell->candidX, cell->omega->vals[solveSet[obs]], lambda->pi[obs], &mubBar, lambda->mu2[obs], lambda->mu3[obs])) {
			errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			goto TERMINATE;
		}


		/*obtain the random sections associated with obs*/
		bOmega->cnt = prob[1]->num->rvbOmCnt; 
		bOmega->col = prob[1]->coord->rvbOmRows;
		bOmega->val = cell->omega->vals[solveSet[obs]] + prob[1]->coord->rvOffset[0];

		COmega->cnt = prob[1]->num->rvCOmCnt; COmega->col = prob[1]->coord->rvCOmCols;
		COmega->row = prob[1]->coord->rvCOmRows;
		COmega->val = cell->omega->vals[solveSet[obs]] + prob[1]->coord->rvOffset[1];

		yuOmega->cnt = prob[1]->num->rvyuOmCnt;
		yuOmega->col = prob[1]->coord->rvyuOmRows;
		yuOmega->val = cell->omega->vals[solveSet[obs]] + prob[1]->coord->rvOffset[2];


		/* put pi.Cbar in sigma (fixed part of beta which does not change when C is deterministic) */
		sigma->vals[obs]->piCar = vxMSparse(lambda->pi[obs]+1,

		prob[1]->Cbar, prob[1]->num->rows);     /*check if pi is from 0,lambda->pi[obs]+1?????*/

		/* finding fixed part of alpha equal to -1/2 lambda Q lambda - xiBar .pi - yl.mu3 + yuBar.mu2 */

		/* first  find -1/2 lambda Q lambda = obj -(beta x  - xi(obs) .pi - yl.mu3 + yu(obs).mu2)*/
		/*- xi(obs) .pi - yl.mu3 + yu(obs).mu2*/
		fixedAlpha = -vXv(prob[1]->bBar, lambda->pi[obs] + 1, NULL, prob[1]->num->rows)
			- vXv(prob[1]->ylbar, lambda->mu3[obs], NULL, prob[1]->num->cols) + vXv(prob[1]->yubar, lambda->mu2[obs], NULL, prob[1]->num->cols);

		alpha1 = fixedAlpha -vXvSparse( lambda->pi[obs] + 1, bOmega , prob[1]->num->rows)	+ vXvSparse( lambda->mu2[obs]+1, yuOmega, prob[1]->num->cols);
		Bx = vXv(sigma->vals[obs]->piCar, cell->candidX + 1, NULL, prob[0]->num->cols);
		lql = getObjective(cell->subprob->model) - Bx - alpha1;

		alpha = alpha1 + lql;
		sigma->vals[obs]->fixed = alpha + Bx ;
		sigma->vals[obs]->interceptBar = alpha;


		}
	
	double estimat;
	for (int rand = 0; rand < cell->omega->cnt - subset; rand++) {
		/*obtain the random sections associated with obs*/
		bOmega->cnt = prob[1]->num->rvbOmCnt;
		bOmega->col = prob[1]->coord->rvbOmRows;
		bOmega->val = cell->omega->vals[notsolve[rand]] + prob[1]->coord->rvOffset[0];

		COmega->cnt = prob[1]->num->rvCOmCnt; COmega->col = prob[1]->coord->rvCOmCols;
		COmega->row = prob[1]->coord->rvCOmRows;
		COmega->val = cell->omega->vals[notsolve[rand]] + prob[1]->coord->rvOffset[1];

		yuOmega->cnt = prob[1]->num->rvyuOmCnt;
		yuOmega->col = prob[1]->coord->rvyuOmRows;
		yuOmega->val = cell->omega->vals[notsolve[rand]] + prob[1]->coord->rvOffset[2];
		double max =- 10 ^ 10; /* set it to a large number */
		double violate;
	/*calculate -delta xi.pi + delta ybar.m2 for a given observation*/

		for (int obs = 0; obs < subset; obs++) {
			max = -10 ^ 10;
			violate = -vXvSparse(lambda->pi[obs] + 1, bOmega, prob[1]->num->rows) + vXvSparse(lambda->mu2[obs] + 1, yuOmega, prob[1]->num->cols);
			estimat = sigma->vals[obs]->fixed + violate;
			if (estimat > max) {
				delta->vals[rand]->interceptBar = sigma->vals[obs]->interceptBar+ violate;
				delta->vals[rand]->piCar = sigma->vals[obs]->piCar;
			}		}
		

	/* calculate maximum delta alpha for the observation and keep it as an inexact solution*/


    /* build the cut calculating linear coefs and return it */


	}
	/*solve the subproblem inexactly for the rest of observations in notsolve*/
	


TERMINATE:
	return NULL;
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
