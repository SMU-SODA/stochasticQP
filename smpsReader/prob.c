/*
 * prob.c
 *
 *  Created on: Sep 30, 2015
 *      Author: Harsha Gangammanavar
 */

#include "prob.h"

/* Decomposes the problem _orig_ into subproblems as well as decomposes the stochastic information _stoc_ into stage stochastic information. The decomposition
 * is carried out using information specified in _tim_. The function also stores stage lower bound information provided in _Lb_. It returns an array of
 * probType structures, each probType corresponds to a particular stage */
probType **newProbwSMPS(cString inputDir, cString probName, stocType **stoc, int *numStages) {
	oneProblem 	*orig = NULL;
	timeType 	*tim = NULL;
	probType 	**prob = NULL; 
	dVector 	meanX = NULL, lb = NULL;
	int		 	i, k, m, t;


	/* Read the SMPS files */
	openSolver();
	if ( readFiles(inputDir, probName, &orig, &tim, stoc) ) {
		errMsg("read", "newProb_SMPS", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}

	(*numStages) = tim->numStages;

	/* allocate memory to elements of probType */
	prob = (probType **) arr_alloc(tim->numStages, probType *); 

	/* allocate memory to members of probType for stagewise subProblems, and allocate values to static fields*/
	for ( t = 0; t < tim->numStages; t++ ) {
		prob[t] = (probType *) mem_malloc(sizeof(probType));
		prob[t]->sp = (oneProblem *) mem_malloc (sizeof(oneProblem));
		prob[t]->sp->model = NULL; prob[t]->lb = 0.0;

		/* Cost coefficient */
		prob[t]->dBar = (sparseVector *) mem_malloc(sizeof(sparseVector));
		prob[t]->dBar->col = (iVector) arr_alloc(orig->mac+1, int);
		prob[t]->dBar->val = (dVector) arr_alloc(orig->mac+1, double);
		prob[t]->dBar->cnt = 0;

		/* Right-hand side */
		prob[t]->bBar = (sparseVector *) mem_malloc(sizeof(sparseVector));
		prob[t]->bBar->col = (iVector) arr_alloc(orig->mar+1, int);
		prob[t]->bBar->val = (dVector) arr_alloc(orig->mar+1, double);
		prob[t]->bBar->cnt = 0;

		/* Variable upper bound */
		prob[t]->uBar = (sparseVector*)mem_malloc(sizeof(sparseVector));
		prob[t]->uBar->col = (iVector)arr_alloc(orig->mac + 1, int);
		prob[t]->uBar->val = (dVector)arr_alloc(orig->mac + 1, double);
		prob[t]->uBar->cnt = 0;

		/* Variable lower bound */
		prob[t]->lBar = (sparseVector*)mem_malloc(sizeof(sparseVector));
		prob[t]->lBar->col = (iVector)arr_alloc(orig->mac + 1, int);
		prob[t]->lBar->val = (dVector)arr_alloc(orig->mac + 1, double);
		prob[t]->lBar->cnt = 0;

		if ( t < tim->numStages - 1 ) {
			prob[t]->sp->mar = prob[t]->sp->marsz = tim->row[t+1] - tim->row[t];
			prob[t]->sp->mac = prob[t]->sp->macsz = tim->col[t+1] - tim->col[t];
		}
		else {
			prob[t]->sp->mar = prob[t]->sp->marsz = orig->mar - tim->row[t];
			prob[t]->sp->mac = prob[t]->sp->macsz = orig->mac - tim->col[t];
		}
		prob[t]->sp->numInt = 0;
		prob[t]->sp->numBin = 0;
		prob[t]->sp->matsz  = 0;
		prob[t]->sp->numnz  = 0;
		prob[t]->sp->objSense = orig->objSense;

		/* stage oneProblem */
		prob[t]->sp->model   = NULL;
		prob[t]->sp->name    = (cString) arr_alloc(NAMESIZE, char);
		prob[t]->sp->objname = (cString) arr_alloc(NAMESIZE, char);
		prob[t]->sp->objx    = (dVector) arr_alloc(prob[t]->sp->macsz, double);
		prob[t]->sp->bdl     = (dVector) arr_alloc(prob[t]->sp->macsz, double);
		prob[t]->sp->bdu     = (dVector) arr_alloc(prob[t]->sp->macsz, double);
		prob[t]->sp->ctype   = (cString) arr_alloc(prob[t]->sp->macsz, char);
		prob[t]->sp->rhsx    = (dVector) arr_alloc(prob[t]->sp->marsz, double);
		prob[t]->sp->senx    = (cString) arr_alloc(prob[t]->sp->marsz, char);
		prob[t]->sp->matbeg  = (iVector) arr_alloc(prob[t]->sp->macsz, int);
		prob[t]->sp->matcnt  = (iVector) arr_alloc(prob[t]->sp->macsz, int);
		prob[t]->sp->matval  = (dVector) arr_alloc(orig->matsz, double);
		prob[t]->sp->matind  = (iVector) arr_alloc(orig->matsz, int);
		prob[t]->sp->cname   = (cString *) arr_alloc(prob[t]->sp->macsz, cString);
		prob[t]->sp->rname   = (cString *) arr_alloc(prob[t]->sp->marsz, cString);

		strcpy(prob[t]->sp->objname, orig->objname);
		sprintf(prob[t]->sp->name, "%s_%d", orig->name, t);

		/* stage transfer matrix */
		if ( t == 0 )
			prob[t]->Cbar = NULL;
		else {
			prob[t]->Cbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix));
			prob[t]->Cbar->row = (iVector) arr_alloc(orig->matsz + 1, int);
			prob[t]->Cbar->col = (iVector) arr_alloc(orig->matsz + 1, int);
			prob[t]->Cbar->val = (dVector) arr_alloc(orig->matsz + 1, double);
			prob[t]->Cbar->cnt = 0;
		}

		/* Allocate memory to all recourse matrices */
		prob[t]->Dbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix));
		prob[t]->Dbar->row = (iVector) arr_alloc(orig->matsz+1, int);
		prob[t]->Dbar->col = (iVector) arr_alloc(orig->matsz+1, int);
		prob[t]->Dbar->val = (dVector) arr_alloc(orig->matsz+1, double);
		prob[t]->Dbar->cnt = 0;

		/* This allocation is done later after confirming that the stage problem has quadratic objective */
		prob[t]->sp->objQ = NULL;

		/* TODO (HG): stage dynamics: include dynamics to the model */
		if ( t == 0) {
			prob[t]->Abar = NULL;
			prob[t]->Bbar = NULL;
			prob[t]->aBar = NULL;
			prob[t]->cBar = NULL;
		}
		else {
			prob[t]->Abar = NULL;
			prob[t]->Bbar = NULL;
			prob[t]->aBar = NULL;
			prob[t]->cBar = NULL;
		}
	}

	t = 0;
	for ( t = 0; t < tim->numStages-1; t++ ) {
		/* lower bound on cost-to-go function */
		prob[t]->lb = 0.0;

		/* copy column information for non-terminal stage */
		for (m = tim->col[t]; m < tim->col[t + 1]; m++) {
			k = m - tim->col[t];


			prob[t]->dBar->val[prob[t]->dBar->cnt + 1] = orig->objx[m];
			prob[t]->dBar->col[prob[t]->dBar->cnt + 1] = m - tim->col[t] + 1;
			prob[t]->dBar->cnt++;


			prob[t]->sp->objx[k] = orig->objx[m];
			prob[t]->sp->bdl[k] = orig->bdl[m];
			prob[t]->sp->bdu[k] = orig->bdu[m];


			if ( orig->ctype[m] == 'I')
				prob[t]->sp->numInt++;
			else if ( orig->ctype[m] == 'B' )
				prob[t]->sp->numBin++;
			prob[t]->sp->ctype[k] = orig->ctype[m];

			prob[t]->sp->cname[k] = (cString) arr_alloc(NAMESIZE, char);
			strcpy(prob[t]->sp->cname[k], orig->cname[m]);

			prob[t]->sp->matcnt[k] = 0;
			for ( i = orig->matbeg[m]; i < orig->matbeg[m]+orig->matcnt[m]; i++ ) {
				if (orig->matind[i] < tim->row[t+1]) {
					/* The coefficient is part of the current stage constraint matrix (recourse)*/
					if ( k == 0 )
						prob[t]->sp->matbeg[k] = 0;
					else
						prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1] + prob[t]->sp->matcnt[k-1];
					prob[t]->sp->matval[prob[t]->sp->matsz] = orig->matval[i];
					prob[t]->sp->matind[prob[t]->sp->matsz] = orig->matind[i] - tim->row[t];
					++prob[t]->sp->matcnt[k];
					++prob[t]->sp->matsz;
					++prob[t]->sp->numnz;
					prob[t]->Dbar->val[prob[t]->Dbar->cnt+1] = orig->matval[i];
					prob[t]->Dbar->col[prob[t]->Dbar->cnt+1] = m-tim->col[t]+1;
					prob[t]->Dbar->row[prob[t]->Dbar->cnt+1] = orig->matind[i] - tim->row[t]+1;
					++prob[t]->Dbar->cnt;
				}
				else {
					/* The coefficient is part of the next stage transfer matrix C matrix, if the row is from the previouse
					stage but the column is from the current one*/
					prob[t+1]->Cbar->val[prob[t+1]->Cbar->cnt+1] = orig->matval[i];
					prob[t+1]->Cbar->row[prob[t+1]->Cbar->cnt+1] = orig->matind[i] - tim->row[t+1]+1;
					prob[t+1]->Cbar->col[prob[t+1]->Cbar->cnt+1] = m-tim->col[t]+1;
					++prob[t+1]->Cbar->cnt;
				}
			}
		}

		/* copy row information for non-terminal stage */
		for ( m = tim->row[t]; m < tim->row[t+1]; m++ ) {
			k = m - tim->row[t];

			prob[t]->sp->rhsx[k] = orig->rhsx[m];
			prob[t]->sp->senx[k] = orig->senx[m];

			prob[t]->sp->rname[k] = (cString) arr_alloc(NAMESIZE, char);
			strcpy(prob[t]->sp->rname[k], orig->rname[m]);

			prob[t]->bBar->val[prob[t]->bBar->cnt + 1] = orig->rhsx[m];
			prob[t]->bBar->col[prob[t]->bBar->cnt + 1] = m - tim->row[t] + 1;
			prob[t]->bBar->cnt++;
		}

		for (int r1 = 0; r1 < orig->objQ->cnt; r1++) {
			if ( (orig->objQ->col[r1] >= tim->col[t] & orig->objQ->col[r1] < tim->col[t + 1]) & (orig->objQ->row[r1] >= tim->col[t]) & (orig->objQ->row[r1] < tim->col[t + 1]) ) {
				/* there is a non-zero element in the quadratic matrix, so allocate memory */
				if ( prob[t]->sp->objQ == NULL ) {
					prob[t]->sp->objQ = (sparseMatrix*) mem_malloc(sizeof(sparseMatrix));
					prob[t]->sp->objQ->col = (iVector) arr_alloc(prob[t]->sp->mac * prob[t]->sp->mac, int);
					prob[t]->sp->objQ->row = (iVector) arr_alloc(prob[t]->sp->mac * prob[t]->sp->mac, int);
					prob[t]->sp->objQ->val = (dVector) arr_alloc(prob[t]->sp->mac * prob[t]->sp->mac, double);
					prob[t]->sp->objQ->cnt = 0;
				}

				prob[t]->sp->objQ->cnt++;
				prob[t]->sp->objQ->col[prob[t]->sp->objQ->cnt] = orig->objQ->col[r1] - tim->col[t] + 1;
				prob[t]->sp->objQ->row[prob[t]->sp->objQ->cnt] = orig->objQ->row[r1] - tim->col[t] + 1;
				prob[t]->sp->objQ->val[prob[t]->sp->objQ->cnt] = orig->objQ->val[r1];
 
			}
		}
	}

	/* Now copy the terminal stage problem */
	for ( m = tim->col[t]; m < orig->mac; m++ ) {
		k = m - tim->col[t];

		prob[t]->dBar->val[prob[t]->dBar->cnt+1] = orig->objx[m];
		prob[t]->dBar->col[prob[t]->dBar->cnt+1] = m - tim->col[t] + 1;
		prob[t]->dBar->cnt++;

		prob[t]->uBar->val[prob[t]->uBar->cnt + 1] = orig->bdu[m];
		prob[t]->uBar->col[prob[t]->uBar->cnt + 1] = m - tim->col[t] + 1;
		prob[t]->uBar->cnt++;

		prob[t]->lBar->val[prob[t]->lBar->cnt + 1] = orig->bdl[m];
		prob[t]->lBar->col[prob[t]->lBar->cnt + 1] = m - tim->col[t] + 1;
		prob[t]->lBar->cnt++;

		prob[t]->sp->objx[k] = orig->objx[m];
		prob[t]->sp->bdl[k] = orig->bdl[m];
		prob[t]->sp->bdu[k] = orig->bdu[m];

		if (orig->ctype[m] != 'C') {
			errMsg("setup", "newProb", "integer variable in non-root stage",0);
			return NULL;
		}
		else
			prob[t]->sp->ctype[k] = orig->ctype[m];

		prob[t]->sp->cname[k] = (cString) arr_alloc(NAMESIZE, char);
		strcpy(prob[t]->sp->cname[k], orig->cname[m]);
		prob[t]->sp->matcnt[k] = 0;
		if ( orig->matcnt[m] > 0 )
			for ( i = orig->matbeg[m]; i < orig->matbeg[m]+orig->matcnt[m]; i++ ) {
				/* The coefficient is part of the current stage constraint matrix */
				if ( k == 0 )
					prob[t]->sp->matbeg[k] = 0;
				else
					prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1] + prob[t]->sp->matcnt[k-1];
				prob[t]->sp->matval[prob[t]->sp->matsz] = orig->matval[i];
				prob[t]->sp->matind[prob[t]->sp->matsz] = orig->matind[i]-tim->row[t];
				++prob[t]->sp->matcnt[k];
				++prob[t]->sp->matsz;
				++prob[t]->sp->numnz;
				prob[t]->Dbar->val[prob[t]->Dbar->cnt+1] = orig->matval[i];
				prob[t]->Dbar->col[prob[t]->Dbar->cnt+1] = m-tim->col[t]+1;
				prob[t]->Dbar->row[prob[t]->Dbar->cnt+1] = orig->matind[i] - tim->row[t]+1;
				++prob[t]->Dbar->cnt;
			}
		else {
			if ( k == 0 )
				prob[t]->sp->matbeg[k] = 0;
			else
				prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1];
		}
	}

	/* copy row information for terminal stage */
	for ( m = tim->row[t]; m < orig->mar; m++ ) {
		k = m - tim->row[t];
		prob[t]->sp->rhsx[k] = orig->rhsx[m];
		prob[t]->sp->senx[k] = orig->senx[m];
		prob[t]->sp->rname[k] = (cString) arr_alloc(NAMESIZE, char);
		strcpy(prob[t]->sp->rname[k], orig->rname[m]);
		prob[t]->bBar->val[prob[t]->bBar->cnt+1] = orig->rhsx[m];
		prob[t]->bBar->col[prob[t]->bBar->cnt+1] = m - tim->row[t]+1;
		prob[t]->bBar->cnt++;
	}

	for (int r1 = 0; r1 < orig->objQ->cnt; r1++) {
		if ((orig->objQ->col[r1] >= tim->col[t]) & (orig->objQ->row[r1] >= tim->col[t]) ) {
			if ( prob[t]->sp->objQ == NULL ) {
				prob[t]->sp->objQ = (sparseMatrix*) mem_malloc(sizeof(sparseMatrix));
				prob[t]->sp->objQ->col = (iVector) arr_alloc(prob[t]->sp->mac * prob[t]->sp->mac, int);
				prob[t]->sp->objQ->row = (iVector) arr_alloc(prob[t]->sp->mac * prob[t]->sp->mac, int);
				prob[t]->sp->objQ->val = (dVector) arr_alloc(prob[t]->sp->mac * prob[t]->sp->mac, double);
				prob[t]->sp->objQ->cnt = 0;
			}

			prob[t]->sp->objQ->cnt++;
			prob[t]->sp->objQ->col[prob[t]->sp->objQ->cnt] = orig->objQ->col[r1]- tim->col[t] + 1;
			prob[t]->sp->objQ->row[prob[t]->sp->objQ->cnt] = orig->objQ->row[r1]- tim->col[t] + 1;
			prob[t]->sp->objQ->val[prob[t]->sp->objQ->cnt] = orig->objQ->val[r1];
		}
	}

	/* if integer or binary variables are encountered, then label the stage problem as a mixed integer LP */
	if ( prob[t]->sp->objQ != NULL ) {
		if (prob[t]->sp->numInt + prob[t]->sp->numBin > 0) {
			prob[t]->sp->type = PROB_MIQP;
		}
		else {
			prob[t]->sp->type = PROB_QP;
		}
	}
	else { 

		if (prob[t]->sp->numInt + prob[t]->sp->numBin > 0) {
			prob[t]->sp->type = PROB_MILP;
		}
		else {
			prob[t]->sp->type = PROB_LP;
		}
	}

#if defined(DECOMPOSE_CHECK)
	/* write stage problems in LP format to verify decomposition */
	char fname[BLOCKSIZE]; 
	for ( t = 0; t < tim->numStages; t++) {
		if ( !(prob[t]->sp->model = setupProblem(prob[t]->sp->name, prob[t]->sp->mac, prob[t]->sp->mar, prob[t]->sp->objSense, 0.0, prob[t]->sp->objx, prob[t]->sp->objQ, prob[t]->sp->senx,
				prob[t]->sp->rhsx, prob[t]->sp->matbeg, prob[t]->sp->matcnt, prob[t]->sp->matind, prob[t]->sp->matval, prob[t]->sp->bdl, prob[t]->sp->bdu, prob[t]->sp->ctype,
				prob[t]->sp->cname, prob[t]->sp->rname)) ) {
			errMsg("solver", "newProb", "failed to setup stage problem in solver", 0);
			goto TERMINATE;
		}
		sprintf(fname, "stageProb%d.lp", t);

		if ( writeProblem(prob[t]->sp->model, fname) ) {
			errMsg("solver", "newProb", "failed to write stage problem", 0);
			goto TERMINATE;
		}
	}
#endif

	/* save size information in numType */
	for ( t = 0; t < tim->numStages; t++ ) {
		if ( !(prob[t]->num = (numType *) mem_malloc(sizeof(numType))) )
			errMsg("allocation", "newProb", "prob[t]->num",0);

		prob[t]->num->cols = prob[t]->sp->mac;
		prob[t]->num->rows = prob[t]->sp->mar;
		prob[t]->num->intCols = prob[t]->sp->numInt;
		prob[t]->num->binCols = prob[t]->sp->numBin;
		prob[t]->num->numRV = prob[t]->num->rvColCnt = prob[t]->num->rvRowCnt = 0;
		prob[t]->num->rvAOmCnt = prob[t]->num->rvBOmCnt = prob[t]->num->rvCOmCnt = prob[t]->num->rvDOmCnt = 0;
		prob[t]->num->rvaOmCnt = prob[t]->num->rvbOmCnt = prob[t]->num->rvcOmCnt = prob[t]->num->rvdOmCnt = 0;
		prob[t]->num->rvylOmCnt = prob[t]->num->rvyuOmCnt = 0;

		if ( t == 0 ) {
			prob[t]->mean = NULL;
			prob[t]->coord = NULL;
			prob[t]->num->prevCols = prob[t]->num->prevRows = prob[t]->num->cntCcols = prob[t]->num->cntCrows = 0;
		}
		else {
			if ( !(prob[t]->coord = (coordType *) mem_malloc(sizeof(coordType))) )
				errMsg("allocation", "newProb", "prob[t]->coord",0);

			prob[t]->num->prevCols = prob[t-1]->num->cols;
			prob[t]->num->prevRows = prob[t-1]->num->rows;

			prob[t]->coord->CCols = findElems(prob[t]->Cbar->col, prob[t]->Cbar->cnt, &prob[t]->num->cntCcols);
			prob[t]->coord->CRows = findElems(prob[t]->Cbar->row, prob[t]->Cbar->cnt, &prob[t]->num->cntCrows);

			prob[t]->coord->allRVCols = prob[t]->coord->allRVRows = prob[t]->coord->rvCols = prob[t]->coord->rvRows = NULL;
			prob[t]->coord->rvCOmCols = prob[t]->coord->rvCOmRows = prob[t]->coord->rvbOmRows = prob[t]->coord->rvdOmCols = NULL;
			prob[t]->coord->rvylOmRows = prob[t]->coord->rvyuOmRows = NULL;
		}
	}

	/* decompose the stochastic elements of the problem. Go through the list of random variable and assign them to
	 * appropriate parts (right-hand side and objective coefficients). */
	int rvOffset = 0;
	int rhs = 0;
	int upper = 0;
	int cost = 0;
	int transfer = 0;

	for ( m = 0; m < (*stoc)->numOmega; m++ ) {
		if ( (*stoc)->col[m] == -1 ) {
			/* randomness in right-hand side */
			t = 0;
			while ( t < tim->numStages ) {
				prob[t]->num->rvyuOmCnt = 0;
				if ( (*stoc)->row[m] < tim->row[t] )
					break;
				t++;
			}
			t--;
		}

		else {
			/* randomness in either objective function coefficients or the transfer matrix */
			t = 0;
			while ( t < tim->numStages ) {
				if ( (*stoc)->col[m] < tim->col[t] )
					break;
				t++;
			}
			t--;
		}

		if ( t == 0 ) {
			errMsg("setup", "newProb", "encountered randomness in root-stage", 0);
			goto TERMINATE;
		}

		if ( prob[t]->num->numRV == 0) {
			if ( !(prob[t]->coord->allRVCols = (iVector) arr_alloc((*stoc)->numOmega+1, int)) )
				errMsg("allocation", "newProb", "prob->coord->allRVCols", 0);
			if ( !(prob[t]->coord->allRVRows= (iVector) arr_alloc((*stoc)->numOmega+1, int)) )
				errMsg("allocation", "newProb", "prob->coord->allRVRows", 0);
			if ( !(prob[t]->coord->rvOffset = (iVector) arr_alloc(5, int)))
				errMsg("allocation", "newProb", "prob->coord->rOffset", 0);
			if ( !(prob[t]->mean = (dVector) arr_alloc((*stoc)->numOmega+1, double)) )
				errMsg("allocation", "newProb", "prob->mean", 0);
			prob[t]->omegaBeg = m;
			rvOffset = 0;
		}

		prob[t]->num->numRV++;
		prob[t]->mean[prob[t]->num->numRV] = (*stoc)->mean[m];

		/* The order of random variables is: (i) right-hand side, (ii) transfer matrix, and (iii) cost-coefficients (iiii) upperbounds. */

		if ( (*stoc)->col[m] == -1 && (*stoc)->row[m] != -1 ) {

			/* Right-hand side */
			if ( prob[t]->num->rvbOmCnt == 0 ) {
				prob[t]->coord->rvbOmRows = (iVector) arr_alloc(prob[t]->num->rows , int);
				prob[t]->coord->rvOffset[0] = rvOffset;

			}

			prob[t]->coord->allRVCols[prob[t]->num->numRV] = -1;
			prob[t]->coord->allRVRows[prob[t]->num->numRV] = (*stoc)->row[m] - tim->row[t] + 1;
			prob[t]->coord->rvbOmRows[++prob[t]->num->rvbOmCnt] = prob[t]->coord->allRVRows[prob[t]->num->numRV];
			rhs++;
		}
		/*upper bound*/
		else if ((*stoc)->col[m] != -1 && (*stoc)->row[m] == -2) {
			upper++;
			if (prob[t]->num->rvyuOmCnt == 0) {
				prob[t]->coord->rvyuOmRows = (iVector) arr_alloc(prob[t]->num->cols + 1, int);
			}
			prob[t]->coord->rvyuOmRows[prob[t]->num->rvyuOmCnt+1] = (*stoc)->col[m] - tim->col[t] + 1;
			prob[t]->num->rvyuOmCnt++;
		}
		/*lower bound*/
		else if ((*stoc)->col[m] != -1 && (*stoc)->row[m] == -3) {
			upper++;
			if (prob[t]->num->rvylOmCnt == 0) {
				prob[t]->coord->rvylOmRows = (iVector) arr_alloc(prob[t]->num->cols + 1, int);
			}
			prob[t]->coord->rvylOmRows[prob[t]->num->rvylOmCnt +1] = (*stoc)->col[m] - tim->col[t] + 1;
			prob[t]->num->rvylOmCnt++;
		}

		else if ( (*stoc)->col[m] != -1 && (*stoc)->row[m] != -1 ) {
			/* Transfer matrix */
			if ( prob[t]->num->rvCOmCnt == 0 ) {
				prob[t]->coord->rvCOmCols = (iVector) arr_alloc((*stoc)->numOmega, int);
				prob[t]->coord->rvCOmRows = (iVector) arr_alloc((*stoc)->numOmega, int);
				prob[t]->coord->rvOffset[1] = rvOffset;
				transfer++;
			}
			prob[t]->coord->allRVCols[prob[t]->num->numRV] = (*stoc)->col[m] - tim->col[t] + 1;
			prob[t]->coord->allRVRows[prob[t]->num->numRV] = (*stoc)->row[m] - tim->row[t] + 1;
			prob[t]->num->rvCOmCnt++;
			prob[t]->coord->rvCOmCols[prob[t]->num->rvCOmCnt] = prob[t]->coord->allRVCols[prob[t]->num->numRV];
			prob[t]->coord->rvCOmRows[prob[t]->num->rvCOmCnt] = prob[t]->coord->allRVRows[prob[t]->num->numRV];

		}
		else {
			/* Cost coefficients */
			if ( prob[t]->num->rvdOmCnt == 0 ) {
				prob[t]->coord->rvdOmCols = (iVector) arr_alloc((*stoc)->numOmega+1,int);
				prob[t]->coord->rvOffset[2] = rvOffset;
			}
			prob[t]->coord->allRVCols[prob[t]->num->numRV] = (*stoc)->col[m] - tim->col[t] + 1;
			prob[t]->coord->allRVRows[prob[t]->num->numRV] = -1;
			prob[t]->coord->rvdOmCols[++prob[t]->num->rvdOmCnt] = prob[t]->coord->allRVCols[prob[t]->num->numRV];
			cost++;
		}
		/*update rvofset*/
		prob[t]->coord->rvOffset[0] = 0;
		prob[t]->coord->rvOffset[1] = rhs ;
		prob[t]->coord->rvOffset[2] = rhs + transfer ;
		prob[t]->coord->rvOffset[3] = rhs + transfer + cost ;
		prob[t]->coord->rvOffset[4] = rhs + transfer + cost + upper;
	}

	for ( t = 1; t < tim->numStages; t++ ) {
		prob[t]->coord->rvCols = findElems(prob[t]->coord->allRVCols, prob[t]->num->numRV, &prob[t]->num->rvColCnt);
		prob[t]->coord->rvRows = findElems(prob[t]->coord->allRVRows, prob[t]->num->numRV, &prob[t]->num->rvRowCnt);
	}

	/* Modify the dBar, bBar and Cbar with mean values computed from stoch file */
	rvOffset = 0;

	for ( t = 1; t < tim->numStages; t++ ) {

		/* Right-hand side */
		for ( m = 1; m <= prob[t]->num->rvbOmCnt; m++ ) {
			i = 1;
			while ( i <= prob[t]->bBar->cnt ) {
				if ( prob[t]->bBar->col[i] == prob[t]->coord->rvbOmRows[m])
					break;
				i++;
			}
			prob[t]->bBar->val[i] = (*stoc)->mean[rvOffset + prob[t]->coord->rvOffset[0]+m-1];
		}



		/* Upper bound  */
		for (m = 1; m <= prob[t]->num->rvyuOmCnt; m++) {
			i = 1;
			while (i <= prob[t]->uBar->cnt) {
				if (prob[t]->uBar->col[i] == prob[t]->coord->rvyuOmRows[m])
					break;
				i++;
			}
			prob[t]->uBar->val[i] = (*stoc)->mean[rvOffset + prob[t]->coord->rvOffset[3] + m -1];
		}


		/* Lower bound */
		for (m = 1; m <= prob[t]->num->rvylOmCnt; m++) {
			i = 1;
			while (i <= prob[t]->lBar->cnt) {
				if (prob[t]->lBar->col[i] == prob[t]->coord->rvyuOmRows[m])
					break;
				i++;
			}
			prob[t]->lBar->val[i] = (*stoc)->mean[rvOffset + prob[t]->coord->rvOffset[4] + m - 1];
		}

		/* Transfer matrix */
		for ( m = 1; m <= prob[t]->num->rvCOmCnt; m++ ) {
			i = 1;
			while ( i <= prob[t]->num->rvCOmCnt ) {
				if ( prob[t]->Cbar->col[i] == prob[t]->coord->rvCOmCols[m] && prob[t]->Cbar->row[i] == prob[t]->coord->rvCOmRows[m] )
					break;
				i++;
			}
			prob[t]->Cbar->val[i] = (*stoc)->mean[rvOffset + prob[t]->coord->rvOffset[1]+m-1];
		}

		/* Cost coefficients */
		for ( m = 1; m <= prob[t]->num->rvdOmCnt; m++ ) {
			i = 1;
			while ( i <= prob[t]->dBar->cnt ) {
				if ( prob[t]->dBar->col[i] == prob[t]->coord->rvdOmCols[m])
					break;
				i++;
			}
			prob[t]->dBar->val[i] = (*stoc)->mean[rvOffset + prob[t]->coord->rvOffset[2]+m-1];
		}
		rvOffset += prob[t]->num->numRV;
	}

	printf("2. Decomposition complete.\n");

	/* Solve the mean value problem */
	if ( (meanX = meanProblem(orig, (*stoc))) == NULL) {
		errMsg("setup", "newProbwSMPS", "failed to solve the mean-value problem", 0);
		goto TERMINATE;	}

	/* Compute the stage-wise lower bounds -
	TO DO: NEEDS TO BE UPDATED TO ACCOUNT RANDOMNESS IN BOUNDS */
	lb = calcLowerBound(orig, tim, (*stoc));

	int offset = 0;
	for ( t = 0; t < (*numStages) - 1; t++ ) {
		prob[t]->meanX = duplicVector(meanX+offset, prob[t]->num->cols);
		prob[t]->lb = lb[t];
		offset += prob[t]->num->cols;
	}
	prob[t]->meanX = NULL;

	if (meanX) mem_free(meanX);

	if (tim) freeTimeType(tim);
	if (orig) freeOneProblem(orig);
	mem_free(lb);
	return prob;

	TERMINATE:
	if (meanX) mem_free(meanX);
	if (prob) freeProbType(prob, (*numStages));
	if (tim) freeTimeType(tim);
	if (orig) freeOneProblem(orig);
	mem_free(lb);

	return NULL;
}//END newProbwSMPS()

/* setup and solve the original problem _orig_ with expected values for all random variables provided in _stoc_. If the problem is an mixed-integer program,
 *  then a relaxed problem is solved. The function returns a dVector of mean value solutions, if there is an error it returns NULL.*/
dVector meanProblem(oneProblem *orig, stocType *stoc) {
	dVector	xk;
	double	obj = 0.0;
	int 	n, status;

	printf("3. Solving the mean-value problem.\n");

	/* The original problem resides in the solver under orig->model. This model instance was created in readCore(). */
	/* change the coefficients and right-hand side to the mean values */
	for (n = 0; n < stoc->numOmega; n++ ) {
		if ( stoc->row[n] == -1) {
			status = changeObjCoeff (orig->model, stoc->col[n], stoc->mean[n]);
			if ( status ) {
				errMsg("setup", "meanProblem", "failed to change the coefficients with mean values", 0);
				return NULL;
			}
		}else if (stoc->row[n] == -2) {
			status = changeBDSelement(orig->model,"UB", stoc->col[n], stoc->mean[n]);
			if (status) {
				errMsg("setup", "meanProblem", "failed to change the coefficients with mean values for upperbounds", 0);
				return NULL;
			}
		}
		else if (stoc->row[n] == -3) {
			status = changeBDSelement(orig->model, "LB", stoc->col[n], stoc->mean[n]);
			if (status) {
				errMsg("setup", "meanProblem", "failed to change the coefficients with mean values for upperbounds", 0);
				return NULL;
			}
		}

		else if (stoc->col[n] == -1 ) {
			status = changeRHSelement (orig->model, stoc->row[n], stoc->mean[n]);
			if ( status ) {
				errMsg("setup", "meanProblem", "failed to change the coefficients with mean values", 0);
				return NULL;
			}
		}
		else {
			int vind[]   = {stoc->col[n]};
			int cind[]   = {stoc->row[n]};
			double val[] = {stoc->mean[n]};
			status = changeCoefficients(orig->model, 1, cind, vind, val);
			if ( status ) {
				errMsg("setup", "meanProblem", "failed to change the constraint coefficients with mean values", 0);
				return NULL;
			}
		}
	}

	/* write the mean value problem */
	status = writeProblem(orig->model, "meanValueProblem.lp");
	if ( status ) {
		errMsg("solver", "meanProblem", "failed to write the problem", 0);
		return NULL;
	}

	/* solve the mean value problem */
	status = solveProblem(orig->model);
	if ( status ) {
		errMsg("setup", "meanProblem", "failed to solve mean value problem", 0);
		return NULL;
	}

	/* print objective functiom value */
	obj = getObjective(orig->model);

	/* obtain the primal solution */
	xk = (dVector) arr_alloc(orig->mac+1, double);
	if ( getPrimal(orig->model, xk, 0, orig->mac) ) {
		errMsg("solver", "meanProblem", "failed to retrieve primal solution for the mean-value problem", 0);
		return NULL;
	}

	printf("\tOptimal objective function value for (relaxed) mean value problem = %lf\n", obj);

	return xk;
}//END meanProblem()

/* Calculate the stage-wise lower bound. */
dVector calcLowerBound(oneProblem *orig, timeType *tim, stocType *stoc) {
	sparseMatrix *Cbar;
	dVector		duals, vals, beta, dj, u, bBar;
	iVector		colIdx, rowIdx;
	double		alpha;
	modelPtr	*lpClone;

	printf("4. Computing stage-wise lower bounds.\n");

	duals = (dVector) arr_alloc(orig->mar+1, double);

	dj = (dVector) arr_alloc(orig->mac+1, double);
	u = (dVector) arr_alloc(orig->mac+1, double);

	vals = (dVector) arr_alloc(orig->mac, double);
	Cbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix));
	Cbar->col = (iVector) arr_alloc(orig->matsz+1, int);
	Cbar->row = (iVector) arr_alloc(orig->matsz+1, int);
	Cbar->val = (dVector) arr_alloc(orig->matsz+1, double);
	beta = (dVector) arr_alloc(orig->mac+1, double);

	dVector lb = (dVector) arr_alloc(tim->numStages-1, double);

	colIdx = (iVector) arr_alloc(orig->mac, int);
	for ( int n = 0; n < orig->mac; n++ ) {
		colIdx[n] = n;
	}

	rowIdx = (iVector) arr_alloc(orig->mar, int);
	for ( int n = 0; n < orig->mar; n++ ) {
		rowIdx[n] = n;
	}

	/* obtain primal and dual solutions, column status, and dual slacks from the mean value solve */
	if ( getPrimal(orig->model, u, 0, orig->mac) ) {
		errMsg("solver", "calcLowerBound", "failed to obtain primal solution", 1);
	}
	if ( getDual(orig->model, duals, 0, orig->mar) ) {
		errMsg("setup", "calcLowerBound", "failed to obtain dual for the mean value problem", 1);
	}
	if (getDualSlack (orig->model, dj, 0, orig->mac) ) {
		errMsg("solver", "calcLowerBound", "failed to obtain dual slacks", 1);
	}

	/* Extract bBar */
	bBar = duplicVector(orig->rhsx-1, orig->mar);
	for ( int n = 0; n < stoc->numOmega; n++ )
		if (stoc->col[n] == -1) {
			bBar[stoc->row[n]] = stoc->mean[n];
		}

	printf("\tLower bounds computed = ");

	for ( int t = 1; t < tim->numStages; t++ ) {
		/* 1. clone to problem to be used for computing the lower bound */
		lpClone = cloneProblem(orig->model);

		/* 2b. extract Cbar */
		Cbar->cnt = 0;
		for (int n = tim->col[t-1]; n < tim->col[t]; n++) {
			for (int m = orig->matbeg[n]; m < orig->matbeg[n] + orig->matcnt[n]; m++) {
				if (orig->matind[m] >= tim->row[t]) {
					/* The coefficient is part of the subproblem's T matrix */
					Cbar->val[Cbar->cnt + 1] = orig->matval[m];
					Cbar->row[Cbar->cnt + 1] = orig->matind[m] - tim->row[t] + 1;
					Cbar->col[Cbar->cnt + 1] = n + 1;
					++Cbar->cnt;
				}
			}
		}

		/* 3a. Compute alpha + mubBar */
		alpha = vXv(bBar, duals+tim->row[t], NULL, orig->mar-tim->row[t]) + vXv(dj+tim->col[t], u, NULL, orig->mac-tim->col[t]);

		/* 3b. Compute beta */
		for ( int n = 0; n <= tim->col[t]; n++ )
			beta[n] = 0.0;
		for ( int n = 1; n <= Cbar->cnt; n++ )
			beta[Cbar->col[n]] += duals[tim->row[t] + Cbar->row[n]] * Cbar->val[n];

		/* 4a. Use mean-cut coefficients to create the lower-bounding problem for previous stage*/
		int startCol = (t == 0) ? 0:tim->col[t-1];
		int numCols  = tim->col[t] - startCol;
		for ( int n = startCol; n < tim->col[t]; n++) {
			vals[n-startCol] = orig->objx[n] - beta[n+1];
		}

		/* 4b. Change the objective in solver to prepare for lower bound calculation */
		if ( changeObjCoeffArray(lpClone, numCols, colIdx+startCol, vals) ) {
			errMsg("setup", "calcLowerBound", "failed to change objective coefficients in solver while computing lower bound", 1);
		}

		/* 5. Remove the columns and rows that corresponding to stages starting from t */
		if ( deleteCols(lpClone, orig->mac - tim->col[t], colIdx+tim->col[t]) ) {
			errMsg("setup", "calcLowerBound", "failed to remove columns corresponding to future stages", 1);
		}

		if ( deleteRows(lpClone, orig->mar - tim->row[t], rowIdx+tim->row[t]) ) {
			errMsg("setup", "calcLowerBound", "failed to remove columns corresponding to future stages", 1);
		}

#ifdef DECOMPOSE_CHECK
		writeProblem(lpClone, "lowerBoundCalc.lp");
#endif

		/* 6a. Solve the relaxed problem */
		if (solveProblem(lpClone) ) {
			errMsg("setup", "calcLowerBound", "failed to solve problem computing lower bound", 1);
		}

		/* 6b. Get lower bound */
		double obj = getObjective(lpClone) + alpha;
		lb[t-1] = minimum(obj, 0.0);

		/* release the problem */
		if ( freeProblem(lpClone) ) {
			errMsg("setup", "calcLowerBound", "failed to free problem", 1);
		}
	}

	printVector(lb-1, tim->numStages-1, NULL);

	mem_free(colIdx); mem_free(rowIdx);
	mem_free(vals);
	mem_free(beta);
	mem_free(duals);
	mem_free(bBar);
	mem_free(dj);
	mem_free(u);
	freeSparseMatrix(Cbar);

	return lb;
}//END calcLowerBound()

/* free up the probType. Needs number of stages as input */
void freeProbType(probType **prob, int T) {
	int t;

	if ( prob ) {
		for ( t = 0; t < T; t++ ) {
			if (prob[t]) {
				if (prob[t]->Abar) freeSparseMatrix(prob[t]->Abar);
				if (prob[t]->Bbar) freeSparseMatrix(prob[t]->Bbar);
				if (prob[t]->Cbar) freeSparseMatrix(prob[t]->Cbar);
				if (prob[t]->Dbar) freeSparseMatrix(prob[t]->Dbar);
				if (prob[t]->aBar) freeSparseVector(prob[t]->aBar);
				if (prob[t]->bBar) freeSparseVector(prob[t]->bBar);
				if (prob[t]->cBar) freeSparseVector(prob[t]->cBar);
				if (prob[t]->dBar) freeSparseVector(prob[t]->dBar);
				if (prob[t]->uBar) freeSparseVector(prob[t]->uBar);
				if (prob[t]->lBar) freeSparseVector(prob[t]->lBar);
				if (prob[t]->num) mem_free(prob[t]->num);
				if (prob[t]->coord) freeCoordType(prob[t]->coord);
				if (prob[t]->mean) mem_free(prob[t]->mean);
				if (prob[t]->meanX) mem_free(prob[t]->meanX);

				/* Elements cannot be freed in freeOneProblem() */
				if (prob[t]->sp) {
					mem_free(prob[t]->sp->name);
					for ( int n = 0; n < prob[t]->sp->mac; n++ ) {
						mem_free(prob[t]->sp->cname[n]);
					}
					for ( int n = 0; n < prob[t]->sp->mar; n++ ) {
						mem_free(prob[t]->sp->rname[n]);
					}
					freeOneProblem(prob[t]->sp);
				}
				mem_free(prob[t]);
			}
		}
		mem_free(prob);
	}

}//END freeProb()

/* free up the coordType */
void freeCoordType (coordType *coord) {

	if (coord->allRVCols) mem_free(coord->allRVCols);
	if (coord->allRVRows) mem_free(coord->allRVRows);
	if (coord->CCols) mem_free(coord->CCols);
	if (coord->CRows) mem_free(coord->CRows);
	if (coord->rvCols) mem_free(coord->rvCols);
	if (coord->rvRows) mem_free(coord->rvRows);
	if (coord->rvbOmRows) mem_free(coord->rvbOmRows);
	if (coord->rvdOmCols) mem_free(coord->rvdOmCols);
	if (coord->rvCOmCols) mem_free(coord->rvCOmCols);
	if (coord->rvCOmRows) mem_free(coord->rvCOmRows);
	if (coord->rvylOmRows) mem_free(coord->rvylOmRows);
	if (coord->rvyuOmRows) mem_free(coord->rvyuOmRows);
	if (coord->rvOffset) mem_free(coord->rvOffset);
	mem_free(coord);

}//END freeCoordType()

void printDecomposeSummary(FILE *fptr, cString probName, timeType *tim, probType **prob) {
	int t;

	fprintf(fptr, "====================================================================================================================================\n");
	fprintf(fptr, "Stage optimization problem for given input s_t = (x_t, omega_t) is in the following form : \n\n");
	fprintf(fptr, "\t\t\t\t h_t(s_t) = c_t*x_t + min d_t*u_t\n");
	fprintf(fptr, "\t\t\t\t                      s.t. D_t u_t = b_t - C_t x_t,\n\n");
	fprintf(fptr, "with the following linear dynamics: x_{t+} = a_{t+} + A_{t+}x_t + B_{t+}u_t.\n\n");
	fprintf(fptr, "------------------------------------------------------------------------------------------------------------------------------------\n");
	fprintf(fptr, "\n====================================================================================================================================\n");
	fprintf(fptr, "-------------------------------------------------------- Problem Information -------------------------------------------------------\n");
	fprintf(fptr, "====================================================================================================================================\n");
	fprintf(fptr, "Problem                            : %s\n", probName);
	fprintf(fptr, "Number of stages                   : %d\n", tim->numStages);
	for ( t = 0; t < tim->numStages; t++ ) {
		fprintf(fptr,  "------------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(fptr,  "Stage %d\n", t);
		fprintf(fptr,  "Number of decision variables (u_t) = %d\t\t", prob[t]->sp->mac);
		fprintf(fptr,  "(Continuous = %d\tInteger = %d\tBinary = %d)\n", prob[t]->sp->mac - prob[t]->sp->numInt - prob[t]->sp->numBin, prob[t]->sp->numInt, prob[t]->sp->numBin);
		fprintf(fptr,  "Number of constraints              = %d\n", prob[t]->sp->mar);
		if ( prob[t]->num->numRV != 0 ) {
			fprintf(fptr,  "Number of random variables (omega) = %d\t\t", prob[t]->num->numRV);
			fprintf(fptr,  "(a_t = %d; b_t = %d; c_t = %d; d_t = %d; A_t = %d; B_t = %d; C_t = %d; D_t = %d)\n", prob[t]->num->rvaOmCnt, prob[t]->num->rvbOmCnt, prob[t]->num->rvcOmCnt, prob[t]->num->rvdOmCnt,
					prob[t]->num->rvAOmCnt, prob[t]->num->rvBOmCnt, prob[t]->num->rvCOmCnt, prob[t]->num->rvDOmCnt);
		}
		else
			fprintf(fptr,  "Number of random variables (omega) = 0\n");
	}
	fprintf(fptr, "====================================================================================================================================\n");

}//printDecomposeSummary()
