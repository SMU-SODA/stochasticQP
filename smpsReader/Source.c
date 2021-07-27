/*
 * algorithm.c
 *  Created on: Jul , 2021
 *      Author: Niloofar Fadavi, Harsha Gangammanavar
 */
#include "../solverUtilities/utilities.h"
#include "../solverUtilities/solver_gurobi.h"
#include "smps.h"
#include "prob.h"
#include "stochasticQP.h"

/*building master problem*/
cellType* buildCell(probType** prob) {
	oneProblem* stage0 = NULL;
    stage0 = (oneProblem*)mem_malloc(sizeof(oneProblem));
	stage0->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
	stage0->objQ->col = (iVector)arr_alloc(prob[0]->sp->mac * prob[0]->sp->mac, int);
	stage0->objQ->row = (iVector)arr_alloc(prob[0]->sp->mac * prob[0]->sp->mac, int);
	stage0->objQ->val = (dVector)arr_alloc(prob[0]->sp->mac * prob[0]->sp->mac, double);
	stage0->model = NULL;
	stage0->name = (cString)arr_alloc(NAMESIZE, char);
	stage0->objname = (cString)arr_alloc(NAMESIZE, char);
	stage0->objx = (dVector)arr_alloc(prob[0]->sp->macsz, double);
	stage0->bdl = (dVector)arr_alloc(prob[0]->sp->macsz, double);
	stage0->bdu = (dVector)arr_alloc(prob[0]->sp->macsz, double);
	stage0->ctype = (cString)arr_alloc(prob[0]->sp->macsz, char);
	stage0->rhsx = (dVector)arr_alloc(prob[0]->sp->marsz, double);
	stage0->senx = (cString)arr_alloc(prob[0]->sp->marsz, char);
	stage0->matbeg = (iVector)arr_alloc(prob[0]->sp->macsz, int);
	stage0->matcnt = (iVector)arr_alloc(prob[0]->sp->macsz, int);
	stage0->matval = (dVector)arr_alloc(prob[0]->sp->matsz, double);
	stage0->matind = (iVector)arr_alloc(prob[0]->sp->matsz, int);
	stage0->cname = (cString*)arr_alloc(prob[0]->sp->macsz, cString);
	stage0->rname = (cString*)arr_alloc(prob[0]->sp->marsz, cString);
	stage0->mac = prob[0]->sp->mac;
	stage0->type= prob[0]->sp->type;			/* type of problem: LP, QP, MIP or MIQP */

	for (int i = 0; i < prob[0]->sp->mac; i++) {
		stage0->objx[i] = prob[0]->sp->objx[i];
	}

	stage0->objQ->cnt = prob[0]->sp->objQ->cnt;

	for (int i = 0; i < prob[0]->sp->objQ->cnt; i++) {
		stage0->objQ->val[i] = prob[0]->sp->objQ->val[i];
		stage0->objQ->col[i] = prob[0]->sp->objQ->col[i];
	}

	stage0->objSense = prob[0]->sp->objSense;		
	stage0->mac = prob[0]->sp->mac;			/* number of columns */
	stage0->mar = prob[0]->sp->mar;			/* number of rows */
	stage0->numBin = prob[0]->sp->numBin;			/* number of binary variables in the problem */
	stage0->numInt = prob[0]->sp->numInt;			/* number of integer variables  (includes both binary and general integer variables) in the problem */
	stage0->numnz = prob[0]->sp->numnz;			/* number of non-zero elements in constraint matrix */
	stage0->macsz = prob[0]->sp->macsz;			/* number of columns */
	stage0->marsz = prob[0]->sp->marsz;			/* number of rows */
	stage0->matsz = prob[0]->sp->matsz;			/* number of rows */
	

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->bdl[i] = prob[0]->sp->bdl[i];
	}

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->bdu[i] = prob[0]->sp->bdu[i];
	}

	for (int i = 0; i < prob[0]->sp->marsz; i++) {
		stage0->rhsx[i] = prob[0]->sp->rhsx[i];
	}

	for (int i = 0; i < prob[0]->sp->marsz; i++) {
		stage0->senx[i] = prob[0]->sp->senx[i];
	}

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->matbeg[i] = prob[0]->sp->matbeg[i];
	}

	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->matcnt[i] = prob[0]->sp->matcnt[i];
	}

	for (int i = 0; i < prob[0]->sp->matsz; i++) {
		stage0->matind[i] = prob[0]->sp->matind[i];
	}


	for (int i = 0; i < prob[0]->sp->matsz; i++) {
		stage0->matval[i] = prob[0]->sp->matval[i];
	}


	for (int i = 0; i < NAMESIZE; i++) {
		stage0->objname[i] = prob[0]->sp->objname[i];
	}

	for (int i = 0; i < NAMESIZE; i++) {
		stage0->name[i] = prob[0]->sp->name[i];
	}

	for (int i = 0; i < prob[0]->sp->marsz; i++) {
		stage0->rname[i] = prob[0]->sp->rname[i];
	}
	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->cname[i] = prob[0]->sp->cname[i];
	}



	for (int i = 0; i < prob[0]->sp->macsz; i++) {
		stage0->ctype[i] = prob[0]->sp->ctype[i];
	}
	cString	ctype;			/* type of decision variables: 'C' continuous, 'B' binary, 'I' general integer, 'S' semi-continuous, 'N' semi-integer */
	
	/* These don't seem necessary for Gurobi. */
	//	int		rstorsz;		/* memory size for storing row names */
	//	cString	rstore;			/* row names cString */
	//	int		cstorsz;		/* memory size for storing column names */
	//	cString	cstore;			/* column name cString */

	cellType* prb;
	prb = (cellType*)mem_malloc(sizeof(cellType));

	prb->master = stage0;
	oneProblem* stage1 = NULL;
	stage1 = (oneProblem*)mem_malloc(sizeof(oneProblem));
	stage1->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));   /*why do we write it?*/
	stage1->objQ->col = (iVector)arr_alloc(prob[1]->sp->mac * prob[1]->sp->mac, int);
	stage1->objQ->row = (iVector)arr_alloc(prob[1]->sp->mac * prob[1]->sp->mac, int);
	stage1->objQ->val = (dVector)arr_alloc(prob[1]->sp->mac * prob[1]->sp->mac, double);
	stage1->model = NULL;
	stage1->name = (cString)arr_alloc(NAMESIZE, char);
	stage1->objname = (cString)arr_alloc(NAMESIZE, char);
	stage1->objx = (dVector)arr_alloc(prob[1]->sp->macsz, double);
	stage1->bdl = (dVector)arr_alloc(prob[1]->sp->macsz, double);
	stage1->bdu = (dVector)arr_alloc(prob[1]->sp->macsz, double);
	stage1->ctype = (cString)arr_alloc(prob[1]->sp->macsz, char);
	stage1->rhsx = (dVector)arr_alloc(prob[1]->sp->marsz, double);
	stage1->senx = (cString)arr_alloc(prob[1]->sp->marsz, char);
	stage1->matbeg = (iVector)arr_alloc(prob[1]->sp->macsz, int);
	stage1->matcnt = (iVector)arr_alloc(prob[1]->sp->macsz, int);
	stage1->matval = (dVector)arr_alloc(prob[1]->sp->matsz, double);
	stage1->matind = (iVector)arr_alloc(prob[1]->sp->matsz, int);
	stage1->cname = (cString*)arr_alloc(prob[1]->sp->macsz, cString);
	stage1->rname = (cString*)arr_alloc(prob[1]->sp->marsz, cString);
	stage1->mac = prob[1]->sp->mac;
	stage1->type = prob[1]->sp->type;			/* type of problem: LP, QP, MIP or MIQP */
	for (int i = 0; i < prob[1]->sp->mac; i++) {
		stage1->objx[i] = prob[1]->sp->objx[i];
	}

	stage1->objQ->cnt = prob[1]->sp->objQ->cnt;

	for (int i = 0; i < prob[1]->sp->objQ->cnt; i++) {
		stage1->objQ->val[i] = prob[1]->sp->objQ->val[i];
		stage1->objQ->col[i] = prob[1]->sp->objQ->col[i];
	}

	stage1->objSense = prob[1]->sp->objSense;
	stage1->mac = prob[1]->sp->mac;			/* number of columns */
	stage1->mar = prob[1]->sp->mar;			/* number of rows */
	stage1->numBin = prob[1]->sp->numBin;			/* number of binary variables in the problem */
	stage1->numInt = prob[1]->sp->numInt;			/* number of integer variables  (includes both binary and general integer variables) in the problem */
	stage1->numnz = prob[1]->sp->numnz;			/* number of non-zero elements in constraint matrix */
	stage1->macsz = prob[1]->sp->macsz;			/* number of columns */
	stage1->marsz = prob[1]->sp->marsz;			/* number of rows */
	stage1->matsz = prob[1]->sp->matsz;			/* number of rows */


	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->bdl[i] = prob[1]->sp->bdl[i];
	}

	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->bdu[i] = prob[1]->sp->bdu[i];
	}

	for (int i = 0; i < prob[1]->sp->marsz; i++) {
		stage1->rhsx[i] = prob[1]->sp->rhsx[i];
	}

	for (int i = 0; i < prob[1]->sp->marsz; i++) {
		stage1->senx[i] = prob[1]->sp->senx[i];
	}

	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->matbeg[i] = prob[1]->sp->matbeg[i];
	}

	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->matcnt[i] = prob[1]->sp->matcnt[i];
	}

	for (int i = 0; i < prob[1]->sp->matsz; i++) {
		stage1->matind[i] = prob[1]->sp->matind[i];
	}


	for (int i = 0; i < prob[1]->sp->matsz; i++) {
		stage1->matval[i] = prob[1]->sp->matval[i];
	}


	for (int i = 0; i < NAMESIZE; i++) {
		stage1->objname[i] = prob[1]->sp->objname[i];
	}

	for (int i = 0; i < NAMESIZE; i++) {
		stage1->name[i] = prob[1]->sp->name[i];
	}

	for (int i = 0; i < prob[1]->sp->marsz; i++) {
		stage1->rname[i] = prob[1]->sp->rname[i];
	}
	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->cname[i] = prob[1]->sp->cname[i];
	}



	for (int i = 0; i < prob[1]->sp->macsz; i++) {
		stage1->ctype[i] = prob[1]->sp->ctype[i];
	}

	/* These don't seem necessary for Gurobi. */
	//	int		rstorsz;		/* memory size for storing row names */
	//	cString	rstore;			/* row names cString */
	//	int		cstorsz;		/* memory size for storing column names */
	//	cString	cstore;			/* column name cString */
	prb->subprob = stage1;

	return prb;
	}


/*This function gets the stoc file and creats a matrix of observations*/


sparseMatrix* buildObs(stocType* stoc) {
	sparseMatrix* obs;


}


