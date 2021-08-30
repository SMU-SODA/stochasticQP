/*
 * solver_gurobi.h
 *
 *  Created on: Dec 29, 2020
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef SOLVERUTILITIES_SOLVER_GUROBI_H_
#define SOLVERUTILITIES_SOLVER_GUROBI_H_

#include "gurobi_c.h"
#include "utilities.h"
#define ENVptr 		GRBenv
#define modelPtr 	GRBmodel
#define GRB_OFF 	0

/* The following are used to determine the status of basis (constraint or variable)
 * The values in Gurobi are negative, i.e., however, we keep the positive values to allow for compatibility with CPLEX and Gurobi.
 * 0 (basic), -1 (non-basic at lower bound), -2 (non-basic at upper bound), and -3 (super-basic)
 *
*/

enum basisStatus {
	BASIC = 0,
	AT_LOWER = 1,
	AT_UPPER = 2,
	SUPER_BASIC = 3,
};

enum modelType {
	PROB_LP,
	PROB_QP,
	PROB_QCP,
	PROB_MILP,
	PROB_MIQP,
	PROB_MIQCP
};

/* Algorithm used to solve continuous models
1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent, 5=deterministic concurrent simplex. */
enum methodType{
	AUTOMATIC = 1,
	PRIMAL_SIMPLEX = 0,
	DUAL_SIMPLEX = 1,
	BARRIER = 2,
	CONCURRENT = 3,
	DET_CONCURRENT = 4,
	DET_CONCURRENT_SIMPLEX = 5
};

/* Environment creation and destruction */
void openSolver();
void closeSolver();

/* Error Handling */
void solverErrMsg();

/* Model creation and modification */
int createProblem(char *probname, modelPtr *model);
int createProblemInit(const char *probname, modelPtr *model, int numvars, double *objx, double *lb, double *ub, char *vtype, char **varnames);
modelPtr *setupProblem(const char *Pname, int numvars, int numconstrs, int objsense, double	objcon, double	*obj, sparseMatrix *objQ, char *sense,
		double *rhs, int *vbeg, int *vlen, int *vind, double *vval, double *lb, double *ub, char *vtype, char **varnames, char **constrnames );
modelPtr *cloneProblem (modelPtr *orig);

int	changeCoefficients ( modelPtr *model, int numchgs, int *cind, int *vind, double *val );
int changeObjCoeff (modelPtr *model, int vind, double val);
int changeObjCoeffArray (modelPtr *model, int numchgs, int *vind, double *val );
int changeRHSelement (modelPtr *model, int vind, double val);
int changeRHSArray (modelPtr *model, int numchgs, int *vind, double *val );
int changeBDSelement (modelPtr *model, const char *attributeName, int vind, double val);
int changeBDSArray (modelPtr *model, const char *attributeName, int numchgs, int *vind, double *val );

int	addQPterms (modelPtr *model, int numqnz, int *qrow, int	*qcol, double *qval);
int	addRows(modelPtr *ptr, int numRows, int numnz, int *rbeg, int *rind, double	*rval, char	*senx, double *rhsx, char **rname);

int addRow(modelPtr *ptr, int numnz, double rhs, char sense, int *rmatind, double *rmatval, char *rowname);

int addCols(modelPtr *ptr, int numvars, int numnz, int *cbeg, int *cind, double *cval, double *objx, double *lb, double *ub, char *ctype, char **cname);
int	addCol (modelPtr *ptr, int numnz, int *cind, double *cval, double objx, double lb, double ub, char ctype, char *cname);
int deleteRows ( modelPtr *ptr, int numRows, int *indices );
int deleteCols ( modelPtr *ptr, int numCols, int *indices );

int freeProblem(modelPtr *model);

/* Model solution */
int solveProblem (modelPtr *model);
double getObjective (modelPtr *model);
int getPrimal ( modelPtr *model, double *X, int start, int len );
int getDual ( modelPtr *model, double *pi, int start, int len );
int getDualSlack (modelPtr *model, double *dj, int start, int len);
int getBasis(modelPtr *ptr, int *cstat, int *rstat, int mac, int mar);
int getBasisHead(modelPtr *ptr, int *head);
int getBasisInvRow(modelPtr *ptr, int idx, int len, double *phi);
int getBasisInvCol(modelPtr *ptr, int idx, int len, double *phi);

/* Model queries */
int getModelType();
int	getConstraints(modelPtr *model, int	*numnzP, int *matbeg, int *matind, double *matval, int start, int len);
int	getVariables (modelPtr *model, int *numnzP, int *vbeg, int *vind, double *vval, int start, int len);
int getObjName(char *srcFile, char **objName);
sparseMatrix *getQmatrix(modelPtr* model, int numvar);

/* Input/Output */
int readProblem(char *probpath, modelPtr **model);
int writeProblem(modelPtr *model, char *fname);

/* Attribute Management */
int getIntAttribute(modelPtr *model, const char *attributeName);
int getIntAttibuteElement (modelPtr *model, const char *attributeName, int element);
int getIntAttributeArray(modelPtr *model, const char *attributeName, int start, int len, int *attr);

int getDoubleAttribute(modelPtr *model, char *attributeName, double *attr);
int setDoubleAttribute(modelPtr *model, char *attributeName, double newvalue);
int getDoubleAttributeElement (modelPtr *model, const char *attributeName, int element, double *attr);
int setDoubleAttributeElement (modelPtr *model, const char *attributeName, int element, double newvalue);
int getDoubleAttributeArray (modelPtr *model, const char *attributeName, int start, int len, double *values);
int setDoubleAttributeArray (modelPtr *model, const char *attributename, int start, int len, double *values);

int getStringAttribute(modelPtr *model, char *attributeName, char **attr);
char *getStringAttributeElement (modelPtr *model, const char *attributeName, int element);
char **getStringAttributeArray (modelPtr *model, const char *attributeName, int start, int len);

char getCharAttrbuteElement(modelPtr *model, const char *attributeName, int element);
char *getCharAttributeArray(modelPtr *model, const char *attributeName, int start, int len );

/* Parameter Management and Tuning */
int setIntParam ( const char *paramname, double newvalue );
int setDoubleParam ( const char *paramname, double newvalue );

/* Monitoring Progress - Logging and Callbacks */

/* Modifying Solver Behavior - Callbacks */

/* Batch Requests */

/* Advanced simplex routines */

#endif /* SOLVERUTILITIES_SOLVER_GUROBI_H_ */
