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

#define ENVptr 		GRBenv
#define modelPtr 	GRBmodel   /*Changing the name of the structure to modify it for different solvers*/
#define GRB_OFF 	0

enum modelType {
	PROB_LP,
	PROB_QP,
	PROB_QCP,
	PROB_MILP,
	PROB_MIQP,
	PROB_MIQCP
};

/* Environment creation and destruction */
void openSolver();
void closeSolver();

/* Error Handling */
void solverErrMsg();

/* Model creation and modification */
int createProblem(char *probname, modelPtr *model);
int createProblemInit(const char *probname, modelPtr *model, int numvars, double *objx, double *lb, double *ub, char *vtype, char **varnames);
modelPtr *setupProblem(const char *Pname, int numvars, int numconstrs, int objsense, double	objcon,  double	*obj, char	*sense, double *rhs, int *vbeg,
		int *vlen, int *vind, double *vval, double *lb, double *ub, char *vtype, char **varnames, char **constrname);
modelPtr *cloneProblem (modelPtr *orig);

int	changeCoefficients ( modelPtr *model, int numchgs, int *cind, int *vind, double *val );
int changeObjCoeff (modelPtr *model, int vind, double val);
int changeObjCoeffArray (modelPtr *model, int numchgs, int *vind, double *val );
int changeRHSvalue (modelPtr *model, int vind, double val);
int changeRHSvalueArray (modelPtr *model, int numchgs, int *vind, double *val );

int freeProblem(modelPtr *model);

/* Model solution */
int solveProblem (modelPtr *model);
double getObjective (modelPtr *model);
int getPrimal ( modelPtr *model, double *X, int start, int len );
int getDual ( modelPtr *model, double *pi, int start, int len );

/* Model queries */
int getModelType();
int	getConstraints(modelPtr *model, int	*numnzP, int *matbeg, int *matind, double *matval, int start, int len);
int	getVariables (modelPtr *model, int *numnzP, int *vbeg, int *vind, double *vval, int start, int len);

/* Input/Output */
int readProblem(char *probpath, modelPtr **model); 
int writeProblem(modelPtr *model, char *fname);

/* Attribute Management */
int getIntAttribute(modelPtr *model, const char *attributeName);
int getIntAttibuteElement (modelPtr *model, const char *attributeName, int element);
int *getIntAttributeArray(modelPtr *model, const char *attributeName, int start, int len);

int getDoubleAttribute(modelPtr *model, char *attributeName, double *attr);
int setDoubleAttribute(modelPtr *model, char *attributeName, double newvalue);
int getDoubleAttributeElement (modelPtr *model, const char *attributeName, int element, double *attr);
int setDoubleAttributeElement (modelPtr *model, const char *attributeName, int element, double newvalue);
int getDoubleAttributeArray (modelPtr *model, const char *attributeName, int start, int len, double *values);
int setDoubleAttributeArray (modelPtr *model, const char *attributename, int start, int len, double *values);

char *getStringAttribute(modelPtr *model, char *attributeName);
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
