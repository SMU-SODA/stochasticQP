/*
 * solver_gurobi.c
 *
 *  Created on: Dec 29, 2020
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "../solverUtilities/solver_gurobi.h"
#include "../solverUtilities/utilities.h"

ENVptr *env;
extern cString outputDir;

/********************************************************** Environment Creation and Destruction ********************************************************/
/* Create a Gurobi solver environment. */
void openSolver() {

	if ( GRBloadenv(&env, "solverErrors.log") ) {
		solverErrMsg();
		exit(1);
	}

	setIntParam("LogToConsole", GRB_OFF);

}//END openSolver()

/* Close the Gurobi solver environment. */
void closeSolver() {

	GRBfreeenv(env);

}//END closeSolver()

/******************************************************************** Error Handling ********************************************************************/
/* Retrieve the error message associated with the most recent error that occurred in an environment. */
void solverErrMsg() {

	fprintf(stderr, "Solver Error: %s\n", GRBgeterrormsg(env));

}//END solverErrMsg()

/************************************************************ Model creation and modification ************************************************************/
/* Create a new optimization model */
int createProblem(char *probname, modelPtr *model) {
	int status;

	status = GRBnewmodel(env, &model, NULL, 0, NULL, NULL, NULL, NULL, NULL);
	if ( status )
		solverErrMsg(status);

	return status;
}// END createProblem()

/* Create a new optimization model while initializing a set of decision variables. However, no constraints are included in the initial model. */
int createProblemInit(const char *probname, modelPtr *model, int numvars, double *objx, double *lb, double *ub, char *vtype, char **varnames) {
	int status;

	status = GRBnewmodel(env, &model, probname, numvars, objx, lb, ub, vtype, varnames);
	if ( status )
		solverErrMsg(status);

	return status;
}//END createProblemInit()

/* Create a new optimization model, using the provided arguments to initialize the model data (objective function, variable bounds, constraint matrix,
etc.). The model is then ready for optimization, or for modification (e.g., addition of variables or constraints, changes to variable types or bounds,
etc.). */
modelPtr *setupProblem(const char *Pname, int numvars, int numconstrs, int objsense, double	objcon,  double	*obj, char	*sense, double *rhs, int *vbeg,
		int *vlen, int *vind, double *vval, double *lb, double *ub, char *vtype, char **varnames, char **constrnames ) {
	modelPtr *model;
	int	status;

	status = GRBloadmodel (env, &model, Pname, numvars, numconstrs, objsense, objcon, obj, sense, rhs, vbeg, vlen, vind, vval, lb, ub, vtype, varnames, constrnames);
	if ( status ) {
		solverErrMsg(status);
		return NULL;
	}

	return model;
}//END setupProblem()

/* Create a copy of an existing model. */
modelPtr *cloneProblem (modelPtr *orig) {

	GRBupdatemodel(orig);
	modelPtr *copy = GRBcopymodel(orig);

	return copy;
}//END cloneProblem()

/* Change a set of constraint matrix coefficients. This routine can be used to set a non-zero coefficient to zero, to create a non-zero coefficient
 * where the coefficient is currently zero, or to change an existing non-zero coefficient to a new non-zero value. */
int	changeCoefficients ( modelPtr *model, int numchgs, int *cind, int *vind, double *val ) {
	int status;

	status = GRBchgcoeffs(model, numchgs, cind, vind, val);
	if ( status )
		solverErrMsg(status);

	return status;
}//END changeCoefficients()

/* Change the objective coefficients. */
int changeObjCoeff (modelPtr *model, int vind, double val) {
	int status;

	status = setDoubleAttributeElement(model, "Obj", vind, val);
	if ( status )
		solverErrMsg(status);

	return status;
}//changeObjCoeff()

int changeObjCoeffArray (modelPtr *model, int numchgs, int *vind, double *val ) {
	int status;

	for ( int n = 0; n < numchgs; n++ ) {
		status = setDoubleAttributeElement(model, "Obj", vind[n], val[n]);
		if ( status )
			solverErrMsg(status);
	}

	return status;
}//END changeObjCoeffArray()

/* Change the right-hand side values. */
int changeRHSvalue (modelPtr *model, int vind, double val) {
	int status;

	status = setDoubleAttributeElement(model, "RHS", vind, val);
	if ( status )
		solverErrMsg(status);

	return status;
}//END changeRHSvalue()

int changeRHSvalueArray (modelPtr *model, int numchgs, int *vind, double *val ) {
	int status;

	for ( int n = 0; n < numchgs; n++ ) {
		status = setDoubleAttributeElement(model, "RHS", vind[n], val[n]);
		if ( status )
			solverErrMsg(status);
	}

	return status;
}//END changeRHSvalueArray()

/* Free a model and release the associated memory. */
int freeProblem(modelPtr *model) {
	int status;

	status = GRBfreemodel(model);
	if (status)
		solverErrMsg(status);

	return status;
}//END freeProblem()

/********************************************************************** Model Solution ********************************************************************/
/* Optimize a model. The algorithm used for the optimization depends on the model type (simplex or barrier for a continuous model; branch-and-cut for
 * a MIP model). Upon successful completion, this method will populate the solution related attributes of the model. */
int solveProblem ( modelPtr *model ) {
	int status, optimstatus;

	status =  GRBoptimize(model);
	if ( status ) {
		solverErrMsg(status);
		return 1;
	}

	/* Capture solution information */
	status = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
	if ( status )
		return 1;

	/* If model is infeasible or unbounded, turn off presolve and resolve */
	if (optimstatus == GRB_INF_OR_UNBD) {
		/* Change parameter on model environment.  The model now has a copy of the original environment, so changing the original will no longer affect
		 * the model.  */
		status = GRBsetintparam(GRBgetenv(model), "PRESOLVE", 0);
		if ( status ) {
			solverErrMsg(status);
			return 1;
		}

		status = GRBoptimize(model);
		if ( status ) {
			solverErrMsg(status);
			return 1;
		}

		status = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
		if ( status ) {
			solverErrMsg(status);
			return 1;
		}
	}

	if (optimstatus == GRB_OPTIMAL) {
		goto TERMINATE;
	} else if (optimstatus == GRB_INFEASIBLE) {
		printf("Model is infeasible\n\n");

		status = GRBcomputeIIS(model);
		if ( status ) {
			solverErrMsg(status);
			return 1;
		}

		status = GRBwrite(model, "model.ilp");
		if ( status ) {
			solverErrMsg(status);
			return 1;
		}
	} else if (optimstatus == GRB_UNBOUNDED) {
		printf("Model is unbounded\n\n");
		return 1;
	} else {
		printf("Optimization was stopped with status = %d\n\n", optimstatus);
		return 1;
	}

	TERMINATE:
	return 0;
}//END solveProblem()

/* Obtain the optimal objective function value. */
double getObjective ( modelPtr *model ) {
    int status;
    double objval;

    status = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if ( status ) {
    	solverErrMsg(status);
    }

    return objval;
}//END getObjective()

/* Obtain the primal solution vector. */
int getPrimal ( modelPtr *model, double *X, int start, int len ) {
	int status;

	status = getDoubleAttributeArray(model, "X", start, len, X);
	if ( status )
		errMsg("solver", "getPrimal", "failed to retrieve the primal solution", 0);

	return status;
}//END getPrimal()

/* Obtain the dual solution vector. */
int getDual ( modelPtr *model, double *pi, int start, int len ) {
	int status;

	status = getDoubleAttributeArray(model, "Pi", start, len, pi+1);
	if ( status )
		errMsg("solver", "getPrimal", "failed to retrieve the primal solution", 0);

	pi[0]= oneNorm(pi+1, len);

	return status;
}//END getPrimal()

/********************************************************************** Model queries *********************************************************************/
int getModelType(modelPtr *model) {
	int status, IsMIP, IsQP, IsQCP, attr = 0, type;

	status = GRBgetintattr(model, "IsMIP", &IsMIP);
	if ( status ) {
		solverErrMsg(status);
	}
	attr = attr << IsMIP;

	status = GRBgetintattr(model, "IsQP", &IsQP);
	if ( status ) {
		solverErrMsg(status);
	}
	attr = attr << IsQP;

	status = GRBgetintattr(model, "IsQCP", &IsQCP);
	if ( status ) {
		solverErrMsg(status);
	}
	attr = attr << IsQCP;

	switch ( attr ) {
	case 0: type = PROB_LP; break;
	case 1: type = PROB_QCP; break;
	case 2: type = PROB_QP; break;
	case 4: type = PROB_MILP; break;
	case 5: type = PROB_MIQCP; break;
	case 6: type = PROB_MIQP; break;
	default: break;
	}

	return type;
}// END getModelType()

int	getConstraints(modelPtr *model, int	*numnzP, iVector cbeg, iVector cind, dVector cval, int start, int len) {
	int status;

	status = GRBgetconstrs(model, numnzP, cbeg, cind, cval, start, len);
	if ( status ) {
		solverErrMsg(status);
	}

	return status;
}//END getConstraints()

int	getVariables (modelPtr *model, int *numnzP, int *vbeg, int *vind, double *vval, int start, int len) {
	int status;

	status = GRBgetvars (model, numnzP, vbeg, vind, vval, start, len);
	if ( status ) {
		solverErrMsg(status);
	}

	return status;
}//END getVariables()

/********************************************************************** Input/Output **********************************************************************/
/* Read a model from a file. */
int readProblem(cString probpath, modelPtr **model) {
	int status;

	status = GRBreadmodel(env, probpath, model);
	if( status ) {
		solverErrMsg(status);
	}

	return status;
}// END readProblem()

int writeProblem(modelPtr *model, cString fname) {
	int status;
	char buffer[2*BLOCKSIZE];

	strcpy(buffer,outputDir);
	strcat(buffer,fname);

	status = GRBwrite(model, buffer);
	if (status)
		solverErrMsg(status);

	return status;
}//END writeProblem()

/******************************************************************* Attribute Management ******************************************************************/
/* Query the value of an integer-valued model attribute. */
int getIntAttribute(modelPtr *model, const char *attributeName) {
	int status, attr;

	status = GRBgetintattr(model, attributeName, &attr);
	if ( status )
		solverErrMsg(status);

	return attr;
}//END getIntAttribute()

/* Set the value of an integer-valued model attribute. */
int	setIntAttribute (GRBmodel *model, const char *attrname, int newvalue) {
	int status;

	status =  GRBsetintattr (model, attrname, newvalue);
	if ( status )
		solverErrMsg(status);

	return status;
}//END setIntAttribute()


/* Query a single value from an integer-valued array attribute. */
int getIntAttibuteElement (modelPtr *model, const char *attributeName, int element) {
	int status, valueP;

	status = GRBgetintattrelement(model, attributeName, element, &valueP);
	if ( status )
		solverErrMsg(status);

	return valueP;
}//END getIntAttibuteElement()

/* Query the values of an integer-valued array attribute. */
int *getIntAttributeArray(modelPtr *model, const char *attributeName, int start, int len) {
	int status;
	int *attr;

	attr = (int *) arr_alloc(len, int);

	status = GRBgetintattrarray(model, attributeName, start, len, attr);
	if ( status ) {
		solverErrMsg(status);
		return NULL;
	}

	return attr;
}//END getIntAttributeArray()

/* Query the value of a double-valued model attribute. */
int getDoubleAttribute(modelPtr *model, char *attributeName, double *attr) {
	int status;

	status = GRBgetdblattr(model, attributeName, attr);
	if ( status )
		solverErrMsg(status);

	return status;
}//END getDoubleAttribute()

/* Set the value of a double-valued model attribute. */
int setDoubleAttribute(modelPtr *model, char *attributeName, double newvalue) {
	int status;

	status = GRBsetdblattr(model, attributeName, newvalue);
	if ( status )
		solverErrMsg(status);

	return status;
}//END getDoubleAttribute()

/* Query a single value from a double-valued array attribute. */
int getDoubleAttributeElement (modelPtr *model, const char *attributeName, int element, double *attr) {
	int status;

	status = GRBgetdblattrelement(model, attributeName, element, attr);
	if ( status )
		solverErrMsg(status);

	return status;
}//END getDoubleAttibuteElement()

int setDoubleAttributeElement (modelPtr *model, const char *attributeName, int element, double newvalue) {
	int status;

	status = GRBsetdblattrelement(model, attributeName, element, newvalue);
	if ( status )
		solverErrMsg(status);

	return status;
}//END getDoubleAttibuteElement()

/* Query the values of a double-valued array attribute. */
int getDoubleAttributeArray (modelPtr *model, const char *attributeName, int start, int len, double *values) {
	int status;

	status = GRBgetdblattrarray(model, attributeName, start, len, values);
	if ( status ) {
		solverErrMsg(status);
	}

	return status;
}//END getIntAttributeArray()

/* Set the values of a double-valued array attribute. */
int setDoubleAttributeArray (modelPtr *model, const char *attributename, int start, int len, double *values) {
	int status;

	status = GRBsetdblattrarray (model, attributename, start, len, values);
	if ( status ) {
		solverErrMsg(status);
	}

	return status;
}//END setDoubleAttributeArray()

/* Query the value of a string-valued model attribute. */
char *getStringAttribute(modelPtr *model, char *attributeName) {
	int status;
	char *attr;

	status = GRBgetstrattr(model, attributeName, &attr);
	if ( status )
		solverErrMsg(status);

	return attr;
}//END getDoubleAttribute()

/* Query a single value from a string-valued array attribute. */
char *getStringAttributeElement (modelPtr *model, const char *attributeName, int element) {
	int status;
	char *attr;

	status = GRBgetstrattrelement(model, attributeName, element, &attr);
	if ( status )
		solverErrMsg(status);

	return attr;
}//END getCharAttibuteElement()

/* Query the values of a string-valued array attribute. */
char **getStringAttributeArray (modelPtr *model, const char *attributeName, int start, int len) {
	int status;
	char **attr;

	attr = (char **) arr_alloc(len, char *);

	status = GRBgetstrattrarray(model, attributeName, start, len, attr);
	if ( status ) {
		solverErrMsg(status);
		return NULL;
	}

	return attr;
}//END getStringAttributeArray()

/* Query a single value from a character-valued array attribute. */
char getCharAttrbuteElement(modelPtr *model, const char *attributeName, int element) {
	int status;
	char attr;

	status = GRBgetcharattrelement(model, attributeName, element, &attr);
	if ( status )
		solverErrMsg(status);

	return attr;
}//END getCharAttrbuteElement()

/* Query the values of a character-valued array attribute. */
char *getCharAttributeArray(modelPtr *model, const char *attributeName, int start, int len ) {
	int status;
	char *attr;

	attr = (char *) arr_alloc(len, char);

	status = GRBgetcharattrarray(model, attributeName, start, len, attr);
	if ( status ) {
		solverErrMsg(status);
		return NULL;
	}

	return attr;
}//END GRBgetcharattrarray()

/******************************************************************* Parameter Management and Tuning ******************************************************************/
/* Modify the value of an integer-valued parameter. */
int setIntParam ( const char *paramname, double newvalue ) {
	int status;

	status = GRBsetintparam (env, paramname, newvalue);
	if (status)
		solverErrMsg(status);

	return status;
}//END setIntParam()

/* Modify the value of a double-valued parameter. */
int setDoubleParam ( const char *paramname, double newvalue ) {
	int status;

	status = GRBsetdblparam (env, paramname, newvalue);
	if (status)
		solverErrMsg(status);

	return status;
}//END setDoubleParam()
