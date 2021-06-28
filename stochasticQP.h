#include "./solverUtilities/utilities.h"
#include "solverUtilities/solver_gurobi.h"

typedef struct{
	int cols;
	int rows;
	
	dVector objx;
	dVector rhsx;
	sparseMatrix *objQ; /*why this one has star but dvector not?*/
	sparseMatrix *consA;
		}oneProblem;

void parseCmdLine(int argc, char* argv[], cString* probName, cString* inputDir);

void printHelpMenu();

int readFile(cString inputDir, cString probname);

