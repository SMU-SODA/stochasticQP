#include "stochasticQP.h"
#include "./smpsReader/smps.h"

cString outputDir;
long MEM_USED;

int main(int argc, char *argv[]) {
	oneProblem* orig;
	timeType* tim;
	stocType* stoc;
	cString inputDir , probname;

	/* Obtain parameter input from the command line */
	parseCmdLine(argc, argv, &probname, &inputDir); /*what are the inputs?*/

	/* Setup a solver environment */
	openSolver();

	/* read the problem */
	if(readFiles(inputDir, probname, &orig, &tim, &stoc)){
		errMsg("read", "main", "failed to read problem main file", 0);
		return 1;
	}

	/* close the solver environment */

	closeSolver();

	return 0;
} /*END main()*/

void parseCmdLine(int argc, char* argv[], cString* probName, cString* inputDir) {

	if ( argc == 1 ) {
		printHelpMenu(); exit(0);
	}

	for (int i = 1; (i < argc); i++) {
		if (argv[i][0] == '-') {
			switch ((argv[i])[1]) {
			case '?': printHelpMenu(); exit(0);
			case 'p': {
				(*probName) = (cString)arr_alloc(2 * BLOCKSIZE, char);
				strcpy((*probName), argv[++i]); break;
			}
			case 'i': {
				(*inputDir) = (char*)arr_alloc(2 * BLOCKSIZE, char);
				strcpy((*inputDir), argv[++i]); break;
			}
			case 'o': {
				outputDir = (char*)arr_alloc(2 * BLOCKSIZE, char);
				strcpy(outputDir, argv[++i]); break;
			}
			}
		}
		else {
			printf("Input options must begin with a '-'. Use '-?' for help.\n"); exit(0);
			printHelpMenu(); exit(0);
		}

	}
}//END printCmdLine()

void printHelpMenu() {

	printf("Required input to the program.\n");
	printf("         -p string  -> problem name.\n");
	printf("         -i string  -> input directory where the problem SMPS files are saved.\n");
	printf("         -o string  -> output directory where the result files will be written.\n");

}//END printHelpMenu()


///* int readFile(cString inputDir, cString probname ) 
//{
//	char Fname[BLOCKSIZE];
//	int       error = 0;
//
//	modelPtr *model = NULL;  /*to define a variable from a str*/
//	double obj;
//	oneProblem* prb = NULL;
//
//	/*Finding the path to the problem*/
//
//	strcpy(Fname, inputDir);
//	strcat(Fname, probname);
//	strcat(Fname, ".lp");
//
//	/* Create environment */
//
//	openSolver();  /*why are we openning the solver?*/
//
//
//	/*read the problem from the file*/
//
//	error = readProblem(Fname , &model); /*GRBreadmodel gets gurobi model so why do we use another type like
//										 prtmodel?*/
//	if (error) {
//		goto QUIT;
//	}
//
//    /*extract the information of the problem*/
//
//	prb = (oneProblem*)mem_malloc(sizeof(oneProblem)); /*what is the meaning?*/
//
//	/*extract problem size*/
//
//	prb->cols = getIntAttribute( model, "NumVars");
//
//	prb->rows = getIntAttribute(model, "NumConstrs");
//
//	/*extract linear objective coefficients*/
//
//	prb->objx = (dVector)arr_alloc(prb->cols,double); /*how memaloc and */
//
//	error = getDoubleAttributeArray(model, "Obj" , 0 , prb->cols , prb->objx);
//
//	if (error) {
//		goto QUIT;
//	}
//
//	/*extract Q matrix*/
//
//
//	prb->objQ = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix)); /*why do we write it?
//																 */
//
//	prb->objQ->col = (iVector)arr_alloc(prb->cols * prb->cols, int);
//
//	prb->objQ->row = (iVector)arr_alloc(prb->cols * prb->cols, int);
//
//	prb->objQ->val = (dVector)arr_alloc(prb->cols * prb->cols, double);
//		
//	GRBgetq(model,		&prb->objQ->cnt,		prb->objQ->row,		prb->objQ->col,		prb->objQ->val);
//
//	prb->objQ->col = (iVector)mem_realloc(prb->objQ->col, prb->objQ->cnt * sizeof(int));
//	prb->objQ->row = (iVector)mem_realloc(prb->objQ->row, prb->objQ->cnt * sizeof(int));
//	prb->objQ->val = (iVector)mem_realloc(prb->objQ->val, prb->objQ->cnt * sizeof(double));
//
//	/*extract the right-hand side*/
//
//	prb->rhsx = (dVector)arr_alloc(prb->rows, double);
//	error = getDoubleAttributeArray(model, "RHS", 0, prb->rows, prb->rhsx);
//
//	if (error) {
//		goto QUIT;
//	}
//
//
//	/*extract matrix A*/
//
//	prb->consA = (sparseMatrix*)mem_malloc(sizeof(sparseMatrix));
//
//	prb->consA->col = (iVector)arr_alloc(prb->cols , int);
//
//	prb->consA->row = (iVector)arr_alloc(prb->rows , int);
//
//	prb->consA->val = (dVector)arr_alloc(prb->cols * prb->cols, double);
//
//	/*solve the problem*/
//
//	error = solveProblem(model);
//	if (error) {
//		goto QUIT;
//	}
//
//	obj = getObjective( model);
//
//	printf("The optimal value is equal to %lf", obj);
//
//QUIT:
//
//	/* Error reporting */
//
//	if (error) {
//		
//		solverErrMsg();
//		exit(1);
//	}
//
//	/*close envirinment*/
//
//	closeSolver();
//	
//	return 0;
//} */
