#include "./solverUtilities/utilities.h"
#include "solverUtilities/solver_gurobi.h"
#include "./smpsReader/smps.h"
#include "./smpsReader/prob.h"


#undef WRITE_FILES
#undef ALGO_CHECK
#define STOCH_CHECK


double* testLA (double* A , int S);

typedef enum {
	FULL,
	DUALLBASED,
	PARTITIONBASED
} algoType;

typedef struct {
	int		numRV;					/* Number of random variables */
	int 	cnt;					/* Number of observations */
	dVector	probs;					/* Probability of observation */
	iVector weights;				/* Number of times a particular observation is encountered (not used, but included to align with sp
									   algorithms repository subroutines). */
	dVector* vals;					/* Observation values */
} omegaType;


typedef struct {
	dVector entries;
	int row;
	int col;
}Mat;
Mat* pdas(sparseMatrix* D, int* partition, sparseMatrix* Q, dVector coef , dVector Rhs  , int cols , int rows , dVector Ubound , dVector Lbound);
typedef struct {
	int		cnt;					/* number of elements in the structure */
	double* 	mubar;					
	double**  y;					/* value of duals the section equal to primal */
	double** pi;					/* value of duals(associated with equality constraints) with random elements in right-hand side */
	double** umu;					/* value of  duals (reduced costs) with random elements in right-hand side(upperbound) */
	double** lmu;					/* value of  duals (reduced costs) with random elements in right-hand side(lowerbound) */
	Mat** pd;
}lambdaType;

typedef struct {
	int         feas ;
	double  	alpha;	           	/* scalar pi x b */
	dVector 	beta;	           	/* dVector pi x C */
} pixbCType;

typedef struct {
	dVector y;
	dVector pi;
	dVector lmu;
	dVector umu;
	double  mubBar;
} solnType;

typedef struct {
	int cnt;
	double*** dy;
	double*** dmu;
	double*** dnu;
	pixbCType ***vals;
} deltaType;

typedef struct {
	double	repTime;
	double 	iterTime;
	double 	masterIter;
	double 	subprobIter;
	double 	optTestIter;
	double  argmaxIter;
	double 	iterAccumTime;
	double 	masterAccumTime;
	double 	subprobAccumTime;
	double 	optTestAccumTime;
	double  argmaxAccumTime;
	double	reduceTime;
}runTime;

typedef struct {
	int 		cnt;				/* Number of elements */
	pixbCType** vals;				/* product terms */
} sigmaType;

typedef struct {
	int         cnt;        /* Number of elements */
	iVector*	part;       /* Storing partitions with 0 inactive 1 lower bound, and 2 upperbound*/
	long long int* basnum;	/* partition ID encoded using the vector part[i] */
	iVector		low;
	iVector 	up;
	iVector 	inact;
}PartitionType;

/* structure for the problem type:
 * c_t^\top x_t + \min  d_t^\top u_t + \expect{h_{t+}(s_{t+})}
 *                 s.t. D_t u_t = b_t - C_tx_t
 * where, x_{t+} = a_{t+} + A_{t+}x_t + B_{t+}u_t.
 */
typedef struct {
	int		ck;					/* Iteration when the cut was generated */
	double  alpha;              /* scalar value for the right-hand side */
	dVector  beta;               /* coefficients of the master problems's primal variables */

	// bool	isIncumb;			/* indicates if the cut is an incumbent cut */
	// double alphaIncumb;		/* right-hand side when using QP master, this is useful for quick updates */

	int 	rowNum;				/* row number for master problem in solver */
	// int	omegaID;			/* the observation ID used when multi-cut option is used */
	// iVector iStar;				/* Holds the ID for the sigma which is associated with each observation, dual index id, which dual we  */
	// int 	form;				/* determines the form of the cut (l-shaped regular, l-shaped callback, MIR, GMI */
	cString	name;
}oneCut;


typedef struct {
	int    	cnt;                    /* number of cuts */
	oneCut** vals;					/* values which define the set of cuts */
}cutsType;




typedef struct {
	int cnt;
	solnType*** sol;
}deltaSolType;



typedef struct {
	int cnt;
	Mat** wt ;
}WT;


typedef struct {
	int cnt;
	solnType** vals;
}solutionSetType;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved.*/

	oneProblem* master;             /* store master information */
	oneProblem* subprob;			/* store subproblem information */
	dVector      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */
	dVector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */

	//double		gamma;				/* Improvement in obejctive function */
	//double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	//bool        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	//iVector     iCutIdx;			/* index of incumbent cuts in cell->cuts structure. If multicut is used, there will be one//							   for each observation. */

	dVector		piM;
	int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType*   cuts;              /* optimality cuts */
	cutsType*   fCuts;             /* feasibility cuts */
     
	omegaType*  omega;				/* all realizations observed during the algorithm */

	bool        optFlag;
	bool		spFeasFlag;			/* Indicates whether the subproblem is feasible */
	int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	bool		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */

	runTime* time;				/* Run time structure */
	lambdaType* lambda;
	sigmaType* sigma;
	deltaType* delta;
	PartitionType* partition;
	deltaSolType* deltaSol;
	WT* wtSet;

	int numit;
	double obj;
	double Tcut; /* Total time required to produce cuts */
	double Tmas; /* Total time required to solve the master problem */
	double Tsub; /* Total time for solving subpoblems*/
	double  Totaltime; /* total time */
	int IterPart; /* the iteration that all the partitions are recognized  */
	double stochupdate; /* average time of iterations */

}cellType;


void newWTset(int structSize, cellType* cell);
void freeWTset(cellType* cell);





typedef struct {

	long long* RUN_SEED;		/* seed used during optimization */
	long long* EVAL_SEED;		/* seed used during evaluation */
	int 	NUM_SEEDS;			/* Number of run seeds provided */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTER_TYPE;		/* type of master problem */
	double	EPSILON;			/* Optimality gap */
	int		MULTICUT;			/* Set to 1 if multicut is to be solved */
	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	int 	MULTIPLE_REP;		/* When multiple replications are needed, set this to an integer number of replications */
	int		ALGOTYPE;			/*0 we do not use the duals and we solve all sebproblems  , 1 reuse the partitions,2 when we use the duals*/
	double	SAMPLE_FRACTION;	/* A fraction (0,1] to determine what fraction of samples are solved */


	int		SAA; 				/* Use SAA when continuous distribution in stoch file (1), or not (0) */
	int		MAX_OBS;			/* Maximum number of iterations before which SAA is invoked */

	double  MAX_TIME;			/* Maximum per replication run time */
}configType;

typedef struct {
	int				ck;			/* The first time the basis was encountered. */
	int				weight;		/* Frequency of observation for each unique basis */
	unsigned long* rCode;		/* Encoded row status in the basis (currently not being used */
	unsigned long* cCode;		/* Encoded column status in the basis */
	int				phiLength;	/* Number of basic columns with random cost coefficients */
	dVector* phi;		/* The phi matrix: the columns of inverse dual basis matrix which have random cost coefficients */
	iVector			omegaIdx;	/* Indices within the random cost coefficient dVector to which the columns of phi matrix correspond to. */
	iVector			sigmaIdx;	/* Indices within the random cost coefficient dVector to which the columns of phi matrix correspond to. */
	dVector			piDet;		/* Deterministic component of the dual solution. This depends only on the basis. */
	double			mubBar;
	dVector			gBar;
	sparseMatrix* psi;		/* The simplex tableau matrix corresponding to the basis. */
	bool			feasFlag;
}oneBasis;

/* The basis type data structure holds all the information regarding the basis identified during the course of the algorithm.
 * This structure will be at the heart of all calculations related to stochastic updates. */
typedef struct {
	int			basisDim;	/* The dimension of the basis matrix */
	int			cnt;		/* Number of unique basis encountered by the algorithm */
	int			rCodeLen;	/* Length of encoded row status */
	int			cCodeLen;	/* Length of encoded column status */
	bool** obsFeasible;
	oneBasis** vals;		/* a structure for each basis */
}basisType;




/* Subroutine stochasticQP_main.c */
void parseCmdLine(int argc, char* argv[], cString* probName, cString* inputDir);
void printHelpMenu();

/* fullSolve.c */
oneCut *fullSolveCut(probType *prob, cellType* cell, stocType* stoch, double* x);

/* dualSolve.c */
int argmax(probType *prob, sigmaType *sigma, deltaType *delta, dVector Xvect, int obs);
int argmaxPart(cellType* cell, double tol, probType* prob, sigmaType* sigma, deltaType* delta, dVector Xvect, int obs, int prev, int* index, bool* flag, pixbCType** deltaxp, double** dy, double** dld, double** dnu, double** dmu, sparseVector* bomega, sparseVector* lomega, sparseVector* uomega);

cellType* buildCell(probType** prob, stocType* stoc);
int solveSubprob(probType* prob, oneProblem* subproblem, dVector Xvect, dVector obsVals,
		sparseVector* bOmega, sparseMatrix* COmega, sparseVector* dOmega, sparseVector* lOmega, sparseVector* uOmega, solnType *dual);

dVector computeRHS(sparseVector *bBar, sparseMatrix *Cbar, sparseVector* bOmega, sparseMatrix* COmega, dVector X, int numRows);
dVector computeCostCoeff(sparseVector *dBar, sparseVector* dOmega, int numCols);
dVector computeBDS(sparseVector* bdsBar, sparseVector* bdsOmega, int numCols);

int getBoundDual(modelPtr *model, int numCols, double* mu_up, double* mu_low);
int computeMUdual(modelPtr* model, int numCols, double* dj);

void buildOmegaCoordinates (probType *prob, sparseVector bOmega, sparseMatrix COmega, sparseVector dOmega, sparseVector uOmega, sparseVector lOmega);

int readConfig(cString configFile);
void freeConfig();
oneProblem* setRhs(oneProblem* subProb, dVector rhs);
oneProblem *newMaster(oneProblem *probSP);
oneProblem* newSubproblem(oneProblem* probSP);

int chgObjxwObserv(modelPtr* lp, numType* num, coordType* coord, dVector cost, iVector indices, dVector observ);
int chgRHSwObserv(modelPtr* lp, numType* num, coordType* coord, dVector observ, dVector spRHS, dVector X);

/* Subroutines in algo.c */
int addCut2Solver(modelPtr *model, oneCut *cut, int lenX);
oneCut *newCut(int numX);
int runAlgo(probType** prob, stocType* stoc, cellType* cell);
int updateRHSwState(numType* num, coordType* coord, sparseVector* bBar, sparseMatrix* Cbar, dVector X,
		dVector obs, dVector *rhs);
void freeCellType(cellType* cell);
void freecut(cutsType* cut);
void freeonecut(oneCut* cut);
void freeLambdaDelta(pixbCType* lambdadelta);
void freeOmegaType(omegaType* omega, bool partial);
void freeSigma(sigmaType* sigma);
void freeLambda(lambdaType* lambda);
void freeDelta(deltaType *delta , int numobs);
void freePartition(PartitionType* partition);
oneCut* dualSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset);
int solveSubprobdual(probType* prob, oneProblem* subproblem, dVector Xvect, dVector obsVals, dVector piS, double*,double* mu2, double* mu3);
int calcSigma(sigmaType* sigma, cellType* cell  ,probType** prob, dVector pi, dVector mu2, dVector mu3 , sparseVector* bOmega, sparseMatrix* COmega,
		sparseVector* yuOmega , int obs);
int stochasticUpdates(probType** prob, cellType* cell, stocType* stoch, lambdaType* lambda,sigmaType* sigma  ,double* x ,int rand);
sigmaType* newSigma(double SigmaSize, probType** prob);
lambdaType* newLambda(double SigmaSize, probType** prob);
deltaType* newDelta(double SigmaSize, probType** prob, cellType* cell);
void freeDelta(deltaType* delta, int numobs);
void freeLambda(lambdaType* lambda);
void freeSigma(sigmaType* sigma);
void sample(int* omegaP, int numsample, int numobs);
void  subtractSample(int * omegaP, int* omegaQ , int numobs, int numsample);
void  SampleGen(int* omegaP, int* omegaQ, int cnt, int solveset);
bool *subsetGenerator(int numObs);

omegaType* newOmega(stocType* stoc);

void addtoLambdaP(cellType* cell, solnType* soln, Mat* W, probType* prob, sparseVector* bOmega, sparseVector* uOmega,
	sparseVector* lOmega, int low, int up, int inact, dVector dx);


 int addtoLambda(lambdaType* lambda, solnType* dual, int numRows, int numCols, bool* newLambdaFlag);
void addtoSigma(cellType* cell, probType* prob, solnType *soln);
void addtoDelta(cellType* cell, probType* prob, sparseMatrix* COmega, sparseVector* bOmega,
		sparseVector* ybar, sparseVector* yund, int obs,int num);

solnType* buildSolnType (numType *num);
void freeSolnType(solnType *soln);

void VsumVsparse(dVector result , dVector v, sparseVector* vs , int len);
int stocUpdateQP(cellType* cell, probType* prob, solnType* dual, sparseMatrix* COmega, sparseVector* bOmega,
		sparseVector* uOmega, sparseVector* lOmega);
void PartCalc(solnType* sol, dVector yund, dVector ybar, int numc, int* part, int* up, int* inact, int* low);
PartitionType *newPartition(int Partsize);
oneCut * partSolve(probType* prob, cellType* cell, stocType* stoch, double* x, double solveset);
int addtoPartition(probType* prob, cellType* cell, sparseVector* uOmega, sparseVector* lOmega, solnType* soln,
	bool* flag, int* up, int* inact, int* low,  dVector lStat, dVector uStat);
void newSolSet(int Partsize, probType* prob, cellType* cell);
void freeSolSet(solutionSetType* SolSet);
void Buildbase(long long int* basis, int cols, int bas);

Mat* CombineWT(int rows , int cols, Mat* W, Mat* T, int low, int up, int inact);

//typedef struct MatList {
	//Mat* mat;
	//MatList* next;
//}MatList;


Mat* newmat(int r, int c, double d);
double det(Mat* M);
Mat* transSparsM(sparseMatrix* M, int col, int row);
void removeRow(Mat* A, int r);
void removeCol(Mat* A, int c);
Mat *shrinkMat_Col(Mat *A, int colIdx);
Mat* shrinkMat_Row(Mat *A, int rowIdx);
Mat* transpose(Mat* A);
Mat* inverse(Mat* A);
Mat* scalermultiply(Mat* M, double c);
Mat* multiply(Mat* A, Mat* B);
Mat* sum(Mat* A, Mat* B);
void freemat(Mat* A);
sparseMatrix* BuildHess(sparseMatrix* M);
void CalC(int* part, sparseMatrix* Q, sparseMatrix* D , Mat** W, Mat** T, int low, int up, int inact , int rows , int cols );
Mat* transSparsM(sparseMatrix* M, int col, int row);
Mat* removerow(Mat* A, int r);
Mat* removecol(Mat* A, int c);
void AddtoSigmaP(cellType* cell, solnType* sol, probType* prob);

void addtoDeltaP(cellType* cell, Mat* W, Mat* T, Mat* WT, probType* prob, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* uOmega, sparseVector* lOmega, int obs, int lambdaIdx, int inact, int up, int low, dVector lStat, dVector uStat, double* dy, double* dld, double* ddnu, double* ddmu);

void newDeltaSol(cellType* cell, int sigmaSize, int obsnum);
Mat* adjoint(Mat* A);
void removecol2(Mat* A, Mat* B, int c);
void showmat(Mat* A);
int StocUpdatePart(cellType* cell, probType* prob, sparseVector* bOmega, sparseMatrix* COmega, sparseVector* lOmega,
	sparseVector* uOmega, solnType* soln,  int* partIndx, dVector dx, pixbCType** deltax, double** dy, double** dld, double** dnu, double** dmu);

void AddtoDettaX(probType* prob, cellType* cell, pixbCType** delta, dVector deltaX, int low, int up, int inact, int partindex, double** deltay, double** deltaLd, double** dnu, double** dmu);
