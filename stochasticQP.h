#include "./solverUtilities/utilities.h"
#include "solverUtilities/solver_gurobi.h"
#include "./smpsReader/smps.h"
#include "./smpsReader/prob.h"

#define WRITE_FILES
#define ALGO_CHECK
#define STOCH_CHECK

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
	int		cnt;					/* number of elements in the structure */
	double** pi;					/* value of duals(associated with equality constraints) with random elements in right-hand side */
	double** umu;					/* value of  duals (reduced costs) with random elements in right-hand side(upperbound) */
	double** lmu;					/* value of  duals (reduced costs) with random elements in right-hand side(lowerbound) */
}lambdaType;

typedef struct {
	int         state;
	double  	alpha;	           	/* scalar pi x b */
	dVector 	beta;	           	/* dVector pi x C */
} pixbCType;

//typedef struct {
//	int 			state; 			/*0 if assignd value 1 ow*/
//	double 			alpha;
//	sparseVector* 	beta;
//} lambdadeltaType;

typedef struct {
	dVector y;
	dVector pi;
	dVector lmu;
	dVector umu;
	double  mubBar;
} solnType;

typedef struct {
	int cnt;
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
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */

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
}cellType;

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

/* Source.c */
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
void freeOmegaType(omegaType* omega, bool partial);
void freeSigma(sigmaType* sigma);
void freeLambda(lambdaType* lambda);
void freeDelta(deltaType *delta , int numobs);
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

int addtoLambda(lambdaType* lambda, solnType*dual, int numRows, int numCols, bool *newLambdaFlag);
void addtoSigma(cellType* cell, probType* prob, solnType *soln);
void addtoDelta(cellType* cell, probType* prob, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* ybar, sparseVector* yund, int obs,int num);

solnType* buildDual (numType *num);
void VsumVsparse(dVector result , dVector v, sparseVector* vs , int len);
int stocUpdateQP(cellType* cell, probType* prob, solnType* dual, sparseMatrix* COmega, sparseVector* bOmega, sparseVector* uOmega, sparseVector* lOmega);
