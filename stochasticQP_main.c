#include "stochasticQP.h"

cString outputDir;
long MEM_USED;
configType config;

int main(int argc, char* argv[]) {
	cString inputDir, probname;
	int numStages;
	stocType*  stoch= NULL;
	probType** prob = NULL;
	cellType*  cell = NULL;
	char configFile[BLOCKSIZE];

	/* Obtain parameter input from the command line */
	parseCmdLine(argc, argv, &probname, &inputDir);

	/* read algorithm configuration file */
#if _WIN64
	strcpy(configFile, "C:\\Users\\Niloofar\\source\\repos\\stochasticQP\\config.sqp");
#else
	strcpy(configFile, "./config.sqp");
#endif
	if (readConfig(configFile))
		goto TERMINATE;

	/* set up output directory: using the outputDir in config file and the input problem name */
	createOutputDir(outputDir, "stochasticQP", probname);

	/*This function reads the problem and decomposes that into stages.*/
	prob = newProbwSMPS(inputDir, probname, &stoch, &numStages);
	if ( prob == NULL ) {
		errMsg("read", "main", "failed to read files or setup the probType", 0);
		goto TERMINATE;
	}

	/*Build the algorithm cell..*/
	cell = buildCell(prob, stoch);
	if ( cell == NULL ) {
		errMsg("setup", "main", "failed to build the cell", 0);
		goto TERMINATE;
	}

	runAlgo(prob, stoch, cell);

	printf("Successfully completed executing the algorithm.\n");

	/* Free all the structures */
	if (prob) freeProbType(prob, 2);
	mem_free(probname);
	mem_free(inputDir);
	if (cell) freeCellType(cell);
	if (stoch) freeStocType(stoch);
	TERMINATE: return 0;
} /*END main()*/ 

void parseCmdLine(int argc, char* argv[], cString* probName, cString* inputDir) {

	if (argc == 1) {
		printHelpMenu(); exit(0);	}


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
			case 'a': {
				config.ALGOTYPE = atoi(argv[++i]);
				break;
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

int readConfig(cString configFile) {
	FILE* fptr;
	char	line[2 * BLOCKSIZE], comment[2 * BLOCKSIZE];
	int 	status, r2 = 0, maxReps = 30;

	fptr = fopen(configFile, "r");

	if (fptr == NULL) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	config.RUN_SEED = (long long*)arr_alloc(maxReps + 1, long long);
	config.EVAL_SEED = (long long*)arr_alloc(maxReps + 1, long long);
	config.MAX_OBS = 0;
	config.NUM_SEEDS = 0;

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED"))) {
			fscanf(fptr, "%lld", &config.RUN_SEED[config.NUM_SEEDS + 1]);
			config.NUM_SEEDS++;
			if (config.NUM_SEEDS > maxReps + 1) {
				maxReps *= 2;
				config.RUN_SEED = (long long*)mem_realloc(config.RUN_SEED, (maxReps + 1) * sizeof(long long));
			}
		}
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTER_TYPE")))
			fscanf(fptr, "%d", &config.MASTER_TYPE);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);
		else if (!(strcmp(line, "MULTICUT")))
			fscanf(fptr, "%d", &config.MULTICUT);
		else if (!(strcmp(line, "ALGOTYPE")))
			fscanf(fptr, "%d", &config.ALGOTYPE);

		else if (!(strcmp(line, "EVAL_SEED"))) {
			fscanf(fptr, "%lld", &config.EVAL_SEED[r2 + 1]);
			r2++;
			if (r2 > maxReps + 1) {
				maxReps *= 2;
				config.EVAL_SEED = (long long*)mem_realloc(config.EVAL_SEED, (maxReps + 1) * sizeof(long long));
			}
		}

		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);

		else if (!(strcmp(line, "MULTIPLE_REP")))
			fscanf(fptr, "%d", &config.MULTIPLE_REP);

		else if (!(strcmp(line, "SAA")))
			fscanf(fptr, "%d", &config.SAA);
		else if (!(strcmp(line, "MAX_OBS")))
			fscanf(fptr, "%d", &config.MAX_OBS);


		else if (!(strcmp(line, "SAMPLE_FRACTION")))
			fscanf(fptr, "%lf", &config.SAMPLE_FRACTION);


		else if (!(strcmp(line, "MAX_TIME")))
			fscanf(fptr, "%lf", &config.MAX_TIME);

		else if (!strcmp(line, "//"))
			fgets(comment, 2 * BLOCKSIZE, fptr);
		else {
			printf("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	config.NUM_SEEDS = minimum(config.NUM_SEEDS, r2);
	if (config.MULTIPLE_REP > config.NUM_SEEDS) {
		printf("Requesting to perform more replications than the number of seeds provided.\n");
		return 1;
	}

#if 0
	if (config.MULTICUT == 1) {
		printf("Version alert:: Multicut version of L-shaped is not supported in this version.\n");
		return 1;
	}
#endif

	return 0;
}//END readConfig()

void freeConfig() {
	if (config.RUN_SEED)
		mem_free(config.RUN_SEED);
	if (config.EVAL_SEED)
		mem_free(config.EVAL_SEED);
}//END freeConfig()


//* Mtrix W caculation  *//
void CalcWT(cellType* cell, probType* prob, sparseMatrix* Q, sparseMatrix* D , Mat** W, Mat** T, int low, int up, int inact) {

	int elm = 0;
	Mat* QII = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QIU = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QLU = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QUU = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QUI = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* QLI = transSparsM(Q, prob->num->cols, prob->num->cols);
	Mat* DMU = transSparsM(D, prob->num->cols, prob->num->rows);
	Mat* DML = transSparsM(D, prob->num->cols, prob->num->rows);
	Mat* DMI = transSparsM(D, prob->num->cols, prob->num->rows);
	Mat* M1;
	Mat* M2;
	Mat* minvM1;
	Mat* DMLT;
	Mat* DMUT;
	Mat* w2;
	Mat* W2W;
	Mat* DMIT;
	Mat* invM1;
	Mat* w1;

	/*Build Q(II) in a mat strcture*/
	int cnt = cell->partition->cnt - 1;
	for (int i = prob->num->cols ; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			QII = removecol(QII, i);
			QII = removerow(QII, i);
		}
	}

	/*Build D(MI)*/
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			DMI = removecol(DMI, i);
		}
	}
	
	/*Build D(MI) transpose*/

	DMIT = transpose(DMI);

	M1 = newmat(prob->num->rows + inact , prob->num->rows + inact , 0);
	
	//* Build Matrix M1 which is equal to [ QII , DMIT ; DMI , 0 ] *//

	// 1. place the first I rows

	for (int i = 0; i < inact; i++) {
		for (int j = 0; j < inact ; j++) {
			M1->entries[elm] = QII->entries[i* inact + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			M1->entries[elm] = DMIT->entries[i* prob->num->rows +j];
			elm++;
		}
	}

	// 2. place the next M rows
	for (int i = 0; i < prob->num->rows; i++) {
		for (int j = 0; j < inact; j++) {
			M1->entries[elm] = DMI->entries[(i) * inact + j];
			elm++;
		}
		for (int j = 1; j <= prob->num->rows; j++) {
			M1->entries[elm] = 0;
			elm++;
		}
	}
	invM1 = inverse(M1);
	
	/* Build Matrix M2 Which is equal to [QIU , 0 ; DMU , -I] */

	M2 = newmat(prob->num->rows + inact, prob->num->rows + up, 0);

	/*Build QIU in a mat strcture*/
	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			removerow(QIU, i);
		}
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
		
			removecol(QIU, i);
		}
	}

	/*Build DMU in a mat strcture*/
	for (int i = 1; i < prob->num->cols; i++) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removecol(DMU, i);
		}
	}

    /* first I rows of M2*/

	elm = 0;
	for (int i = 0; i < inact; i++) {

		for (int j = 0; j < up; j++) {
			M2->entries[elm] = QIU->entries[i * up + j];
			elm++;
		}

		for (int j = 0; j < prob->num->rows; j++) {
			M2->entries[elm] = 0;
			elm++;
		}

	}

	/* next M rows of M2*/
	for (int i = 0; i < prob->num->rows; i++) {
		for (int j = 0; j < up; j++) {
			M2->entries[elm] = DMU->entries[i*up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			if (i  == j ) {
				M2->entries[elm] = -1;
				elm++;
			}
			else {
				M2->entries[elm] = 0;
				elm++;
			}			
		}
	}

	/*Calculate 4 components of the W*/
	
	 minvM1 = scalermultiply(invM1 , -1);
   	(*W) = multiply(minvM1, M2);

	//* Mtrix T caculation  *//

	/* QLU */

	for (int i = prob->num->cols ; i >= 1  ; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 2)
		{
			removerow(QLU, i);
		}
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removecol(QLU, i);
		}
	}

	/* QUU */

	for (int i = prob->num->cols; i >=  1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removerow(QUU, i);
			removecol(QUU, i);
		}
	}

	/* QUI */

	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removerow(QUI, i);
		}
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			removecol(QUI, i);
		}
	}

	/* QLI */

	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 2)
		{
			removerow(QLI, i);
		}
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			removecol(QLI, i);
		}
	}

	/* QUI */

	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 1)
		{
			removerow(QLI, i);
		}
		if (cell->partition->part[cnt][i] == 1 || cell->partition->part[cnt][i] == 2)
		{
			removecol(QLI, i);
		}
	}

	/*DML*/

	for (int i = prob->num->cols; i >= 1; i--) {
		if (cell->partition->part[cnt][i] == 0 || cell->partition->part[cnt][i] == 2)
		{
			removecol(DML, i);
		}
	}

	/* Build T=  w1 - w2 * w*/

	w1 = newmat(low+up,prob->num->rows + up ,0);

	/*Build w1*/
	/* first L rows of W1*/

	elm = 0;
	for (int i = 0; i < low; i++) {
		for (int j = 0; j < up; j++) {
			w1->entries[elm] = -QLU->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			w1->entries[elm] = 0;
			elm++;
		}
	}

	/* next U rows of W1*/

	for (int i = 0; i < up; i++) {
		for (int j = 0; j < up; j++) {
			w1->entries[elm] = -QUU->entries[i*up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
	
				w1->entries[elm] = 0;
				elm++;
			
		}
	}

	/*Build w2*/

	w2 = newmat(low + up, prob->num->rows + inact, 0);
	elm = 0;
	DMLT = transpose(DML);
	DMUT  = transpose(DMU);

	/* first L rows of W2*/
	for (int i = 0; i < low; i++) {
		for (int j = 0; j < inact; j++) {
			w2->entries[elm] = -QLU->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			w2->entries[elm] = -DMLT->entries[i * prob->num->rows + j];
			elm++;
		}
	}

	/* next U rows of W2*/

	for (int i = 0; i < up; i++) {
		for (int j = 0; j < inact; j++) {
			w2->entries[elm] = -QUI->entries[i * up + j];
			elm++;
		}
		for (int j = 0; j < prob->num->rows; j++) {
			w2->entries[elm] = -DMUT->entries[i * prob->num->rows + j];
			elm++;
		}
	}

	/*Build T*/

	W2W = multiply(w2 , (*W));
	(*T) = sum(w1, W2W);

	freemat(QII);
	freemat(QIU);
	freemat(QLU);
	freemat(QUU);
	freemat(QUI);
	freemat(QLI);
	freemat(DMU);
	freemat(DML);
	freemat(DMI);
	freemat(M1);
	freemat(M2);
	freemat(invM1);
	freemat(DMLT);
	freemat(DMUT);
	freemat(w2);
	freemat(w1);
	freemat(W2W);
	freemat(minvM1);
	freemat(DMIT);
}



int StocUpdatePart(cellType* cell, probType* prob, sparseVector* bOmega, sparseMatrix* COmega , sparseVector* lOmega, sparseVector* uOmega, solnType* soln , int* basis , int* partIndx) {
	int up = 0, inact = 0, low = 0; /*Number of variables on their bounds*/
	bool newPartFlag = false;
	Mat* W ;
	Mat* T ;
	/* 4b. Calculate the partition */
	 (*partIndx) = AddtoPart(prob, cell, uOmega, lOmega, soln, &newPartFlag, &up, &inact, &low, basis);

	/* 4d. Store the fixed parts of current partition if needed*/

	if (newPartFlag) {

		/* 4d.1 Extract the WT matrices*/

		CalcWT(cell, prob, prob->sp->objQ, prob->Dbar, &W, &T, low, up, inact);

		/* 4d.2 Add the obtained solution to the lambda structure*/

		/* 4.Build [W;T] */

		Mat* WT = CombineWT(prob, W, T, low, up, inact);

		addtoLambdaP(cell, soln, WT, prob, bOmega, uOmega, lOmega, low, up, inact);

		/* 4d.2 Add to alpha and beta Bar*/

		AddtoSigmaP(cell, soln, prob);

		/* 4d.3 add to  delta sol and complete a row*/

		for (int obs = 0; obs < cell->omega->cnt; obs++) {
			/* Add a new row to the delta structure for all observations and the latest lambda (lambdaIdx) */
			bOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[0] - 1;
			COmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[1] - 1;

			uOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[3] - 1;
			lOmega->val = cell->omega->vals[obs] + prob->coord->rvOffset[4] - 1;

			addtoDeltaP(cell, soln, W, T, WT, prob, COmega, bOmega, uOmega, lOmega, obs, (*partIndx), inact, up, low);
		}
		freemat(WT);
		freemat(W);
		freemat(T);
	}
	return partIndx;
}; //EndStocUpdatePart