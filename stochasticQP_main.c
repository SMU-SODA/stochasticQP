#include "stochasticQP.h"

cString outputDir;
long MEM_USED;
configType config;

int main(int argc, char* argv[]) {
	cString inputDir, probname;
	int numStages;
	stocType* stoch;
	probType** prob;
	cellType* cell = NULL;
	char configFile[BLOCKSIZE];

	/* Obtain parameter input from the command line */
	parseCmdLine(argc, argv, &probname, &inputDir);

	/* read algorithm configuration file */
#if _WIN64
	strcpy(configFile, "C:\\Users\\Niloofar\\source\\repos\\stochasticQP\\config.sd");
#else
	strcpy(configFile, "./config.sd");
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
