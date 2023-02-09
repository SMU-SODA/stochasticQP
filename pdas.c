#include "../stochasticQP.h"
extern configType config;

Mat* pdas(sparseMatrix* D, int* partition, sparseMatrix* Q, dVector coef , dVector Rhs , int cols , int rows , dVector Ubound , dVector Lbound){

	int flage = 0;
	int up;
	int low;
	int inact;
	sparseMatrix* H = BuildHess(Q);
	Mat* M1;
	Mat* M2;
	Mat* M3;
	Mat* M4;
	Mat* M5;
	Mat* M6;
	Mat* M7;
	Mat* M8;
	Mat* M9;
	Mat* M10;
	Mat* DMLT;
	Mat* DMUT;
	Mat* DMIT;
	Mat* invM1;
	Mat* QII;
	Mat* QIU;
	Mat* QIL;
	Mat* QLU;
	Mat* QUU;
	Mat* QUI;
	Mat* QLI;
	Mat* QLL;
	Mat* QUL;
	Mat* DMU;
	Mat* DML;
	Mat* DMI;
	Mat* temp1;
	Mat* temp2;
	Mat* temp3;
	Mat* rhs1;
	Mat* ypi;
	Mat* temp4;
	Mat* temp5;
	Mat* temp6;
	Mat* temp7;
	Mat* temp8 ;
	Mat* munu ;
	int* VPU;
	int* VPL;
	int* VDU;
	int* VDL;
	int* part = arr_alloc(cols+1 , int);
	for(int i = 1 ; i <= cols ; i++){
		part[i] = partition[i];
	}
	while(flage == 0){
		inact = 0;
		low = 0;
		up = 0;

		for(int i = 1 ; i <= cols; i++){
			if(part[i]==0){
				inact++;
			}
		}
		for(int i = 1 ; i <= cols; i++){
			if(part[i]==2){
				up++;
			}
		}
		for(int i = 1 ; i <= cols; i++){
			if(part[i]==1){
				low++;
			}
		}
		QII = transSparsM(H, cols, cols);
		QIU = transSparsM(H, cols, cols);
		QIL = transSparsM(H, cols, cols);
		QLU = transSparsM(H, cols, cols);
		QUU = transSparsM(H, cols, cols);
		QUI = transSparsM(H, cols, cols);
		QLI = transSparsM(H, cols, cols);
		QLL = transSparsM(H, cols, cols);
		QUL = transSparsM(H, cols, cols);
		DMU = transSparsM(D, cols, rows);
		DML = transSparsM(D, cols, rows);
		DMI = transSparsM(D, cols, rows);

		/*A. [QII , DMIT ; DMI , 0]*/

		/*A1. Build Q(II) */
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 1 || part[i] == 2) {
				removeCol(QII, i);
				removeRow(QII, i);
			}
		}

		/*A2. Build D(MI)*/
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 1 || part[i] == 2) {
				removeCol(DMI, i);
			}
		}
		DMIT = transpose(DMI);

		M1 = newmat(rows + inact, rows + inact, 0);

		/* A3. Build Matrix M1 which is equal to [ QII , DMIT ; DMI , 0 ] */
		/* A3.1. place the first I rows */

		int elm = 0;
		for (int i = 0; i < inact; i++) {
			for (int j = 0; j < inact; j++) {
				M1->entries[elm] = QII->entries[i * inact + j];
				elm++;
			}
			for (int j = 0; j < rows; j++) {
				M1->entries[elm] = DMIT->entries[i * rows + j];
				elm++;
			}
		}
		/* A3.2. place the next M rows */

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < inact; j++) {
				M1->entries[elm] = DMI->entries[i *inact + j];
				elm++;
			}
			for (int j = 1; j <= rows; j++) {
				M1->entries[elm] = 0;
				elm++;
			}
		}

		invM1 = inverse(M1);
		// showmat(M1);
		/* B. -[QIL ; DML] */

		/*B1. DML*/
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 2) {
				removeCol(DML, i);
			}
		}

		/*B2. QIL */
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 1 || part[i] == 2) {
				removeRow(QIL, i);
			}
			if (part[i] == 0 || part[i] == 2) {
				removeCol(QIL, i);
			}
		}

		/*B3. -[QIL;DML]*/
		M2 = newmat( inact + rows , low , 0);
		/* B3.1. place the first I rows */

		elm = 0;
		for (int i = 0; i < inact; i++) {
			for (int j = 0; j < low; j++) {
				M2->entries[elm] = -QIL->entries[i * low + j];
				elm++;
			}
		}
		/* B3.2. place next M rows */
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < low; j++) {
				M2->entries[elm] = -DML->entries[i *low + j];
				elm++;
			}
		}

		/* C. -[QIU;DMU] */
		/*C1. Build DMU in a mat strcture*/
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 1) {
				removeCol(DMU, i);
			}
		}
		/*C2. Build QIU in a mat structure*/
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 1 || part[i] == 2) {
				removeRow(QIU, i);
			}
			if (part[i] == 0 || part[i] == 1) {
				removeCol(QIU, i);
			}
		}
		M3 = newmat( inact + rows , up , 0);
		/* C3.1. place the first I rows */
		elm = 0;
		for (int i = 0; i < inact; i++) {
			for (int j = 0; j < up; j++) {
				M3->entries[elm] = -QIU->entries[i * up + j];
				elm++;
			}
		}
		/* C3.2. place next M rows */
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < up; j++) {
				M3->entries[elm] = -DMU->entries[i *up + j];
				elm++;
			}
		}
		/* D. -[dI; -rho] */
		elm = 0;
		M4 = newmat( inact + rows , 1 , 0);
		for(int i = 1; i <= cols; i++){
			if(part[i]==0){
				M4->entries[elm] = -coef[i];
				elm++;}
		}
		for(int i = 1; i <= rows; i++){
			M4->entries[elm] = Rhs[i];
			elm++;
		}

		/* F. Yundlow , Ybarup*/
		M5 = newmat( low , 1 , 0);
		M6 = newmat( up , 1 , 0);
		elm= 0;
		for(int i = 1; i <= cols; i++){
			if(part[i]==1){
				M5->entries[elm] = Lbound[i];
				elm++;}
		}
		elm= 0;
		for(int i = 1; i <= cols; i++){
			if(part[i]==2){
				M6->entries[elm] = Ubound[i];
				elm++;}
		}
		/* E. Calculate the right-hand side*/

		temp1;
		temp2;
		temp1 = multiply(M2,M5);
		temp2 = multiply(M3,M6);
		temp3 = sum(temp1 , temp2);
		rhs1 = sum(temp3,M4);
		ypi = multiply(invM1,rhs1);


		/* -[QLL;QUL] */

		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 2) {
				removeRow(QLL, i);
			}
			if (part[i] == 0 || part[i] == 2) {
				removeCol(QLL, i);
			}
		}
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 1) {
				removeRow(QUL, i);
			}
			if (part[i] == 0 || part[i] == 2) {
				removeCol(QUL, i);
			}
		}
		M7 = newmat( low + up , low , 0);
		/* G.1. place the first low rows */
		elm = 0;
		for (int i = 0; i < low; i++) {
			for (int j = 0; j < low; j++) {
				M7->entries[elm] = QLL->entries[i * low + j];
				elm++;
			}
		}
		/* G.2. place next M rows */
		for (int i = 0; i < up; i++) {
			for (int j = 0; j < low; j++) {
				M7->entries[elm] = -QUL->entries[i *low + j];
				elm++;
			}
		}
		/* -[QLU;QUU] */
		/* QLU */
		for (int i = cols; i >= 1; i--) {
			if ( part[i] == 0 || part[i] == 2)
			{
				removeRow(QLU, i);
			}
			if (part[i] == 0 || part[i] == 1)
			{
				removeCol(QLU, i);
			}
		}

		/* QUU */
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 1)
			{
				removeRow(QUU, i);
				removeCol(QUU, i);
			}
		}

		M8 = newmat( low + up , up , 0);
		/* G.1. place the first low rows */
		elm = 0;
		for (int i = 0; i < low; i++) {
			for (int j = 0; j < up; j++) {
				M8->entries[elm] = QLU->entries[i * up + j];
				elm++;
			}
		}
		/* G.2. place next M rows */
		for (int i = 0; i < up; i++) {
			for (int j = 0; j < up; j++) {
				M8->entries[elm] = -QUU->entries[i *up + j];
				elm++;
			}
		}
		/* -[QLI,DMLT;QUI,DMUT] */
		/* QLI */
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 2)
			{
				removeRow(QLI, i);
			}
			if (part[i] == 1 || part[i] == 2)
			{
				removeCol(QLI, i);
			}
		}
		/* QUI */
		for (int i = cols; i >= 1; i--) {
			if (part[i] == 0 || part[i] == 1)
			{
				removeRow(QUI, i);
			}
			if (part[i] == 1 || part[i] == 2)
			{
				removeCol(QUI, i);
			}
		}
		DMLT = transpose(DML);
		M9 = newmat( low  + up, rows + inact, 0);

		elm = 0;
		for (int i = 0; i < low; i++) {
			for (int j = 0; j < inact; j++) {
				M9->entries[elm] = QLI->entries[i * inact + j];
				elm++;
			}
			for (int j = 0; j < rows; j++) {
				M9->entries[elm] = DMLT->entries[i * rows + j];
				elm++;
			}
		}
		/* place the next Up rows */
		DMUT = transpose(DMU);
		for (int i = 0; i < up; i++) {
			for (int j = 0; j < inact; j++) {
				M9->entries[elm] = -QUI->entries[i *inact + j];
				elm++;
			}
			for (int j = 0; j < rows; j++) {
				M9->entries[elm] = -DMUT->entries[i * rows + j];
				elm++;
			}
		}
		/* D. -[dL; dU] */
		elm = 0;
		M10 = newmat( up + low , 1 , 0);
		for(int i = 1; i <= cols; i++){
			if(part[i]==1){
				M10->entries[elm] = coef[i];
				elm++;}
		}

		for(int i = 1; i <= cols; i++){
			if(part[i]==2){
				M10->entries[elm] = -coef[i];
				elm++;}
		}

		/* F. */
		/* E. Calculate the right-hand side*/


		temp4 = multiply(M9,ypi);
		temp5 = multiply(M7,M5);
		temp6 = multiply(M8,M6);
		temp7 = sum(temp4 , temp5);
		temp8 = sum(temp6 , temp7);
		munu = sum(temp8,M10);



		/* obtain set VD and VP*/
		int VPcntU = 0;
		int VPcntL = 0;
		VPU = (int*)arr_alloc(cols + 1,int);
		VPL = (int*)arr_alloc(cols + 1,int);
		elm = 0;
		for(int i = 1; i <= cols ; i++){
			if(part[i] == 0){
				if( ypi->entries[elm] > Ubound[i]+config.TOLERANCE )
				{VPU[VPcntU+1] = i;
				VPcntU++;}
				if( ypi->entries[elm] < Lbound[i]-config.TOLERANCE )
				{VPL[VPcntL+1] = i;
				VPcntL++;}
				elm++;
			}
		}

		int VDcntU = 0;
		VDU = (int*)arr_alloc(cols + 1 , int);
		int VDcntL = 0;
		VDL = (int*)arr_alloc(cols + 1 , int);
		elm = 0;

		for(int i = 1; i <= cols ; i++){
			if(part[i] == 1){
				if((munu->entries[elm] < -config.TOLERANCE))
				{VDL[VDcntL+1] = i;
				VDcntL++;}
				elm++;
			}}

		for(int i = 1; i <= cols ; i++){
			if(part[i] == 2){
				if((munu->entries[elm] < -config.TOLERANCE))
				{VDU[VDcntU+1] = i;
				VDcntU++;}
				elm++;
			}}


		if(VDcntL + VDcntU + VPcntL + VPcntU != 0){

			//    	   int numExchg = minimum(VDcntL + VDcntU, VPcntL + VPcntU);
			//
			//    	   for ( int i = 0; i < numExcg; i++ ) {
			//
			//    	   }

			/* Create a new partition*/
			for(int i = 1; i <= VPcntL ; i++){
				part[VPL[i]] = 1;
			}
			for(int i = 1; i <= VPcntU ; i++){
				part[VPU[i]] = 2;
			}
			for(int i = 1; i <= VDcntL ; i++){
				part[VDL[i]] = 0;
			}
			for(int i = 1; i <= VDcntU ; i++){
				part[VDU[i]] = 0;
			}

			freemat(M1);freemat(M2);freemat(M3);freemat(M4);freemat(M5);freemat(M6);freemat(M7);freemat(M8);
			freemat(M9);freemat(M10);
			freemat(DMLT);
			freemat(DMUT);
			freemat(DMIT);
			freemat(invM1);
			freemat(QII);
			freemat(QIU);
			freemat(QIL);
			freemat(QLU);
			freemat(QUU);
			freemat(QUI);
			freemat(QLI);
			freemat(QLL);
			freemat(QUL);
			freemat(DMU);
			freemat(DML);
			freemat(DMI);
			freemat(temp1); freemat(temp2); freemat(temp3); freemat(rhs1);freemat(ypi);
			freemat(temp4);freemat(temp5);freemat(temp6);freemat(temp7);
			freemat(temp8);freemat(munu);
			mem_free(VPL);           mem_free(VPU);
			mem_free(VDL);           mem_free(VDU);
		}
		else {flage = 1;}
	}
	freemat(M1);freemat(M2);freemat(M3);freemat(M4);freemat(M5);freemat(M6);freemat(M7);freemat(M8);
	freemat(M9);freemat(M10);
	freemat(DMLT);
	freemat(DMUT);
	freemat(DMIT);
	freemat(invM1);
	freemat(QII);
	freemat(QIU);
	freemat(QIL);
	freemat(QLU);
	freemat(QUU);
	freemat(QUI);
	freemat(QLI);
	freemat(QLL);
	freemat(QUL);
	freemat(DMU);
	freemat(DML);
	freemat(DMI);
	freemat(temp1); freemat(temp2); freemat(temp3); freemat(rhs1);//freemat(ypi);
	freemat(temp4);freemat(temp5);freemat(temp6);freemat(temp7);
	freemat(temp8);freemat(munu);
	mem_free(VPL);           mem_free(VPU);
	mem_free(VDL);           mem_free(VDU);
	freeSparseMatrix(H);
	mem_free(part);
	return ypi;

}
