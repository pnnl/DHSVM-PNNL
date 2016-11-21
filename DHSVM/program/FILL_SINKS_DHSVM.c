/*
 * SUMMARY:      FILL_SINKS_DHSVM.c - Fill dem sinks in four directions
 *               and impose a slope on all flat areas
 * USAGE:        Preprocessing for DHSVM
 *
 * AUTHOR:       Laura Bowling
 * ORG:          Purdue University, Department of Agronomy
 * ORIG-DATE:    March 15, 2004
 * DESCRIPTION:  This program is for pre-processing DEM files for use in DHSVM.          
 * It performs the following functions: 
 * 1) Fills sinks in 4 directions (arc/info assumes 8 flow directions).
 * 2) Forces flat areas to have known drainage directions by adding incremental 
 *    elevation adjustments.
 * 
 * Usage: <input DEM> <output DEM> <rows> <columns> <NODATA>
 * Dems should be binary floats, as needed for DHSVM input.
 * DESCRIP-END.
 * FUNCTIONS: equal.c from DHSVM
 * COMMENTS: compile with: gcc FILL_SINKS_DHSVM.c -lm -o FILL_SINKS_DHSVM
 * Probably not the most efficient way to do this, but it seems to work.
 */

/******************************************************************************/
/*				    INCLUDES                                  */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/******************************************************************************/
/*				GLOBAL VARIABLES                              */
/******************************************************************************/
#define NDIR 4
#define SINK_HUGE 1e6

typedef struct {
  float Rank;
  int   x;
  int   y;
} ITEM;

int xneighbor[NDIR] = {0, 1, 0, -1};
int yneighbor[NDIR] = {-1, 0, 1, 0};
int DirIndex[NDIR] = {1, 2, 3, 4};

/******************************************************************************/
/*				    FUNCTION DECLARATIONS                     */
/******************************************************************************/
void find_flowdir(int xmin, int ymin, int xmax, int ymax, int ncols, int nrows, float **Dem, int **Dir, float NODATA, int SetOutlets);
float find_outlet(int ncols, int nrows, float **Dem, int **Dir, float NODATA, int **FlowAcc, int *MaxAccum);
int check_sinks(int ncols, int nrows, float **Dem, int **Dir, 
		int *NumUndefined, float NODATA);
int check_outlets(int ncols, int nrows, float **Dem, int **Dir, float NODATA);
void RadialSearch(int ncols, int nrows, float **Dem, int **Dir, int **NumxL, int **NumxR, int **NumyN, int **NumyS, float NODATA);
void assign_pour_point(int **NumxL, int **NumxR, int **NumyN, int **NumyS, 
		       float **Dem, int **Dir, int nrows, int ncols, int max, 
		       int min, float NODATA);
void AllocateArrays(int ***NumxL, int ***NumxR, int ***NumyS, int ***NumyN,
		    int ***Dir, float ***Dem, int ***FlowAcc, unsigned char ***Mask, int ncols, int nrows);
float average(float **Dem,int i, int j, int ncols, int nrows, float NODATA);
unsigned char fequal(float a, float b);
unsigned char fless(float a, float b);
void FlowAccumulation(float **Dem, int ncols, int nrows, float NODATA, int **Dir, int **FlowAcc);
void quick(ITEM *OrderedCells, int count);
void qs(ITEM *item, int left, int right);

/******************************************************************************/
/*				    MAIN PROGRAM                              */
/******************************************************************************/

int main(int argc, char **argv)
{
  int i, j, n;
  int x, y;
  int xn, yn;
  float **Dem;
  int **Dir;
  int **NumxL, **NumxR, **NumyS, **NumyN;
  int **FlowAcc;
  unsigned char **Mask;
  FILE *fi, *fo;
  int nrows, ncols;
  float xmin, xmax, ymin, ymax;
  float NODATA, cellsize;
  float xll, yll;
  int NumSinks, NumUndefined, NumOutlets;
  char InFile[100], OutFile[100], MaskFile[100];
  int NElements;
  float *Matrix;
  unsigned char *MaskArray;
  float OutletElevation;
  int MaxAccum;
  float min;
  int steepestdirection;
  int SetOutlets = 1;

  if(argc != 7) {
    fprintf(stderr, "%s <input dem> <mask> <output dem> <rows> <columns> <NODATA>\n",
	    argv[0]);
    fprintf(stderr, "Dems should be binary float grids, as used by DHSVM.\n");
    fprintf(stderr, "The mask file should be a binary unsigned char grid, as used by DHSVM.\n");
    exit(0);
  }
  strcpy(InFile, argv[1]);
  strcpy(MaskFile, argv[2]);
  strcpy(OutFile, argv[3]);
  nrows = atoi(argv[4]);
  ncols = atoi(argv[5]);
  NODATA = atof(argv[6]);

  if((fi=fopen(InFile,"rb")) == NULL) {
    fprintf(stderr, "Could not open %s\n", InFile);
    exit(0);
  } 

  AllocateArrays(&NumxL,&NumxR,&NumyS,&NumyN,
		 &Dir, &Dem, &FlowAcc, &Mask, ncols, nrows);

  if (!(Matrix = (float *) calloc(ncols*nrows,
				  sizeof(float)))) {
    fprintf(stderr, "Error allocating matrix.\n");
    exit(0);
  }
  NElements = fread(Matrix, sizeof(float), ncols*nrows, fi);

  if(NElements != nrows*ncols) {
    fprintf(stderr, "Problem reading in %s\n",InFile);
    fprintf(stderr, "NElements = %d\n", NElements);
    exit(0);
  }
  fclose(fi);

   if((fi=fopen(MaskFile,"rb")) == NULL) {
    fprintf(stderr, "Could not open %s\n", MaskFile);
    exit(0);
  } 
  if (!(MaskArray = (unsigned char *) calloc(ncols*nrows,
				  sizeof(unsigned char)))) {
    fprintf(stderr, "Error allocating matrix.\n");
    exit(0);
  }
  NElements = fread(MaskArray, sizeof(unsigned char), ncols*nrows, fi);

  if(NElements != nrows*ncols) {
    fprintf(stderr, "Problem reading in %s\n",MaskFile);
    fprintf(stderr, "NElements = %d\n", NElements);
    exit(0);
  }
  fclose(fi);

  for (y = 0, i = 0; y < nrows; y++){
    for (x = 0; x < ncols; x++, i++) {
      Dem[y][x] = Matrix[i];
      Mask[y][x] = MaskArray[i];

      /* Temporary hack to make program deal with mask files.*/
      if(Mask[y][x] == 0)
	Dem[y][x] = NODATA;
    }
  }
 
  /* Fill in flow direction grid. */
  find_flowdir(0, 0, ncols, nrows, ncols, nrows, Dem, Dir, NODATA, SetOutlets);  SetOutlets = 0;



  NumOutlets = check_outlets(ncols, nrows, Dem, Dir, NODATA);

  /* Begin processing DEM */
  NumSinks = check_sinks(ncols, nrows, Dem, Dir, &NumUndefined, NODATA);
  fprintf(stderr, "NumSinks = %d, NumUndefined=%d\n", NumSinks, NumUndefined);

   /* First fill sinks with undefined drainage directions. */
  while(NumUndefined > 0 || NumSinks > 0 || NumOutlets > 1) {

    while(NumUndefined > 0 ) {
      RadialSearch(ncols, nrows, Dem, Dir, NumxL, NumxR, NumyN, NumyS, NODATA);
  
      assign_pour_point(NumxL, NumxR, NumyN, NumyS, Dem, Dir, nrows, 
			ncols, -8, -10, NODATA);
    
   
      /* Fill in flow direction grid. */
      find_flowdir(0, 0, ncols, nrows, ncols, nrows, Dem, Dir, NODATA, SetOutlets);
     
    
      NumSinks = check_sinks(ncols, nrows, Dem, Dir, &NumUndefined, NODATA);
      //    fprintf(stderr, "NumSinks = %d, NumUndefined=%d\n", NumSinks, NumUndefined);
    }
    
    if(NumSinks > 0) {
      RadialSearch(ncols, nrows, Dem, Dir, NumxL, NumxR, NumyN, NumyS, NODATA);
  
      assign_pour_point(NumxL, NumxR, NumyN, NumyS, Dem, Dir, nrows, 
			ncols, 0, -8, NODATA);

 
      /* Fill in flow direction grid. */
      find_flowdir(0, 0, ncols, nrows, ncols, nrows, Dem, Dir, NODATA, SetOutlets);

      NumSinks = check_sinks(ncols, nrows, Dem, Dir, &NumUndefined, NODATA);
      fprintf(stderr, "NumSinks = %d, NumUndefined=%d\n", NumSinks, NumUndefined);
    }
    if(NumSinks==0 && NumUndefined == 0 && NumOutlets > 1) {
        /* Identify true basin outlet. */
      FlowAccumulation(Dem, ncols, nrows, NODATA, Dir, FlowAcc);

      OutletElevation = find_outlet(ncols, nrows, Dem, Dir, NODATA, FlowAcc, &MaxAccum);

      NumOutlets = check_outlets(ncols, nrows, Dem, Dir, NODATA);
      fprintf(stderr, "NumOutlets = %d, OutletElevation=%f, FlowAcc=%d\n", NumOutlets, OutletElevation, MaxAccum);
    
      NumSinks = check_sinks(ncols, nrows, Dem, Dir, &NumUndefined, NODATA);
      fprintf(stderr, "End of outlet check: NumSinks = %d, NumUndefined=%d\n", NumSinks, NumUndefined);
    }
  }
  
   NumOutlets = check_outlets(ncols, nrows, Dem, Dir, NODATA);
   fprintf(stderr, "NumOutlets = %d\n", NumOutlets);
  /* Perform Final Check. */
  

  for (y = 0; y < nrows; y++){
     for (x = 0; x < ncols; x++) {
       if(Dem[y][x] != NODATA){
	 
	 min = SINK_HUGE;
	 for (n = 0; n < NDIR; n++) {
	   xn = x + xneighbor[n];
	   yn = y + yneighbor[n];
	   
	   if (xn >=0 && xn <ncols && yn>=0 && yn<ncols) {
	     if (Dem[yn][xn] != NODATA) {
	       if(Dem[yn][xn] < min)
		 { 
		   min = Dem[yn][xn];
		   steepestdirection = n;
		 }
	     }
	   }
	 }
	 if(min < Dem[y][x]) {
	   Dir[y][x] = DirIndex[steepestdirection];
	 }
	 else{
	   	   fprintf(stderr, "Assigning invalid flow direction, elev=%f.\n",Dem[y][x]);
	     fprintf(stderr, "min=%f, ",min);
	     if (y-1 >=0) fprintf(stderr, "%f, ",Dem[y-1][x]);
	    else fprintf(stderr,"Out of basin, ");
	    if (x+1 < ncols) fprintf(stderr, "%f, ",Dem[y][x+1]);
	    else fprintf(stderr,"Out of basin, ");
	    if (y+1 <nrows) fprintf(stderr, "%f, ",Dem[y+1][x]);
	    else fprintf(stderr,"Out of basin, ");
	    if (x-1 >=0) fprintf(stderr, "%f, ",Dem[y][x-1]);
	    else fprintf(stderr,"Out of basin, ");
	     fprintf(stderr, "\n");
	 }
       }
     }
  }


   if((fo=fopen("FlowAcc.bin","wb")) == NULL) {
    fprintf(stderr, "Could not open %s\n", OutFile);
    exit(0);
  } 

   for (y = 0; y < nrows; y++){
     for (x = 0; x < ncols; x++) {
       Matrix[y * ncols + x] = FlowAcc[y][x];
     }


/******************************************************************************/ 
/* Added the following section so that the flow direction grid would be saved */ 
/******************************************************************************/ 

if((fo=fopen("Dir.bin","wb")) == NULL) { 
fprintf(stderr, "Could not open %s\n", OutFile); 
exit(0); 
} 

for (y = 0; y < nrows; y++){ 
for (x = 0; x < ncols; x++) { 
Matrix[y * ncols + x] = Dir[y][x]; 
} 
} 

/******************************************************************************/ 
/* End of section added so that the flow direction grid would be saved */ 
/******************************************************************************/ 
   }


  NElements = 0;
  NElements = fwrite(Matrix,sizeof(float), nrows*ncols,fo);
 
  if(NElements != nrows*ncols) {
    fprintf(stderr, "Problem writing in %s\n",OutFile);
    fprintf(stderr, "NElements = %d\n", NElements);
    exit(0);
  }
  fclose(fo);

  if((fo=fopen(OutFile,"wb")) == NULL) {
    fprintf(stderr, "Could not open %s\n", OutFile);
    exit(0);
  } 

   for (y = 0; y < nrows; y++){
     for (x = 0; x < ncols; x++) {
       Matrix[y * ncols + x] = Dem[y][x];
     }
   }

  NElements = 0;
  NElements = fwrite(Matrix,sizeof(float), nrows*ncols,fo);
 
  if(NElements != nrows*ncols) {
    fprintf(stderr, "Problem writing in %s\n",OutFile);
    fprintf(stderr, "NElements = %d\n", NElements);
    exit(0);
  }
  fclose(fo);
  free(Matrix);

  for(i=0; i<nrows; i++){
    free(NumxL[i]);
    free(NumxR[i]);
    free(NumyN[i]);
    free(NumyS[i]);
  }
  free(NumxL);
  free(NumxR);
  free(NumyN);
  free(NumyS);

} /* End of Main. */

void find_flowdir(int xmin, int ymin, int xmax, int ymax, int ncols, int nrows, float **Dem, int **Dir, float NODATA, int SetOutlets)
{
  int x, y, n;
  int xn, yn;
  float min;
  int Numoutlets;
  int Numbounding;
  int OldDir;

  for(y=ymin; y<ymax; y++) {
    for(x=xmin; x<xmax; x++) {
      OldDir = Dir[y][x];

      if(Dem[y][x] != NODATA) {
      min = SINK_HUGE;
      Dir[y][x] = 0;
     
      /* First find the pour point for current cell */

      for(n=0; n < NDIR; n++) {
	yn = y + yneighbor[n];
	xn = x + xneighbor[n];

	if(yn >= 0 && yn < nrows && xn >= 0 && xn < ncols) {
	  if(Dem[yn][xn] < min && Dem[yn][xn] != NODATA) {
	    min = Dem[yn][xn];
	    if(min < Dem[y][x]) 
	      Dir[y][x] = DirIndex[n];
	    else
	      Dir[y][x] = -1*DirIndex[n];
	  }
	}
      }
      
      /* Check to see if there are multiple pour points. */
      Numoutlets = 0;
      Numbounding =0;
      for(n=0; n < NDIR; n++) {
	yn = y + yneighbor[n];
	xn = x + xneighbor[n];
	
	if(yn >= 0 && yn < nrows && xn >= 0 && xn < ncols ) {
	  if(Dem[yn][xn] != NODATA) {
	    Numbounding +=1;
	    if(fequal(Dem[yn][xn],min)) {
	      Numoutlets++;
	    }
	  }
	}
      }
      if(Numoutlets > 1) {
	/* Same change in z value in multiple directions. */
	    
	if(Dem[y][x] <= min && Numbounding == NDIR) {
	  /* Part of a sink. */
	  Dem[y][x] = min;
	  Dir[y][x] = -9;
	}
	else if(Dem[y][x] <= min && Numbounding < NDIR) {
	  if(SetOutlets) {
	    /* Basin outlet. */
	    Dir[y][x] = -99;
	  }
	  else {
	    if(OldDir == -99)
	      Dir[y][x] = -99;
	    else
	      Dir[y][x] = -9;
	  } 
	}
	else{
	  // Dir[y][x] = 99;
	}
      }
      else {
	if(Dem[y][x] <= min && Numbounding == NDIR) {
	  /* Part of a sink. */
	  Dem[y][x] = min;
	}
	else if(Dem[y][x] <= min && Numbounding < NDIR) {
	  if(SetOutlets) {
	    /* Basin outlet. */
	    Dir[y][x] = -99;
	  }
	  else {
	    if(OldDir == -99)
	      Dir[y][x] = -99;
	    else
	      Dir[y][x] = -9;
	  } 
	}
      }
      }
      else {
	Dir[y][x] = NODATA;
      }
      if(Dir[y][x] == 0){
	fprintf(stderr, "Dir = 0\n");
	exit(0);
      }
    }
  }
} /* End of function. */


int check_sinks(int ncols, int nrows, float **Dem, int **Dir, 
		int *NumUndefined, float NODATA)
{
  int x, y, n;
  int xn, yn;
  int count;

  count = 0;
  *NumUndefined = 0;
  for(y=0; y<nrows; y++) {
    for(x=0; x<ncols; x++) {
      
      if(Dem[y][x] != NODATA) {
	if(Dir[y][x] > -99 && Dir[y][x] < 99) {
	  
	  if(Dir[y][x] < 0) {
	    count++;
	  }
	  if(Dir[y][x] == -9) {
	    *NumUndefined+=1;
	  }
	}
      }
    }
  }
  //  fprintf(stderr, "\n");
  return count;
} /* End of function. */

int check_outlets(int ncols, int nrows, float **Dem, int **Dir, float NODATA)
{
  int x, y, n;
  int xn, yn;
  int count;

  count = 0;
  for(y=0; y<nrows; y++) {
    for(x=0; x<ncols; x++) {
      
      if(Dem[y][x] != NODATA) {
	if(Dir[y][x] == -99) {
	    count++;
	}
      }
    }
  }
  return count;
} /* End of function. */

void RadialSearch(int ncols, int nrows, float **Dem, int **Dir, 
		  int **NumxL, int **NumxR, int **NumyN, int **NumyS, float NODATA)
{
  int x, y;
  int i, j;

  /* For each cell, find the number of pixels in each direction with
     the same elevation. */
  for(y=0; y<nrows; y++) {
    for(x=0; x<ncols; x++) {
      NumxR[y][x]=0;
      NumxL[y][x]=0;
      NumyN[y][x]=0;
      NumyS[y][x]=0;
    }
  }

  for(y=0; y<nrows; y++) {
    for(x=0; x<ncols; x++) {

      if(Dem[y][x] != NODATA) {
	if(x!=ncols-1) {
	  j=x;
	  while(fequal(Dem[y][j+1],Dem[y][x]) && Dem[y][j+1] != NODATA) {
	    NumxR[y][x] += 1;
	    j++;
	    if(j==ncols-1)
	      break;
	  }
	}
	if(x!=0) {
	  j=x;
	  while(fequal(Dem[y][j-1],Dem[y][x]) && Dem[y][j-1] != NODATA) {
	    NumxL[y][x] += 1;
	    j--;
	    if(j==0)
	      break;
	  }
	}

	if(y!=nrows-1) {
	  j=y;
	  while(fequal(Dem[j+1][x],Dem[y][x]) && Dem[j+1][x] != NODATA) {
	   
	    NumyS[y][x] += 1;
	    j++;
	    if(j==nrows-1)
	      break;
	  }
	}
	if(y!=0) {
	  j=y;
	  while(fequal(Dem[j-1][x], Dem[y][x]) && Dem[j-1][x] != NODATA) {
	    NumyN[y][x] += 1;
	    j--;
	    if(j==0)
	      break;
	  }
	}
    }
    }
  }
}

void assign_pour_point(int **NumxL, int **NumxR, int **NumyN, int **NumyS, 
		  float **Dem, int **Dir, int nrows, int ncols, int max, 
		  int min, float NODATA)
{
  int x, y, n;
  int xn, yn;
  int UNDEFINED;
  int i, j;
  float pour_pt, rim;
  float ave;


  UNDEFINED = 0;
  for(y=0; y<nrows; y++) {
    for(x=0; x<ncols; x++) {

      if(Dir[y][x] < max && Dir[y][x] > min) 
	UNDEFINED = 1;
      
      if(UNDEFINED && Dem[y][x] != NODATA) {

	/* The pour_pt is the lowest boundary cell, it may be less than,
	   greater, or equal to the center cell. The rim is the lowest boundary cell 
	   greater than the center cell.      
	*/
	pour_pt = Dem[y][x];
	rim = SINK_HUGE;
	
	  for(j=x-NumxL[y][x]; j<= x+NumxR[y][x]; j++ ) {
	    for(i=y-NumyN[y][j]; i<=y+NumyS[y][j]; i++ ) {
	    
	      if(j!=0  && Dem[i][j-1] != NODATA) {
		/* Check left boundary. */
		if(fless(Dem[i][j-1],pour_pt)) {
		  pour_pt = Dem[i][j-1];
		}
		if(fless(Dem[i][j-1],rim) && fless(Dem[y][x],Dem[i][j-1])) {
		  rim = Dem[i][j-1];
		}
	      }
	    
	      if(j!=ncols-1 && Dem[i][j+1] != NODATA) {
		/* Check right boundary. */
		if(fless(Dem[i][j+1],pour_pt)) {
		  pour_pt = Dem[i][j+1];
		}
		if(fless(Dem[i][j+1],rim) && fless(Dem[y][x],Dem[i][j+1])) {
		  rim = Dem[i][j+1];
		}
	      }
	    
	      if(i==y-NumyN[y][j] && i !=0 && Dem[i-1][j] != NODATA) {
		/* Check above. */
		if(fless(Dem[i-1][j],pour_pt)) {
		  pour_pt = Dem[i-1][j];
		}
		if(fless(Dem[i-1][j],rim) && fless(Dem[y][x],Dem[i-1][j])) {
		  rim = Dem[i-1][j];
		}
	      }

	      if(i==y+NumyS[y][j] && i !=nrows-1 && Dem[i+1][j] != NODATA) {
		/* Check below. */
		if(fless(Dem[i+1][j],pour_pt)) {
		  pour_pt = Dem[i+1][j];
		}
		if(fless(Dem[i+1][j],rim) && fless(Dem[y][x],Dem[i+1][j])) {
		  rim = Dem[i+1][j];
		}
	      }   
	    }
	  }  /* End of loop around flat cells. */
	    

	  for(j=x-NumxL[y][x]; j<= x+NumxR[y][x]; j++ ) {
	    for(i=y-NumyN[y][j]; i<=y+NumyS[y][j]; i++ ) {

	      if(Dir[i][j] != 0 && Dir[i][j] != -99) {
	
		if(Dem[i][j] != NODATA) {
		  if(pour_pt < Dem[i][j]) {
		    /* Flat area, not a true sink, just need to add gradient. */
		    
		    if(Dir[i][j] < 0) {
		      yn = i + yneighbor[-1*Dir[i][j]-1];
		      xn = j + xneighbor[-1*Dir[i][j]-1];
		    }
		    else {
		      yn = i + yneighbor[Dir[i][j]-1];
		      xn = j + xneighbor[Dir[i][j]-1];
		    }

		    /* Rim defined, these are flat hollows. */
		    if(fless(rim,SINK_HUGE)) {
		      
		      /* Cells with defined flow directions towards cells
			 with the same elevation are "headwaters",
			 raise these first. */
		      if(Dir[i][j] != -9 && fequal(Dem[yn][xn],Dem[i][j])) {
			Dem[i][j] = Dem[i][j] + 0.5 * (rim - Dem[i][j]);
			Dir[i][j] = 0;
		      }
		      else if(Dir[i][j] == -9) {
			ave = average(Dem, i, j, ncols, nrows, NODATA);
			if(ave > 0.5 * (rim + Dem[i][j]) || ave <= Dem[i][j])
			  Dem[i][j] = 0.5 * (rim + Dem[i][j]);
			else
			  Dem[i][j] = ave;
			Dir[i][j] = 0;
		      }
		    }
		    else {
		      /* Rim undefined, these are flat peaks. */
		       /* Cells with defined flow directions are "headwaters",
			 lower these first. */
		      if(Dir[i][j] != -9) {
			Dem[i][j] = Dem[i][j] - 0.5 * (Dem[i][j] - pour_pt);
			Dir[i][j] = 0;
		      }
		    }
		  }
		  else if(fequal(pour_pt, Dem[i][j]) && fless(rim,SINK_HUGE)) {
		    /* Area is a sink, need to raise all elevations
		       and rerun. */
		    Dem[i][j] = rim;
		    Dir[i][j] = 0;
		  }
		  else  { 
		    /*pour_pt > Dem[y][x] || pour_pt == Dem[y][x] && rim undefined */
		    fprintf(stderr, "Search radius can't identify if sink or peak.\n");
		  }
		}
	      }
	    }
	
	  }
      }
      UNDEFINED=0;

    }
  }
 
}

void AllocateArrays(int ***NumxL, int ***NumxR, int ***NumyS, int ***NumyN,
		    int ***Dir, float ***Dem, int ***FlowAcc, unsigned char ***Mask,
		    int ncols, int nrows)
{
  int i;


  if (!((*NumxL) = (int **) calloc(nrows, sizeof(int *))))
    {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
    }
  for(i=0; i<nrows; i++){
    if(!((*NumxL)[i] = (int *) calloc(ncols, sizeof(int))))
      {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
      }
    }

  if (!((*NumxR) = (int **) calloc(nrows, sizeof(int *))))
    {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
    }
  for(i=0; i<nrows; i++){
    if(!((*NumxR)[i] = (int *) calloc(ncols, sizeof(int))))
      {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
      }
    }

   if (!((*NumyN) = (int **) calloc(nrows, sizeof(int *))))
    {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
    }
  for(i=0; i<nrows; i++){
    if(!((*NumyN)[i] = (int *) calloc(ncols, sizeof(int))))
      {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
      }
    }

  if (!((*NumyS) = (int **) calloc(nrows, sizeof(int *))))
    {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
    }
  for(i=0; i<nrows; i++){
    if(!((*NumyS)[i] = (int *) calloc(ncols, sizeof(int))))
      {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
      }
    }

    if (!((*Dem) = (float **) calloc(nrows, sizeof(float *))))
    {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
    }
  for(i=0; i<nrows; i++){
    if(!((*Dem)[i] = (float *) calloc(ncols, sizeof(float))))
      {
      printf("Cannot allocate memory for DEM.\n");
      exit(0);
      }
    }

  if (!((*Dir) = (int **) calloc(nrows, sizeof(int *))))
	{
	  printf("Cannot allocate memory for DEM.\n");
	  exit(0);
	}
  for(i=0; i<nrows; i++){
    if(!((*Dir)[i] = (int *) calloc(ncols, sizeof(int))))
      {
	printf("Cannot allocate memory for DEM.\n");
	exit(0);
      }
  }

  if (!((*FlowAcc) = (int **) calloc(nrows, sizeof(int *))))
	{
	  printf("Cannot allocate memory for DEM.\n");
	  exit(0);
	}
  for(i=0; i<nrows; i++){
    if(!((*FlowAcc)[i] = (int *) calloc(ncols, sizeof(int))))
      {
	printf("Cannot allocate memory for DEM.\n");
	exit(0);
      }
  }

    if (!((*Mask) = (unsigned char **) calloc(nrows, sizeof(unsigned char *))))
	{
	  printf("Cannot allocate memory for DEM.\n");
	  exit(0);
	}
  for(i=0; i<nrows; i++){
    if(!((*Mask)[i] = (unsigned char *) calloc(ncols, sizeof(unsigned char))))
      {
	printf("Cannot allocate memory for DEM.\n");
	exit(0);
      }
  }
}  

float average(float **Dem,int i, int j, int ncols, int nrows, float NODATA)
{
  float ave, count;

  ave = 0.;
  count = 0.;

  if(i != 0) {
    if(Dem[i-1][j] != NODATA) {
      ave += Dem[i-1][j];
      count+=1.;
    }
  }

  if(i!=nrows-1) {
    if(Dem[i+1][j] != NODATA) {
       ave += Dem[i+1][j];
      count+=1.;
    }
  }

    if(j != 0) {
    if(Dem[i][j-1] != NODATA) {
      ave += Dem[i][j-1];
      count+=1.;
    }
  }

  if(j!=ncols-1) {
    if(Dem[i][j+1] != NODATA) {
       ave += Dem[i][j+1];
      count+=1.;
    }
  }

  ave = ave/count;

  return ave;
}

/* -------------------------------------------------------------
   find_outlet()
   ------------------------------------------------------------- */
float find_outlet(int ncols, int nrows, float **Dem, int **Dir, float NODATA, int **FlowAcc, int *MaxAccum)
{
  float minimum;
  int maximum;
  int i, j;

  /* Outlet should have minimum elevation and maximum flow accumulation.*/
  minimum = SINK_HUGE;
  maximum = 0;
  for(i=0; i<nrows; i++) {
    for(j=0; j<ncols; j++) {
      if(Dem[i][j] != NODATA && Dem[i][j] < minimum) {
	if(Dir[i][j] == -99)
	  minimum = Dem[i][j];
      }
      if(FlowAcc[i][j] != NODATA && FlowAcc[i][j] > maximum) {
	if(Dir[i][j] == -99)
	  maximum = FlowAcc[i][j];
      }
    }
  }

  if(minimum < SINK_HUGE && maximum > 0) {
    for(i=0; i<nrows; i++) {
      for(j=0; j<ncols; j++) {
	if(Dem[i][j] != NODATA && Dir[i][j] == -99) {
	  if(Dem[i][j] == minimum)
	    {
	      fprintf(stdout, "Dem = %f, i=%d, j=%d\n",Dem[i][j], i, j);
	      fprintf(stdout, "FlowAcc = %d, Dem = %f, i=%d, j=%d\n",FlowAcc[i][j],Dem[i][j], i, j);
	      if(FlowAcc[i][j] == maximum) {
		fprintf(stdout, "FlowAcc = %d, Dem = %f, i=%d, j=%d\n",FlowAcc[i][j],Dem[i][j], i, j);
		Dir[i][j] = -99;
	      }
	      else
		Dir[i][j] = -9;
	    }
	  else
		Dir[i][j] = -9;
	}
      }
    }
  }
  *MaxAccum = maximum;
  return minimum;
}
  
  
  
/* -------------------------------------------------------------
   FlowAccumulation
   ------------------------------------------------------------- */
void FlowAccumulation(float **Dem, int ncols, int nrows, float NODATA, int **Dir, int **FlowAcc)
{
  int x, xn;
  int y, yn;
  int count, k, n;
  ITEM *OrderedCells;
  
  count = 0;
  for (y = 0; y < nrows; y++) {
    for (x = 0; x < ncols; x++) {

      if(Dem[y][x] != NODATA) {
	count +=1;
	FlowAcc[y][x] =1;
      }
      else
	FlowAcc[y][x] = NODATA;
    }
  }

  /* Create a structure to hold elevations of all cells within the mask 
     and the y,x of those cells.*/
  
  if (!(OrderedCells = (ITEM *) calloc(count, sizeof(ITEM))))
    {
      fprintf(stderr, "Error allocating memory in FlowAccumulation().\n");
      exit(0);
    }

  k = 0;
  for (y = 0; y < nrows; y++) {
    for (x = 0; x < ncols; x++) {
    
      /* Save the elevation, y, and x in the ITEM structure. */
      if (Dem[y][x] != NODATA) {
	OrderedCells[k].Rank = Dem[y][x];
	OrderedCells[k].y = y;
	OrderedCells[k].x = x;
	k++;
      }
    }
  }
 
  /* Sort Elev in descending order-- Elev.x and Elev.y hold indices. */
  
  quick(OrderedCells, count);


  /* Use the Dir grid to find flow accumulation. */
  for(k=count-1; k>=0; k--) {
    for(n=0; n < NDIR; n++) {
      y = OrderedCells[k].y;
      x = OrderedCells[k].x;
      if(Dir[y][x] == n+1) {
	yn = OrderedCells[k].y + yneighbor[n];
	xn = OrderedCells[k].x + xneighbor[n];
	FlowAcc[yn][xn] += FlowAcc[y][x];
      }
    }
  }
  free(OrderedCells);
  return;
}

/* -------------------------------------------------------------
   QuickSort
   ------------------------------------------------------------- */

/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/

void quick(ITEM *OrderedCells, int count)
{
  qs(OrderedCells,0,count-1);
}

void qs(ITEM *item, int left, int right)
/**********************************************************************
        this is the quick sort subroutine - it returns the values in
        an array from high to low.
**********************************************************************/
{
  register int i,j;
  ITEM x,y;

  i=left;
  j=right;
  x=item[(left+right)/2];

  do {
    while(item[i].Rank<x.Rank && i<right) i++;
    while(x.Rank<item[j].Rank && j>left) j--;

    if (i<=j) {
      y=item[i];
      item[i]=item[j];
      item[j]=y;
      i++;
      j--;
    }
  } while (i<=j);

  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right);

}

/*
 * SUMMARY:      equal.c - Determine whether two float or two doubles are equal
 *                         to within machine precision
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Wed Feb  3 08:39:13 1999
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     Adapted from netcdf library
 * $Id: equal.c,v 1.4 2003/07/01 21:26:29 olivier Exp $     
 */

/******************************************************************************/
/*				    INCLUDES                                  */
/******************************************************************************/
#ifdef TEST_EQUAL
#include <stdio.h>
#include <stdlib.h>
#endif

#ifndef ABSVAL
#define ABSVAL(x)  ( (x) < 0 ? -(x) : (x) )
#endif

/* if float.h is not available on your system, comment out the following line,
   and uncomment the line after that */
#include <float.h>
/* #define NO_FLOAT_H */
#include <math.h>

/******************************************************************************/
/*				GLOBAL VARIABLES                              */
/******************************************************************************/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef FLT_DIGITS
#define FLT_DIGITS 7		/* default sig. digits for float data */
#endif
#ifndef DBL_DIGITS
#define DBL_DIGITS 15		/* default sig. digits for double data */
#endif
static double double_eps;
static float float_eps;
static int initeps = 0;

/******************************************************************************/
/******************************************************************************/
/*			           FUNCTIONS                                  */
/******************************************************************************/
/******************************************************************************/

static float float_epsilon(void);
static void init_epsilons(void);

/******************************************************************************/
/*				     fequal                                   */
/******************************************************************************/
unsigned char fequal(float a, float b)
{
  if (!initeps) {		/* make sure epsilons get initialized */
    init_epsilons();
    initeps = 1;
  }

  /* Two float values anly need to be equal to within machine precision */
  if ((a > 0) == (b > 0) &&	/* prevents potential overflow */
      (ABSVAL(a - b) <= ABSVAL(float_eps * b)))
    return TRUE;
  else
    return FALSE;
}

/******************************************************************************/
/*				     fless                                    */
/******************************************************************************/
unsigned char fless(float a, float b)
{
  if (!initeps) {		/* make sure epsilons get initialized */
    init_epsilons();
    initeps = 1;
  }

  /* Two float values anly need to be equal to within machine precision */
  if ((a > 0) == (b > 0) &&	/* prevents potential overflow */
      (ABSVAL(a - b) <= ABSVAL(float_eps * b)))
    return FALSE;               /* Values are equal, not less than. */
  else if(a < b)
    return TRUE;
  else
    return FALSE;
}

/******************************************************************************/
/*				 float_epsilon                                */
/******************************************************************************/
static float float_epsilon(void)
{
  float float_eps;
#ifndef NO_FLOAT_H
  float_eps = FLT_EPSILON;
#else /* NO_FLOAT_H */
  {
    float etop, ebot, eps;
    float one = 1.0;
    float two = 2.0;
    etop = 1.0;
    ebot = 0.0;
    eps = ebot + (etop - ebot) / two;
    while (eps != ebot && eps != etop) {
      float epsp1;

      epsp1 = one + eps;
      if (epsp1 > one)
	etop = eps;
      else
	ebot = eps;
      eps = ebot + (etop - ebot) / two;
    }
    float_eps = two * etop;
  }
#endif /* NO_FLOAT_H */
  return float_eps;
}

/******************************************************************************/
/*				 init_epsilons                                */
/******************************************************************************/
static void init_epsilons(void)
{
  float_eps = float_epsilon();
}

