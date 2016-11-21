/*
 * SUMMARY:      calculate skyview factor from a dem
 * USAGE:        
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    March-2000
 * Last Change: 
 *               
 * DESCRIP-END.cd
 * COMMENTS:
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int GetNumber(char *numberStr);

int main(int argc, char **argv)
{
  FILE  *demfile,*outfile;
  char   demfilename[255],outfilename[255];
  int    nRows;                    /* Number of rows */
  int    nCols;                    /* Number of columns */
  float *temp;
  float **elev;
  float **skyview;
  int    i;
  int    ny,nx;
  float  lx,ly;
  double max_angle,angle;
  int    nLook;
  double theta;
  float  dx;
  float  x,y,sx,sy,dz,dist;
  float  mx,my;
  float  start_elev;
  float  max_elev;
  int    stop_flag;


  if(argc<7) {
    printf("usage is: skyview:  \n");
    printf("demfilename, outfilename, # of look direction, nrows, ncols, cellsize\n");
    printf("the last 4 variable should all be entered as integers \n");
    exit(-1);
  }

  strcpy(demfilename, argv[1]);      /* name of the binary float dem input file - no header */
  strcpy(outfilename, argv[2]);      /* name of the binary float skyview output file        */
  nLook = GetNumber(argv[3]);        /* number of directions to look in 8 or 16 are adequate */
  nRows = GetNumber(argv[4]);        /* number of rows in the input dem */
  nCols = GetNumber(argv[5]);        /* number of columns in the input dem */
  dx    = (float)GetNumber(argv[6]); /* the cellsize of the dem (program assumes that */
                                     /* x and y are the same and that the units of dx */
                                     /* are the same units as in the dem)*/
  
  temp = calloc(nRows*nCols, sizeof(float));
  if (temp == NULL)
    exit(-1);

  if (!(demfile = fopen(demfilename, "rb"))){
    printf("dem file not found \n");
    exit(-1);
  }

  if (!(outfile = fopen(outfilename, "wb"))){
    printf("output file not opened \n");
    exit(-1);
  }

  fread(temp, sizeof(float), nCols*nRows, demfile); 

  if (!((elev) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((elev)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }

  if (!((skyview) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((skyview)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }

  max_elev = 0.0;
  for (ny = 0; ny < nRows; ny++) {
    for (nx = 0; nx < nCols; nx++) {
      elev[ny][nx] = temp[ny*nCols + nx]; 
      if(elev[ny][nx]>max_elev) max_elev = elev[ny][nx];
    }
  }

  ly=(float)(nRows*dx-dx);
  lx=(float)(nCols*dx-dx);

  printf("beginning skyview calculations \n");
   
  for (ny = 0; ny < nRows; ny++) {
    for (nx = 0; nx < nCols; nx++) {
      skyview[ny][nx]=0.0;
      start_elev=elev[ny][nx];

      if(start_elev>0) {
		  for (i = 0; i < nLook; i++) {
		  theta=6.283185/((double)nLook)*(double)i;
		  sx=(float)nx*dx+0.5*dx;
		  sy=(float)ny*dx+0.5*dx;
		  x=sx;
	      y=sy;
	      max_angle=0.0;
		  
		  while(x>dx && x < lx && y>dx && y<ly) {
			  x=x+((float)cos(theta))*dx;
			  y=y+((float)sin(theta))*dx;
			  dz=elev[(int)(y/dx)][(int)(x/dx)]-start_elev;
	          dist=sqrt((x-sx)*(x-sx)+(y-sy)*(y-sy));
	          if(dz>0) {
				  angle=atan((double)(dz/dist));
				  if(angle > max_angle) 
					  max_angle = angle;
			  }
		  }
		  skyview[ny][nx]+=(cos(max_angle)*cos(max_angle));
		  }
		  skyview[ny][nx]=skyview[ny][nx]/(float)nLook;
	  }
	}
  }
  for (ny = 0; ny < nRows; ny++) {
    fwrite(skyview[ny],sizeof(float),nCols,outfile); 
  }
}

/*****************************************************************************
  GetNumber()
*****************************************************************************/
int GetNumber(char *numberStr) 
{
  char *endPtr;
  int number = 0;

  number = (int) strtol(numberStr, &endPtr, 0);
  if (*endPtr != '\0')
    exit(-1);

  return number;
}


