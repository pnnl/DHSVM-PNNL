/*
 * SUMMARY:      average shadow images from fine to coarse times
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
  FILE  *infile,*outfile;
  char   infilename[255],outfilename[255];
  int nRows;                    /* Number of rows */
  int nCols;                    /* Number of columns */
  int nIn,nOut;
  unsigned char *input;
  unsigned char  *output;
  float *temp;               /* to trap overflows of unsigned char *output */
  int i,j,k;
  int compress;

  if(argc<7) {
    printf("usage is: average_shadow:  \n");
    printf("inputfile, outputfile, # in, # out, nrows, ncols \n");
    printf("the last 4 variable should all be entered as integers \n");
    exit(-1);
  }

  strcpy(infilename, argv[1]);    /* name of the binary uchar input file - no header */
  strcpy(outfilename, argv[2]);   /* name of the binary uchar output file        */
  nIn = GetNumber(argv[3]);       /* number of input images (typically 24 for hourly) */
  nOut = GetNumber(argv[4]);      /* number of output images (typically 8 for 3 hourly) */
  nRows = GetNumber(argv[5]);     /* number of rows in a single image */  
  nCols = GetNumber(argv[6]);     /* number of columns in a single image */
 
  if(fmod((float)nIn,(float)nOut)!=0.0) {
    printf("Number of input images not wholly divisible by number of output images \n");
    exit(-1);
  }

  compress=nIn/nOut;
  
  input = calloc(nRows*nCols, sizeof(unsigned char));
  if (input == NULL)
    exit(-1);

  output = calloc(nRows*nCols, sizeof(unsigned char));
  if (output == NULL)
    exit(-1);

  temp = calloc(nRows*nCols, sizeof(float));
  if (temp == NULL)
    exit(-1);

  if (!(infile = fopen(infilename, "rb"))){
    printf("input file not found \n");
    exit(-1);
  }

  if (!(outfile = fopen(outfilename, "wb"))){
    printf("output file not opened \n");
    exit(-1);
  }

  for(i=0;i<nOut;i++){
    for(k=0;k<nRows*nCols;k++) 
      temp[k]=0.0;
      for(j=0;j<compress;j++){
      fread(input, sizeof(unsigned char), nCols*nRows, infile);
      for(k=0;k<nCols*nRows;k++)
	temp[k]+=(float)input[k]/(float)compress;
    }
    for(k=0;k<nCols*nRows;k++){
	output[k]=(unsigned char)temp[k];   
	if(temp[k]>255.0) output[k]=255;
    }
  fwrite(output,sizeof(unsigned char),nCols*nRows,outfile);   
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


