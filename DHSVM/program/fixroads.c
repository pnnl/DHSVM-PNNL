#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXSIZE 255
#define NUMROADS 5
#define MAXSTREAMS 10  /* Max. number of stream crossings per road id. */
#define MAXLINES 200 /* Maximum number of road segments. */
//      program fix_the_roads

//     this program will read in the raw dhsvm input road files and fix them

//     inputs to this program are
//     2. dem of area to examine (or entire domain)
//     3. map of stream locations (i.e. which columns and rows have channels)
//     4. map of road segment locations (including id #, length,height,width,aspect,sink)
//        assume that the dem and mask files have beens stipped of headers
main(int argc, char *argv[]) {
  FILE *fi, *fo;
  int col,row,ncol,nrow,ntot,id,maxid;
  float xll, yll, nodata, cellsize;
  int gotroad[MAXLINES], sscol[MAXLINES], ssrow[MAXLINES], icountsinks[MAXLINES]; 
  int gotstream;
  int mecol,merow;
  float minelev;
  int stcol[MAXSTREAMS],strow[MAXSTREAMS],stelev[MAXSTREAMS];
  int minstid,minstelev;
  float **elev;
  int **stream;
  int ***road;
  int **nroad;
  int **sink;
  int **output;
  int igotsink;
  int scts,sctrf;
  char cjunk[MAXSIZE];
  char inroad[50],instream[50], outroad[50], inmask[50], indem[50];
  float length,height,width,aspect;
  int nskips,nskipr,nlsf,nlrf;
  int i, j, ik;
  int totalsinks;
      
  if(argc!=7) {
    fprintf(stderr,"USAGE:  fixroads <dem> <stream map> <road map in> <road map out> <stream header> <road header>\n",argv[0]);
    exit(0);
  }

  strcpy(indem,argv[1]);
  strcpy(instream,argv[2]);
  strcpy(inroad,argv[3]);
  strcpy(outroad,argv[4]);
  nskips = atoi(argv[5]);
  nskipr = atoi(argv[6]);

 
 
  scts=0;
  sctrf=0;
  igotsink=0;
      
  fprintf(stderr, "reading in elev data and inverting nodata value\n");

  /** OPEN INPUT FILE **/
  if((fi=fopen(indem,"r")) == NULL) {
    printf("Error Opening Input File, %s.\n", indem);
    exit(0);
  }
  
  fscanf(fi,"%*s %d\n",&ncol);
  fscanf(fi,"%*s %d\n",&nrow);
  fscanf(fi,"%*s %f\n",&xll);
  fscanf(fi,"%*s %f\n",&yll);
  fscanf(fi,"%*s %f\n",&cellsize);
  fscanf(fi,"%*s %f\n",&nodata);

  /*******************************************************************/
  /* Allocate memory: */
  /*******************************************************************/

  if (!(elev = (float **) calloc(nrow, sizeof(float *))))
    {
      printf("Cannot allocate memory for precip.\n");
      exit(0);
    }
  for(i=0; i<nrow; i++){
    if(!(elev[i] = (float *) calloc(ncol, sizeof(float))))
      {
	printf("Cannot allocate memory for precip.\n");
	exit(0);
      }
  }

  if (!(stream = (int **) calloc(nrow, sizeof(int *))))
    {
      printf("Cannot allocate memory for precip.\n");
      exit(0);
    }
  for(i=0; i<nrow; i++){
    if(!(stream[i] = (int *) calloc(ncol, sizeof(int))))
      {
	printf("Cannot allocate memory for precip.\n");
	exit(0);
      }
  }

  if (!(road = (int ***) calloc(nrow, sizeof(int **))))
    {
      printf("Cannot allocate memory for precip.\n");
      exit(0);
    }
  for(i=0; i<nrow; i++){
    if(!(road[i] = (int **) calloc(ncol, sizeof(int *))))
      {
	printf("Cannot allocate memory for precip.\n");
	exit(0);
      }
  }
  for(i=0; i<nrow; i++){
    for(j=0; j<ncol; j++) {
      if(!(road[i][j] = (int *) calloc(NUMROADS, sizeof(int))))
	{
	  printf("Cannot allocate memory for precip.\n");
	  exit(0);
	}
    }
  }

  if (!(nroad = (int **) calloc(nrow, sizeof(int *))))
    {
      printf("Cannot allocate memory for precip.\n");
      exit(0);
    }
  for(i=0; i<nrow; i++){
    if(!(nroad[i] = (int *) calloc(ncol, sizeof(int))))
      {
	printf("Cannot allocate memory for precip.\n");
	exit(0);
      }
  }

  if (!(sink = (int **) calloc(nrow, sizeof(int *))))
    {
      printf("Cannot allocate memory for precip.\n");
      exit(0);
    }
  for(i=0; i<nrow; i++){
    if(!(sink[i] = (int *) calloc(ncol, sizeof(int))))
      {
	printf("Cannot allocate memory for precip.\n");
	exit(0);
      }
  }

   if (!(output = (int **) calloc(nrow, sizeof(int *))))
    {
      printf("Cannot allocate memory for precip.\n");
      exit(0);
    }
  for(i=0; i<nrow; i++){
    if(!(output[i] = (int *) calloc(ncol, sizeof(int))))
      {
	printf("Cannot allocate memory for precip.\n");
	exit(0);
      }
  }



  ntot=ncol*nrow;
	
  for(row=0; row<nrow; row++) {
    for(col=0; col<ncol; col++) {
      fscanf(fi, "%f", &elev[row][col]);
      
      if(elev[row][col] == -9999) 
	elev[row][col]=9999;
    }
  }
  fclose(fi);
  
  scts=0;
  fprintf(stderr, "Initializing stream and road locations and sinks to blank\n");

  for(row=0; row<nrow; row++) {
    for(col=0; col<ncol; col++) {
      stream[row][col]=0;
      sink[row][col]=0;
      for(ik=0; ik<NUMROADS; ik++) 
	road[row][col][ik]=0;
    }
  }
  
  fprintf(stderr, "We assume that the road and stream files have 0,0 as their top left corner \n");

  fprintf(stderr, "\nReading in stream locations...\n");

  if((fi=fopen(instream,"r")) == NULL) {
    printf("Error Opening Input File, %s.\n", instream);
    exit(0);
  }

  for(i=0; i<nskips; i++) 
    fgets(cjunk, MAXSIZE, fi);

  nlsf = 0;
  while(fscanf(fi, "%d %d", &col,&row) != EOF) {
    nlsf +=1;
    stream[row][col]=1;
    fgets(cjunk, MAXSIZE, fi);
  }

  fclose(fi);
  fprintf(stderr,"Number of lines in the input streams files = %d\n",nlsf+nskips);


  fprintf(stderr, "reached end of stream locations file\n");
  fprintf(stderr, "final location was col: %d, row: %d\n", col, row);
  
  fprintf(stderr, "reading in road locations\n");
   
  if((fi=fopen(inroad,"r")) == NULL) {
    printf("Error Opening Input File, %s.\n", inroad);
    exit(0);
  }

  for(i=0; i<nskipr; i++) 
    fgets(cjunk, MAXSIZE, fi);

  maxid=0;
  
  nlrf = 0;
  while(fscanf(fi, "%d %d %d", &col,&row, &id) != EOF) {
    nlrf +=1;
    fgets(cjunk, MAXSIZE, fi);
	
    
    if(nroad[row][col] >= NUMROADS) {
      fprintf(stderr, "The number of roads per grid cell exceeds the defined maximum %d\n", nroad[row][col]);
      fprintf(stderr, "Change #define NUMROADS and recompile.\n");
      exit(0);
    }
    
    road[row][col][nroad[row][col]] = id;
    nroad[row][col]+=1;

    gotroad[id-1]=1;

    if (id > maxid) {
      maxid=id;
      if(maxid >= MAXLINES) {
	fprintf(stderr, "Number of road segments exceeds defined maximum.\n");
	fprintf(stderr, "Edit #define MAXLINES and recompile.\n");
	exit(0);
      }
    }
  }
  fclose(fi);
  fprintf(stderr,"Number of lines in the input roads files = %d\n",nlrf+nskipr);
  
  fprintf(stderr, "reached end of road locations file\n");
  fprintf(stderr, "final location was col %d row %d id %d\n",col,row,id);
  fprintf(stderr, "maximum road id number is %d\n",maxid);

  for(i=0; i<maxid; i++) 
    {
      if (gotroad[i] == 0) 
	fprintf(stderr, "missing segment id #%d\n",i+1);
    }

  fprintf(stderr, "Defining sinks over entire basin\n");
  fprintf(stderr, "If road segment intersects stream channel then sink is at stream channel;\n");
  fprintf(stderr, "otherwise sink is at minimum elevation\n");

  for(i=0; i< maxid; i++) {
    /* Must find sink for each road id. */

    gotstream=0;
    minelev=200000.;
           
    for(row=0; row < nrow; row++) {
      for(col=0; col < ncol; col++) {
	
	if(nroad[row][col] > 0) {
	  for(ik=0; ik<nroad[row][col]; ik++) {

	    /* for each grid cell that includes the current road id, save the elevation
	       if it is minimum, and the stream location, if there is one. */

	    if(road[row][col][ik] == i+1) {
	      if(elev[row][col] < minelev) {
		minelev=elev[row][col];
		mecol=col;
		merow=row;
	      }
	      /*	      if(stream[row][col] == 1) {
		stcol[gotstream]=col;
		strow[gotstream]=row;
		stelev[gotstream]=elev[row][col];
		gotstream+=1;
 		if(gotstream >= MAXSTREAMS) {
		  fprintf(stderr, "The number of streams per segment exceeds the defined maximum.\n");
		  fprintf(stderr, "Edit #define MAXSTREAMS and recompile. \n");
		  exit(0);
		}
		}*/
	    }
	  }
	}
      }
    }
                    

    if(gotstream > 0) {
      
      /* There is at least one stream crossing for this road id. */

      if(gotstream == 1) {

	sscol[i]=stcol[0];
	ssrow[i]=strow[0];
	igotsink=igotsink+1;
      }
      else /* There are multiple stream crossings (shouldn't really happen); find lowest. */
	{
	  minstelev=999999;
	  for(j=0; j<gotstream; j++) {
	    if(stelev[j] < minstelev) {
	      minstelev = stelev[j];
	      minstid = gotstream-1;
	    }
	  }

	  sscol[i]=stcol[minstid];
	  ssrow[i]=strow[minstid];
	  igotsink=igotsink+1;
	}
    
      scts=scts+1;       
      //     fprintf(stderr, "This segment has sink at a stream crossing: col %d,row %d\n",sscol[i], ssrow[i]);

    }
  
    if(gotstream == 0) {
      //     fprintf(stderr, "This segment has sink at min elev %.1f: col %d,row %d\n",minelev,mecol,merow);
      sscol[i]=mecol;
      ssrow[i]=merow;
      igotsink=igotsink+1;

      sctrf=sctrf+1;
    }

    // at this point we have found the sink, remember it for raster output
    sink[ssrow[i]][sscol[i]]=1;   
  }   

  fprintf(stderr, "Got the sinks for %d out of %d road segments.\n", igotsink, maxid);

  fprintf(stderr, "Total number of sinks = %d\n", scts+sctrf);

  fprintf(stderr, "Writing new output file...\n ");

  if((fi=fopen(inroad,"r")) == NULL) {
    printf("Error Opening Input File, %s.\n", inroad);
    exit(0);
  }

  if((fo=fopen(outroad,"w")) == NULL) {
    printf("Error Opening Input File, %s.\n", outroad);
    exit(0);
  }

  for(i=0; i<maxid; i++) 
    icountsinks[i]=0;

  for(i=0; i<nskipr; i++) {
    fgets(cjunk, MAXSIZE, fi);
    fprintf(fo, "%s",cjunk);
  }

  for(i=0; i<nlrf; i++) {
    fscanf(fi, "%d %d %d %f %f %f %f", &col,&row, &id, &length, &height, &width, &aspect);

    if(col == sscol[id-1] && row == ssrow[id-1]) {
      icountsinks[id-1]+=1;
      fprintf(fo, "%4d%4d%4d%12.4f%10.4f%10.4f%10.4f\tSINK\n",col,row,id,length,height,width,aspect);
    }
    else
      fprintf(fo, "%4d%4d%4d%12.4f%10.4f%10.4f%10.4f\n",col,row,id,length,height,width,aspect);
    
  }
   
  totalsinks = 0;
  for(i=0; i<maxid; i++) {
    if(icountsinks[i] == 1) 
      totalsinks +=1;
    if (icountsinks[i] == 0) 
      fprintf(stderr, "no sink located for segment %d\n", i);
    
    if (icountsinks[i] > 1) 
      fprintf(stderr, "more than one sink for segment %d\n", i);
  }   
    
  fclose(fi);
  fclose(fo);

  if((fo=fopen("raster.dat","w")) == NULL) {
    printf("Error Opening Input File, raster.dat.\n");
    exit(0);
  }

  fprintf(fo,"ncols %d\n", ncol);
  fprintf(fo,"nrows %d\n", nrow);
  fprintf(fo,"xllcorner %f\n", xll);
  fprintf(fo,"yllcorner %f\n", yll);
  fprintf(fo,"cellsize %f\n", cellsize);
  fprintf(fo,"NODATA_value 0\n");
    
  for(row=0; row < nrow; row++) {
    for(col=0; col < ncol; col++) {

      output[row][col]=0;
      
      if(stream[row][col] == 1) 
	output[row][col]=2;
      if(nroad[row][col] > 0) 
	output[row][col]=3;
      if(sink[row][col] == 1 && stream[row][col] != 1)
	output[row][col]=4;
      if(sink[row][col] == 1 && stream[row][col] == 1)
	output[row][col]=5;               
      
      fprintf(fo,"%d ", output[row][col]);   
	
	}
    fprintf(fo, "\n");
  }
}



