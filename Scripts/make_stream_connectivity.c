/******************************************************************************************************
Author: Ning Sun September 2013
Program: computes radiation attenuation for each stream segment at each time step

  - connectivity.txt : segment.id , array of upstream segments converging to this segment
	-> headwater segments have no upstream segment
	-> most segments have one upstream segment
	-> confluence segments have an array of upstream segment
  
  - azimuth.txt : ( all in same file)
	-> estimated azimuth, use ArcInfo if need exact azimuth
	-> segments are straight lines, get the 2 extremities x and y and
	    derive the angle from north. The angle is +/- 180 because does 
	    not take into account the flow direction.
            then get the 180-atan(dx/dy)
	->NEVERMIND, azimuth provided in the map file, get the idea anyway 
        -> see SlopeAscpect.c
        -> NOTE that azimuth was derived by Lan but does not seem the default ( default is aspect?)

Modified: Mar 12, 2014
******************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// to be updated for large basins as needed
#define MAXCONV 15000	/* maxinum number of segments converging */ 
#define MAX_XY  15000	/* max number of cell[x,y] crossed by one segment */

int main (int argc, char** argv)
{
 int x, y, id, NX, NY, Nseg, seg;
 float length;
 int i;             /* counter */
 int *next; 		/* downstream seg id for present id*/
 int **before;		/* upstream seg id for present id */
 int *ct_before; 	/* number of upstream seg for present id */
 int *ct_id;		/* number of cells[x,y] crossed by seg id */
 int *segid;		/* seg id for seg number */
 int **seg_id;		/* seg id present in [x,y] */
 float *seg_length; /* segment length */
 int *x_id;		    /* x of [x,y] number ct_id in segment id; */
 int *y_id;		    /* y of [x,y] number ct_id in segment id; */
 int data;          /* temporarily store the destination cell */
 float *length_id;	        /* lenght of the entire cell id */
 float **length_xy;	        /* lenght if segment in gridcell*/ 
 float **azim_xy, dataf;    /* azimuth of the segment */
 float **elev_xy, hxy;      /* elevation of the segment */
 float lid,lxy, hid;
 char convergence[150], outputdir[100];
 int junki;
 float junkf;
 char junkc[350];
 int nskip;        /* header line number in stream map file */
 FILE *fmap,*fnetwork,*fconv;
 int NsegDummy;

if (argc != 6 ){
  printf("Command line arguments: enter <mapfile> <networkfile> <output directory> <Nseg> <skip>\n");
  printf("skip = lines of the header in stream map file\n");
  exit (0);
 }
 sscanf(argv[4],"%d", &Nseg);
 sscanf(argv[5],"%d", &nskip);

 Nseg++; // Nseg is 1 greater than the original hereafter
 NsegDummy = Nseg + 400;  //Increase the dimension of the matrix in case that streamnetwork is mannually 
                          //modified to remove the redundant outlet "-1".
 
 /*handle file names */
 fmap = fopen(argv[1], "r");
 if (fmap == NULL) {
   fprintf(stderr,"NULL  %s \n", argv[1]); 
   exit(-1);
 }
 fprintf(stdout, " %s opened for reading\n", argv[1]);

 fnetwork = fopen(argv[2], "r");
 if (fmap == NULL) {
   fprintf(stderr,"NULL  %s \n", argv[2]); 
   exit(-1);
 }
 fprintf(stdout, " %s opened for reading\n", argv[2]);

 sprintf(outputdir,argv[3]);
 strcpy(convergence, outputdir);
 strcat(convergence, "convergence.txt");
 fprintf(stdout, " %s opened for writing \n", convergence);

 /* allocate */
 if (!(segid = (int*) calloc(NsegDummy, sizeof(int)))) {
   fprintf(stderr,"Failed to allocate variable 'segid'\n"); 
   exit(-1);
 }

 if (!(seg_length = (float*) calloc(NsegDummy, sizeof(float)))) {
   fprintf(stderr,"Failed to allocate variable 'segid'\n"); 
   exit(-1);
 }

 if (!(next = (int*) calloc(NsegDummy, sizeof(int)))) {
   fprintf(stderr,"Failed to allocate variable 'next'\n"); 
   exit(-1);
 }
 if (!(before = (int**) calloc(MAXCONV,sizeof(int*)))) {
   fprintf(stderr,"Failed to allocate variable 'before'\n"); 
   exit(-1);
 }
 for (seg = 0; seg < MAXCONV; seg++){
   if (!(before[seg] = (int*) calloc(NsegDummy, sizeof(int)))) {
	 fprintf(stderr, "Failed to allocate variable 'before'\n");
	 exit(-1);
   }
 }
 if (!(ct_before = (int*) calloc(NsegDummy,sizeof(int)))) {
   fprintf(stderr, "Failed to allocate variable 'ct_before'\n");
   exit(-1);
 }
 if (!(ct_id = (int*) calloc(NsegDummy,sizeof(int)))) {
   fprintf(stderr, "Failed to allocate variable 'ct_id'\n");
   exit(-1);
 }
 if (!(length_id = (float*)calloc(NsegDummy,sizeof(float)))) {
   fprintf(stderr, "Failed to allocate variable 'length_id'\n");
   exit(-1);
 }
 if (!(azim_xy = (float**)calloc(MAX_XY,sizeof(float*)))) {
   fprintf(stderr, "Failed to allocate variable 'azim_xy'\n");
   exit(-1);
 }
 if (!(elev_xy = (float**)calloc(MAX_XY,sizeof(float*)))) {
   fprintf(stderr, "Failed to allocate variable 'elev_xy'\n");
   exit(-1);
 }
 if (!(length_xy = (float**)calloc(MAX_XY,sizeof(float*)))) {
   fprintf(stderr, "Failed to allocate variable 'length_xy'\n");
   exit(-1);
 }
 for (x = 0; x < MAX_XY; x++){
   if (!(length_xy[x] = (float*) calloc(NsegDummy,sizeof(float)))) {
	 fprintf(stderr, "Failed to allocate variable 'length_xy'\n");
     exit(-1);
   }
   if (!(elev_xy[x] = (float*) calloc(NsegDummy,sizeof(float)))) {
	 fprintf(stderr, "Failed to allocate variable 'elev_xy'\n");
     exit(-1);
   }
   if (!(azim_xy[x] = (float*) calloc(NsegDummy,sizeof(float)))) {
	 fprintf(stderr, "Failed to allocate variable 'azim_xy'\n");
     exit(-1);
   }
  }

 /***** read network file *****/
 seg = -1;
 while (!feof(fnetwork)){
   seg++;
   if (seg >= Nseg) {
     fprintf(stderr,
		 "Increase the maximum number of segment in the basin in the code (Nseg)\n");
     exit(-1);
   }  

   /* the 6th column of the network file is the destination segment of 
   the present segment */
   if (fscanf(fnetwork, "%d %d %f %f %d %d \n",&id,&junki,&junkf,&length,&junki,&data)!= 6){
	 /* if the segment is the outlet, then it has 7 columns.
	 Here just to skip reading the last save flags */
	 if (data == 0 || data == -1 ) {
		fgets(junkc, 100, fnetwork);
		seg--;
	 }
	 else {
       fprintf(stderr, "error reading %s at segment %d\n", argv[2], seg);
       exit(-1);
	 }
   }
   segid[seg] = id;
   seg_length[seg] = length;
   
   /* assign the destination cell to the present channel id */
   next[id] = data;
   if (data == -1 ) {
	 data = 0; // outlet case of -1;
	 printf("next-id = %d\n", next[51]);
   }

   ct_before[data]++;
   if ( ct_before[data] >= MAXCONV ) {
     fprintf(stderr,"Increase max num of stream that can converge at once %d\n",ct_before[next[id]]);
     fprintf(stderr,"  seg %d , id, %d, next_id %d\n", seg, id, next[id]);
     exit(-1);
   }
   before[ct_before[data]][data] = id;
 }
 fclose(fnetwork);
 fprintf(stdout,"Read %d segments\n", seg + 1);
 if ((seg+1) != (Nseg-1)) {
   fprintf(stderr,"Error in the number of segment expected\n");
   exit(-1);
 }
 printf("next-id = %d\n", next[51]);
 /************ read map file **************/

 /* skip the headers */
 for (i = 0; i < nskip; i++)
	fgets(junkc, 100, fmap);

 while (!feof(fmap)){
   /* initialize the variables */
   id = -99;
   lxy = -99;
   
   /* read in cell location, channel id, length and azimuth in sequence */
   if (fscanf(fmap,"%d %d %d %f %f %f %f \n",&x, &y, &id, &lxy, &hxy, &junkf, &dataf)!= 7){
     fprintf(stderr, "error reading %s at cell[%d][%d]\n", argv[1], x, y);
     exit(-1);
   }
   length_id[id] += lxy;
   ct_id[id]++;

   if (ct_id[id] >= MAX_XY) {
     fprintf(stderr,"Increase MAX_XY: seg_id %d has (%d+1) upstream cells largert than MAX_XY %d\n",id,ct_id[id],MAX_XY);
     exit(-1);
   }
   length_xy[ct_id[id]][id] = lxy;
   azim_xy[ct_id[id]][id] = dataf;
   elev_xy[ct_id[id]][id] = hxy;
 }
 fclose(fmap);

 /******** write convergence file *****/
 fconv = fopen(convergence, "w");
 if (fconv == NULL) {
   fprintf(stderr,"NULL  %s \n", convergence); 
   exit(-1);
 }

 for (seg = 0; seg < (Nseg-1); seg++){
   id = segid[seg];

   /* compute the weighted avg azimuth for the diff id based on azim_xy **/
   y = ct_id[id] ;
   lid = -99;
   hid = -99;
   if (y == 0 && id > 0) {
     fprintf(stderr,"Error zero xy in seg id %d\n",id);
     exit(-1);
   }
   else {
     lid = 0;
	 hid = 0.;
     for (x = 1; x <= y; x++){
       lid += azim_xy[x][id] * length_xy[x][id] / length_id[id];
	   hid += elev_xy[x][id] * length_xy[x][id] / length_id[id];
     }
     if ( lid < 0 || lid > 360 ) {
       fprintf(stderr,"Error azimuth %f\n",lid);
       exit(-1);
     }
	 if (hid < 0) {
       fprintf(stderr,"Error azimuth %f\n",hid);
       exit(-1);
     }
   }
   fprintf(fconv,"%d %d %.3f %.2f %.2f ",id, next[id], seg_length[seg], hid, lid);
   y = ct_before[id] ;
   if (y > 0)  {
     for (x = 1; x <= y; x++) {
       fprintf(fconv," %d ", before[x][id]);
     }
   }
   fprintf(fconv,"\n");
 }
 fclose(fconv);  

 /***** clean up *****/
 free(segid);
 free(next);
 free(ct_before);
 for (x = 0; x < MAXCONV; x++) 
   free(before[x]);
 free(before);

 free(length_id); 
 for (x = 0; x < MAX_XY; x++) 
   free(length_xy[x]);
 free(length_xy);
 for (x = 0; x < MAX_XY; x++) 
   free(azim_xy[x]);
 free(azim_xy);

 free(ct_id);

 printf("completed .........................\n");

 return 0;
}

