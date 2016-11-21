import java.io.*;
class AddAat2 {
// Copyleft Harvey Greenberg, UW hgreen@u.washington.edu  Oct 12, 1998
// Calls class Aat.class

// This is based on AddAat.java, but the methodology has been changed,
// and it also calculates stream order, downstream arc, and principal
// upstream arc.
// AddAat assumes that you have a river network represented as one or more
// topological trees, with all segments pointing downstream.  The strahler
// calculation may be inaccurate if more than two major downstream rivers
// meet at a junction.  There is a limit of 7 downstream arcs at a node.

// The input/output file is an AAT (arc attribute table) associated with
// an ARC/INFO "coverage"("vector" data set).  There are a number of
//requirements for this coverage:
//	All arcs point downstream.
//	There is always one  (or zero) downstream arc
//	Nodes are correctly numbered.  (Use "renode" to be safe.)
//	For "uparc" and "downarc" to be meaningful, arcs must have
//		uniques IDs (Try cover-id = cover#.)
//	There should be an item "local", which contains the local
//		contributing area in cells.
//	There should be an item "maxmsq", which will receive the contributing
//		area, in square meters, as measured at the downstream end.
//	There should be an item "meanmsq", which will receive the mean of the
//		upstream area and downstream area.
//      There should be an item "shreve", which will receive the stream order.
//		This is the count of upstream sources.
//      There should be an item "strahler", which will receive the stream order,
//		as calculated by Strahler's method.
//      There should be an item "segorder", which will receive the stream order,
//		as calculated by DHSVM's method.
//      There should be an item "downarc", which will receive the id (cover-id)
//		of the downstream arc.
//      There should be an item "uparc", which will receive the ID of the arc
//		upstream arc with the largest area.
//	All the above items should be numeric.  Binary works well.
//		this class does not check types.
// The area of a cell (used to convert from cell counts to contributing area)
// is hardwired below.

// Some of the error handling is crude.  It will simply exit if it does not
// find the items it wants.
// The AAt class completely rewrites the AAT file.  It has been tested for
// safety, but it is not guaranteed to not corrupt the AAT file.
// Remember, java defaults to a memory limit of 16 megs.  Use the mx switch,
// e.g. "java AddAat -mx200m mycover" for large data. 
// version 1.01: added segorder


  public static void main (String args[])  throws IOException {
    Aat myatt;
    int maxid,ilocal,imax,imean,ishreve,shreve=1,istrahler,isegorder,idown,
	iup,ndone=99,maxdone=0,it=0,maxmax=0,nup,ii;
    int area=22500, halfarea = 11250;  // CRUDE WAY TO SET CELL SIZE
    String covername;
    boolean verbose = false;

    System.out.println("AddAat2, version 1.01 Wed Oct 12 15:36:46 PDT 1998, using a cellsize of " + area + " square meters.");

    if (args.length < 1){
	System.out.println("The cover name is a required argument.");
        return;
    }
    else
	covername = args[0];
//    System.out.println("Cover: " + covername );
    try{
	myatt  = new Aat(covername);
    }
    catch (AatError e){
	System.out.print("Error loading AAT file:");
	System.out.println(e.getMessage());
	return;
    }
//    myatt.report();

    maxid = myatt.nrecords;
	ilocal = myatt.findItem("local");	// input (contib. area)
	imax = myatt.findItem("maxmsq");	// output (tot. contib at outlet)
	imean = myatt.findItem("meanmsq");	// output (maxmsq - half of local)
	ishreve = myatt.findItem("shreve");	//output (NEW)
	istrahler = myatt.findItem("strahler");	//output (NEW)
	isegorder = myatt.findItem("segorder");	//output (NEWER)
	idown = myatt.findItem("downarc");	//output (NEW)
	iup = myatt.findItem("uparc");		//output (NEW)
	if(ilocal < 0) return;
	if(imax < 0) return;
	if(imean < 0) return;
	if(ishreve < 0) return;
	if(istrahler < 0) return;
	if(isegorder < 0) return;
	if(idown < 0) return;
	if(iup < 0) return;
    int flags[] = new int[maxid];
    int down[] = new int[maxid];
    int upi[] = new int[8];
// first pass: id sources and record downarc
      for(int i = 0;i < maxid;i++){
	  myatt.ivalues[imax][i] = myatt.ivalues[ilocal][i];
	  flags[i] = 1;  // tentative
	  myatt.ivalues[idown][i] = -1;
	  down[i] = -1;
	  for(int j = 0;j< maxid;j++){
	    if(myatt.fnode[i] == myatt.tnode[j]){  // found upstream arc
	      flags[i] = 0;
	    }
	    if(myatt.fnode[j] == myatt.tnode[i]){  // found downstream arc
	      myatt.ivalues[idown][i] = myatt.id[j];
	      down[i] = j; // use handy array within this program
	    }
	  } // search for connected arcs
	  if(flags[i] == 1){
	    myatt.ivalues[ishreve][i] = 1;
	    myatt.ivalues[istrahler][i] = 1;
	    myatt.ivalues[isegorder][i] = 1;
	    myatt.ivalues[iup][i] = -1;
	    maxdone++;
	  } else {
	    myatt.ivalues[ishreve][i] = 0;
	  }
      } // for each arc i
System.out.println(maxdone + " source arcs");

// multiple passes until no more are found to be done.
    while(ndone > 0){
      it++;
      System.out.println("iteration " + it);
      ndone = 0;
      for(int i = 0;i < maxid;i++){  // for each arc
	if(verbose) System.out.println("***Arc " + i + ", flag = " + flags[i]);
	if(flags[i] == it){ // just done: look down for candidate
	if(verbose) System.out.println(i + ", flag=" +flags[i] + ",it=" + it);
	 ii = down[i]; //index of downstream arc
	 if(verbose) System.out.println(i + ", id=" + myatt.id[i] + ",ii=" + ii);
	 if((ii >= 0) && (flags[ii] == 0)){ // not a mouth, and not done
	  nup = 0;
	  for(int j = 0;j< maxid;j++){ // looking for upstream arc
	    if(myatt.fnode[ii] == myatt.tnode[j]){  // found upstream arc
	      if(flags[j] > 0){  // upstream arc is ready
		nup++;
		upi[nup] = j;
	      } else {
		nup = 0;
		break;
	      }
	    } // found upstream arc
	  } // checking every upstream arc
	  if (nup > 0) { // all upstream arcs are ready
	  if(verbose) System.out.println(" Doing " + ii + ", id=" + myatt.id[ii] + ",nup =" + nup);
	    for (int j=1;j <= nup;j++){ // for each upstream arc
	      if(verbose) System.out.println(" j=" + j + ",index=" + upi[j] + ", upid=" + myatt.id[upi[j]]);
	      myatt.ivalues[ishreve][ii] += myatt.ivalues[ishreve][upi[j]];
	      myatt.ivalues[imax][ii] += myatt.ivalues[imax][upi[j]];
	      if(j == 1){
		myatt.ivalues[istrahler][ii] = myatt.ivalues[istrahler][upi[1]];
		myatt.ivalues[isegorder][ii] = myatt.ivalues[isegorder][upi[1]] + 1;
		maxmax = myatt.id[upi[1]];
		myatt.ivalues[iup][ii] = myatt.id[upi[1]];
	      } else {
		if(myatt.ivalues[istrahler][upi[j]] == myatt.ivalues[istrahler][ii])
		  myatt.ivalues[istrahler][ii]++; //not rigorous if 3 upstream
		else if(myatt.ivalues[istrahler][upi[j]] > myatt.ivalues[istrahler][ii])
		  myatt.ivalues[istrahler][ii] = myatt.ivalues[istrahler][upi[j]];
		if(myatt.ivalues[isegorder][upi[j]] >= myatt.ivalues[isegorder][ii])
		  myatt.ivalues[isegorder][ii] = myatt.ivalues[isegorder][upi[j]] + 1;
		if(myatt.ivalues[imax][upi[j]] > maxmax){
		  maxmax = myatt.ivalues[imax][upi[j]];
		  myatt.ivalues[iup][ii] = myatt.id[upi[j]];
		}
	      }
	    } // for each upstream arc
	    flags[ii] = it + 1;
	    ndone++;
	    maxdone++;
	  } // all upstream arcs are ready
	 } // not a mouth
	} //just done: look down for candidate
      }			// each arc
      System.out.println(ndone + " more arcs, " + maxdone + " total. ");
    }			// while

// None completed this round.  Should be done.
    if(maxdone != maxid)
      System.out.println("-------------------Hey, " + maxdone + " of " + maxid + " done!");


    for(int i = 0;i < maxid;i++){
	myatt.ivalues[imax][i] *= area;	// cells to square meters
	myatt.ivalues[imean][i] = myatt.ivalues[imax][i] - (myatt.ivalues[ilocal][i] * halfarea);
    }
    System.out.println("Writing to the AAT file.");
    myatt.writeAat();
  }
}
