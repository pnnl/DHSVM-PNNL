import java.io.*;
/* The 1.1.4 compiler wants this in a separate file
class AatError extends Exception {
	public AatError(String s) { super(s); }
}
*/

class Aat {
/** Copyleft Harvey Greenberg, UW hgreen@u.washington.edu Oct 3  1997
   This class is equivalent to an ARC/INFO AAT (arc attribute table)
   It handles the default AAT items, and any other items as well.
   It contains a description of the nonstandard items, as well as arrays
   of actual objects.
   The constructor method reads the binary version of the item list.
   Three sets of arrays are maintained to store data:
                     string	integer	double
     character		X
     Integer		X	X
     Numeric		X		X
     Binary int			X
     Float				X
     date	NOT SUPPORTED
     The user has the choice of accessing Integer (in the INFO sense of the word)
      and Numeric data as string or binary, but it is the binary value that gets
      written back to the AAT.
     Reviewing the 1.1.1 docs, I wonder if it's worth implementing BigInteger for I
     and BigDecimal for N.
    VERSION 1.0 Oct 3  1997
    version 1.1 Mar 11 1998 Fixed a serious bug.  The old version scrambled
			    with an odd number of data bytes.

   Caveats:
    Error handling is neither complete nor elegant
    Max character item = 80 bytes.
    long integers are not handled well
    Although this has passed initial testing, it is not guaranteed to
    not trash your AAT file.
    This was compiled and checked under JDK 1.1.1 under Solaris.
   Public Methods:
	Aat (String covername) throws IOException, AatError {  // constructor
        void writeAat() throws IOException {
        void report(){       // first calling method
        void report(int rec){        // second method.  INFO-style numbering starts at 1
        int findItem(String itmName)
*/

  // information on AAT
    String covname ; // File f;
    int nrecords,recordLength,defaultItems=7;
    boolean isDouble, odd,verbose=false;
  // information on extra items
    int nitems,extrabytes;
    short itemLength[],itemType[],totalitems;
    String itemName[];
  // default items
    int fnode[];
    int tnode[];
    int lpoly[];
    int rpoly[];
    double length[];
    int recno[];
    int id[];
  // extra items
    String strings[][];
    int ivalues[][];
    double dvalues[][];
 // This has to be declared here or the throw statement won't work
    private int rec;

    Aat (String covername) throws IOException, AatError {  // constructor
	covname = covername;
	int nread,testlength,anint;
	byte bytebuff[] = new byte[80],testbytes[] = new byte[32];
	byte filler;	// INFO pads the record with a null
	boolean matched,minus;
	short skipflag;
	String itemFile,tmpstring,teststring;
	teststring = (new String(covername + ".AAT")).toUpperCase();
	testbytes = teststring.getBytes();
	testlength = teststring.length();
	File f0 = new File ("info/arc.dir");
	DataInputStream s0 = new DataInputStream(new FileInputStream(f0));
	try { // read file, handling EOF
	    matching: while (true){
		nread = s0.read(bytebuff,0,32); // read a filename
	    if(nread != 32)  // It doesn't seem to throw an exception on EOF
		throw new AatError("Could not find itemfile name in arc.dir");
		matched = true;
		for(int i=0;i < testlength;i++){ // compare 2 bytearrays;
						 //strings give me trouble
		    if(testbytes[i] != bytebuff[i]){
			matched = false;
			break;
		    }
		}
		if (matched){
		    s0.skipBytes(3); // skip over "ARC"
		    nread = s0.read(bytebuff,0,4);
		    itemFile = new String("info/arc" + new String(bytebuff,0,4) + ".nit");
		    // Check that this is an active entry
		    s0.skipBytes(1); // blank that evens out count
		    totalitems = s0.readShort();
		    recordLength = s0.readShort();
		    s0.skipBytes(18); // 16 blanks and short 132?
		    skipflag = s0.readShort();
		    nrecords = s0.readInt();
		    if(skipflag != 0) // This is a deleted record
		      s0.skipBytes(312);  // to next record
		    else{
		      if(verbose)System.out.println("Arc item data is in " + itemFile);
		      break matching;
		    }
		}
			/* for testing
			System.out.print(new String(bytebuff,0,24).substring(0,22));
			s0.skipBytes(8);
			for(int j=0;j < 20;j++){
			  short ll = s0.readShort();
			  System.out.print(ll + ",");
			}
			System.out.println("and");
			s0.skipBytes(300);  // to next record
			*/
		else // not suspected match
		  s0.skipBytes(348);  // to next record
	    } // matching
	} //try
	catch (IOException e) {
	    throw new AatError("Could not find itemfile name in arc.dir");
	}
	s0.close();

 	File f1 = new File (itemFile);
//	nitems = (int)(f1.length() / 144) - defaultItems;
	nitems = (int)(totalitems) - defaultItems;
	itemName = new String[nitems];

	DataInputStream s1 = new DataInputStream(new FileInputStream(f1));
	s1.skipBytes(4*144+16); // through rpoly to length of length
	short ll   = s1.readShort();
	isDouble = ll == 8;
	extrabytes = 0;
	if(nitems > 0){
	  itemName = new String[nitems];
	  itemLength = new short[nitems]; 
	  itemType= new short[nitems];
	  s1.skipBytes(126+(2*144)); // end of default items;
	  for(int i=0;i < nitems;i++){
	    nread = s1.read(bytebuff,0,16);	// read a line of binary item list
	    itemName[i] = new String(bytebuff);
	    itemLength[i] = s1.readShort();
	    extrabytes += itemLength[i];
	    s1.skipBytes(12);
	    itemType[i] = s1.readShort();  // 2=C,3=I,4=N,5=B,6=F
	    s1.skipBytes(112);
	  }
	}
	odd = extrabytes % 2 == 1;
	s1.close();
/* We figured out how to find this in info/arc.dir
	recordLength = (defaultItems * 4) + extrabytes;
	if(isDouble)
	  recordLength += 4;
*/

// We have general info; now read the actual AAT
	File f = new File (covername + "/aat.adf");
	if(verbose)System.out.println("aat file f is " + f);
	DataInputStream s = new DataInputStream(new FileInputStream(f));
/* We figured out how to find this in info/arc.dir
	nrecords = (int)f.length() / recordLength;
*/
	fnode = new int[nrecords];
	tnode = new int[nrecords];
	lpoly = new int[nrecords];
	rpoly = new int[nrecords];
	length = new double[nrecords];
	recno = new int[nrecords];
	id = new int[nrecords];
	for (int itno=1;itno < nitems;itno++){
	    strings = new String[nitems][nrecords];
	    ivalues = new int[nitems][nrecords];
	    dvalues = new double[nitems][nrecords];
	}
	try {
	   for(rec=0;rec < nrecords;rec++){  // unlike arcinfo, number from 0
		fnode [rec] = s.readInt();
		tnode [rec] = s.readInt();
		lpoly [rec] = s.readInt();
		rpoly [rec] = s.readInt();
		if(isDouble)
		  length [rec] = s.readDouble();
		else
		  length [rec] = (double)s.readFloat();
		recno [rec] = s.readInt();
		id [rec] = s.readInt();
		for (int itno=0;itno < nitems;itno++){
		  switch (itemType[itno]) {
		    case 2:	//character
		    	nread = s.read(bytebuff,0,itemLength[itno]);
		    	strings[itno][rec] = new String(bytebuff,0,itemLength[itno]);
			break;
		    case 3:	// I
			nread = s.read(bytebuff,0,itemLength[itno]); 
			strings[itno][rec] = new String(bytebuff,0,itemLength[itno]);
// (chokes on blanks!)		ivalues[itno][rec] = Integer.parseInt(strings[itno][rec]);
			anint = 0; // I have trouble with parseInt // CAVEAT, how about negatives?
			minus = false;
			for (int j = 0; j < nread; j++){
			  if(bytebuff[j] != 32)	// blank
			     if(bytebuff[j] == 45)
				minus = true;
			     else
			    	anint = anint * 10 + (bytebuff[j] - 48);
			}
			if(minus)
			   anint = 0 - anint;
			ivalues[itno][rec] = anint;
			break;
		    case 4:	// Numeric
			nread = s.read(bytebuff,0,itemLength[itno]); 
			strings[itno][rec] = new String(bytebuff,0,itemLength[itno]);
			dvalues[itno][rec] = Double.valueOf(strings[itno][rec]).doubleValue();
			break;
		    case 5:   // Binary int
			if(itemLength[itno] == 4)
			    ivalues[itno][rec] = s.readInt();
			else
			    ivalues[itno][rec] = (int)s.readShort();
			break;
		    case 6:	// Float
			if(itemLength[itno] == 4)
			    dvalues[itno][rec] = (double)s.readFloat();
			else
			    dvalues[itno][rec] = s.readDouble();
			break;
		    default:
			break;
		  } // switch
		}
		if(odd) s.skipBytes(1);  // fixed 3/11/98
	   }
	}
	catch (IOException e) {
	    throw new AatError("Error reading data record " + rec);
	}
    }	// end constructor

    public void writeAat() throws IOException {
	byte bytebuff[] = new byte[80],bytebuf2[] = new byte[16],filler=0;
	int tmp,stringLength,padding;
      DataOutputStream t = new DataOutputStream(new FileOutputStream(new File(covname + "/aat.adf")));
      for(int rec=0;rec < nrecords;rec++){
        t.writeInt(fnode[rec]);
        t.writeInt(tnode[rec]);
        t.writeInt(lpoly[rec]);
        t.writeInt(rpoly[rec]);
	if(isDouble)
	  t.writeDouble(length[rec]);
	else
	  t.writeFloat((float) length[rec]);
        t.writeInt(recno[rec]);
        t.writeInt(id[rec]);
	for (int itno=0;itno < nitems;itno++){ // changed to 0 9/30/97
	  switch (itemType[itno]) {
	    case 2:	//character
	  	bytebuff = strings[itno][rec].getBytes();
        	t.write(bytebuff,0,itemLength[itno]);
	    	break;
            case 3:     // I
            case 4:     // N
		if(itemType[itno] == 3)
		  bytebuff = String.valueOf(ivalues[itno][rec]).getBytes();  //  Use integer value
		else
		  bytebuff = String.valueOf(dvalues[itno][rec]).getBytes();  //  Use double value
		stringLength = bytebuff.length;
		padding = itemLength[itno] - stringLength;
		if(padding > 0){	// right-justify the number
		  for (int i=0;i < stringLength;i++){
		    bytebuf2[i + padding] = bytebuff[i];
		  }
		  for (int i=0;i < padding;i++){
		    bytebuf2[i] = 32;
		  }
        	  t.write(bytebuf2,0,itemLength[itno]);
		}
        	else
		  t.write(bytebuff,0,itemLength[itno]);
	    	break;
	    case 5:   // Binary int
                if(itemLength[itno] == 4)	// int
                    t.writeInt(ivalues[itno][rec]);
                else	// short
                    t.writeShort((short)ivalues[itno][rec]);
                break;
            case 6:     // Float
                if(itemLength[itno] == 4)	// float
		    t.writeFloat((float)dvalues[itno][rec]);
                else	// double
		    t.writeDouble(dvalues[itno][rec]);
                break;
	  }
	}
        if(odd)
	   t.writeByte(filler);
      } // for
    }

    public void report(){	// first calling method
        String typeName[] = {" "," ","Character","I (BCD integer)","Numeric","Binary integer","Float (or double)"};
	System.out.println("AAT has " + nrecords + " records, " + nitems + " extra items, " + extrabytes + " extra bytes");
	if(isDouble)
	   System.out.println("It is double precision.");
	else
	   System.out.println("It is single precision.");
	for(int i=0;i < nitems;i++){
	    System.out.println(i + itemName[i] + " is type " + typeName[itemType[i]] + ", width =" + itemLength[i]);
	}
    }

    public void report(int rec){	// second method.  INFO-style numbering starts at 1
        String typeName[] = {" "," ","Character","I (BCD integer)","Numeric","Binary integer","Float (or double)"};
	rec--;	// We strat our array at zero
	for(int itno=0;itno < nitems;itno++){
	  switch (itemType[itno]) {
	    case 2:     //character
	    	System.out.println(itemName[itno] + ": " + strings[itno][rec]);
		break;
            case 3:     // I
            case 5:   // Binary int
		System.out.println(itemName[itno] + ": " + ivalues[itno][rec]);
		break;
            case 4:     // N
            case 6:   // Float
                System.out.println(itemName[itno] + ": " + dvalues[itno][rec]);
		break;
	  }
	}
    }

    public int findItem(String itmName){
	String typeName[] = {" "," ","Character","I (BCD integer)","Numeric","Binary integer","Float (or double)"};
// If a non-default data item exists, what's its array index
	if(verbose) System.out.println("looking for " + itmName);
	for(int i=0;i < nitems;i++){
	  if((itemName[i].substring(0,itmName.length())).equalsIgnoreCase(itmName)){
	    System.out.println(itmName + " is item " + i + " (indexed from 0), it is type " + typeName[itemType[i]]);
	    return i;
	  }
	}
	System.out.println("The item "+itmName+" is not present in this AAT.");
	return -1;
    }

}
