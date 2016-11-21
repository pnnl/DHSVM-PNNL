c
c      PROGRAM RMAIN
C
C     Dynamic river basin model for simulating water quality in
C     branching river systems with freely-flowing river segments. 
c
c     This version uses Reverse Particle Tracking in the Lagrangian
c     mode and Lagrangian interpolation in the Eulerian mode.
c
c     Topology and routing is set up to be consistent with output
c     from the Distributed Hydrologic Soil and Vegetation Model (DHSVM)
c     model developed by the Land Surface Hydrology Group at 
c     the University of Washington.
c
C     For additional information visit:
c
c     http://www.hydro.washington.edu/Lettenmaier/Models/DHSVM/
c
c     or contact:
c
C     John Yearsley
C     Land Surface Hydrology Group
C     Dept. of Civil and Environmental Engineering
C     Box 352700
C     University of Washington
C     Seattle, Washington
C     98195-2700
C     yearsley@hydro.washington.edu
C
      character*8  start_data,end_data     
      character*200 Prefix
      integer iargc
      integer numarg
 
c     Command line input
c
      numarg = iargc ( )
      if (numarg .lt. 1) then
        write (*,*) 'Too few arguments were given'
        write (*,*) ' '
        write (*,*) 'First:  Location and prefix of input file'
        write (*,*) '        (networkfile)'
        write (*,*) 'eg: $ <program-name> <Project Name>'
        write (*,*) ' '
        stop
      end if
      call getarg ( 1, Prefix )
c
c     Identify and open necessary files
c
c     open the output file 
      open(unit=20,file=TRIM(Prefix)//'.temp',status='unknown')
c
c     Open file with weather and inflow data
      write(*,*) 'Forcing file -  ', TRIM(Prefix)//'.forcing'
      open(unit=30,file=TRIM(Prefix)//'.forcing',STATUS='old')
C
c     open Mohseni file 
      open(40,file=TRIM(Prefix)//'.Mohseni',STATUS='old')    
c
c     open Leopold file 
      open(50,file=TRIM(Prefix)//'.Leopold',STATUS='old')    
c
c     Open network file
      OPEN(UNIT=90,FILE=TRIM(Prefix)//'.net',STATUS='OLD')
c
c     Read header information from control file
      do n=1,5
        read(90,*)
      end do
c
C     Call systems programs to get started
C
C     SUBROUTINE BEGIN reads control file, sets up topology and
C     important properties of reaches
      write(*,*) 'Calling BEGIN'
      CALL BEGIN
C
C     SUBROUTINE SYSTMM performs the simulations
C
      CALL SYSTMM
C
C     Close files after simulation is complete
C
      write(*,*) ' Closing files after simulation'
      CLOSE(30)
      CLOSE(90)
      STOP
      END
      SUBROUTINE BEGIN
      character*11 end_time,start_time
      character*5 Dummy_B
      character*10 Dummy_A
      integer head_name,trib_cell,first_cell
      dimension ndmo(12)
      logical Test
      INCLUDE 'RBM.fi'
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
      ndelta=2
      delta_n=ndelta
c
c     Read the starting and ending times and the number of
c     periods per day of weather data from the forcing file
c
      read(30,*) start_time,end_time,nwpd,nd_start
      write(*,*) start_time,'  ',end_time,nwpd,nd_start
c 
      write(*,*) 'Number of simulations per day - ',nwpd
c
      read(start_time,'(i4,2i2,1x,i2)') start_year,start_month
     &                                 ,start_day,start_hour
      read(end_time,'(i4,2i2,1x,i2)') end_year,end_month
     &                               ,end_day,end_hour
      nyear1=start_year
      nyear2=end_year
      write(*,*) start_year,start_month,start_day
     &             ,end_year,end_month,end_day
c
c     Establish the Julian day for which simulations begin
c
      jul_start=julian(start_year,start_month,start_day)
c
      read(90,*) no_rch
      write(*,*) "Number of stream reaches - ",no_rch
      read(40,*) Test,a_smooth
      b_smooth=1.- a_smooth
      do nr=1,no_rch
        read(40,*) nnr,alf_Mu(nr),beta(nr),gmma(nr),mu(nr)
      end do
c
c     Start reading the reach date and initialize the reach index, NR
c     and the cell index, NCELL
c
      ncell=0
      nreach=0
      ns_total=0
  100 continue
C
C     Card Group IIb. Reach characteristics
C
      nreach=nreach+1
      if(nreach.gt.no_rch) go to 500
c
c     Initialize NSEG, the total number of segments in this reach
      nseg=0
      write(*,*) ' Starting to read reach ',nreach
c
c     Read the number of cells in this reach, the headwater #,
c     the number of the cell where it enters the next higher order stream,
c     the headwater number of the next higher order stream it enters, and
c     the river mile of the headwaters.
c
c
      read(90,*) Dummy_A,no_cells(nreach),
     &              Dummy_A,main_stem(nreach),Dummy_A,trib_cell                
c
c     If this is reach that is tributary to cell TRIB_CELL, give it the
c     pointer TRIB(TRIB_CELL) the index of this reach for further use.
c     Also keep track of the total number of tributaries for this cell
c
      if (trib_cell.gt.0) then
         no_tribs(trib_cell)=no_tribs(trib_cell)+1
         trib(trib_cell,no_tribs(trib_cell))=nreach
      end if
c
c     Reading Reach Element information
c
c
      first_cell=1
c
      do nc=1,no_cells(nreach)
        ncell=ncell+1
        read(50,*) ncll,U_a(ncell),U_b(ncell),U_min(ncell)
        if (ncll .ne.ncell) then
          write(*,*) 'Mismatch in Leopold file'
        end if  
        read(50,*) D_a(ncell),D_b(ncell),D_min(ncell)
        write(*,*) ' Starting to read reach ',nreach
     &      ,U_a(ncell),U_b(ncell),U_min(ncell)

c
c     The headwaters index for each cell in this reach is given
c     in the order the cells are read*dt_calc/(z*rfac)
c
C     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
C

        read(90,*) Dummy_B,nnc,Dummy_B,node(nnc)
     &           ,Dummy_B,x_0,Dummy_B,x_1
     &           ,Dummy_B,elev(nreach)
        if (first_cell.eq.1) then
          first_cell=0
          head_cell(nreach)=ncell
          x_dist(nreach,0)=x_1
        end if
        dx(ncell)=(x_1-x_0)/ndelta
        nndlta=0
  200 continue
c
        nndlta=nndlta+1
        nseg=nseg+1
        ns_total=ns_total+1
        segment_cell(nreach,nseg)=ncell
        x_dist(nreach,nseg)=x_dist(nreach,nseg-1)-dx(ncell)
c 
        if (nndlta.lt.ndelta) go to 200   
        no_celm(nreach)=nseg
        segment_cell(nreach,nseg)=ncell
        x_dist(nreach,nseg)=x_0
C
      end do
c
c     Define last segment for purposes of specifying tributary
c     temperatures 11/14/2008
c
c
	go to 100
  500	continue
      nreach=no_rch
      xwpd=nwpd
      dt_comp=86400./xwpd
C
C     ******************************************************
C                         Return to RMAIN
C     ******************************************************
C
c
      RETURN
  900 END
      SUBROUTINE SYSTMM
      real*4 xa(4),ta(4),T_head(1000),T_smth(1000)
     *      ,dt_part(1000),x_part(1000)
      real*8 day_fract,hr_fract,sim_incr,year,prnt_time
      integer no_dt(1000),nstrt_elm(1000)
     .     ,ndltp(4),nterp(4),nptest(4),ndmo(12,2)
      logical DONE

      INCLUDE 'RBM.fi'
      data ndltp/-2,-1,-2,-2/,nterp/4,3,2,3/
      data lat/47.6/,pi/3.14159/,rfac/304.8/
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334
     &         ,0,31,60,91,121,152,182,213,244,274,305,335/
c
c
      hour_inc=1./nwpd
      do nr=1,1000
         T_head(nr)=mu(nr)
         T_smth(nr)=mu(nr)
      end do
      n1=1
      n2=2
      nobs=0
c
c     Initialize the day counter used for calculating the
c     record position for direct access files
c
c
      write(*,*) 'Start_year ',start_year,start_month,start_day
      write(*,*) 'End_year   ',end_year,end_month,end_day
      ndays=-Julian(start_year,start_month,start_day)
     &     +Julian(end_year,end_month,end_day)+1
      write(*,*) 'Number of days ',ndays
c      
c     Setup the timing of the simulation
c 
      lp_year=1
      yr_days=365.
      if (mod(start_year,4).eq.0) then
        lp_year=2
        yr_days=366.
      end if
      day_fract=ndmo(start_month,lp_year)+start_day-1
      write(*,*) 'day_fact day ',day_fract
      day_fract=day_fract/yr_days
      hr_fract=start_hour/(24.*yr_days)
      sim_incr=dt_comp/(365.*86400.)
      year=start_year
      time=year+day_fract+hr_fract
c
c     Year loop starts
      write(*,*) start_year,start_day,start_hour,time,day_fract,hr_fract
      do nyear=start_year,end_year
         write(*,*) ' Simulation Year - ',nyear
         nd_year=365
         lp_year=1
         if (mod(nyear,4).eq.0) then 
           nd_year=366
           lp_year=2
         end if
         xd_year=nd_year
c
c        Day loop starts
         if (nyear.eq.start_year) then
           nd1=ndmo(start_month,lp_year)+start_day
         else
           nd1=1
         end if
         if (nyear.eq.end_year) then
           nd2=ndmo(end_month,lp_year)+end_day
         else
           nd2=nd_year
         end if
c
c
         DO ND=nd1,nd2
c
c     Start the numbers of days-to-date counter
           ndays=ndays+1
c
c     Daily period loop starts
           DO ndd=nd_start,nwpd 
c
c     Begin reach computations
c      
             ind=nd
             ipd=ndd
             day=nd
             period=ndd
             time=year+(day-1.+hour_inc*period)/xd_year
c
c     Read advected energy and meteorology data      
             l_seg = 0
             do nr=1,nreach
c
c     Hardwire annual average temperature for headwaters 
c
               do nc=1,no_cells(nr)
                 l_seg=l_seg+1
                 read(30,*,end=900) l1
     &                      ,press(l_seg),dbt(l_seg)
     &                      ,qna(l_seg),qns(l_seg),ea(l_seg),wind(l_seg)
     &                      ,qin(l_seg),qout(l_seg)
                 if (qin(l_seg) < 0.5) then
                     qin(l_seg)=qout(l_seg)
                 end if
                 qavg=0.5*(qin(l_seg)+qout(l_seg))
c
c    Stream speed estimated with Leopold coefficients
c 
                 u(l_seg)=U_a(l_seg)*(qavg**U_b(l_seg)) 
                 u(l_seg) = amax1(u_min(l_seg),u(l_seg))
c 
                 qdiff(l_seg)=(qout(l_seg)-qin(l_seg))/delta_n
c 
                 dt(l_seg)=dx(l_seg)/u(l_seg)
c
c    Depth estimated with Leopold coefficients
c 
                 depth(l_seg)=D_a(l_seg)*(qavg**D_b(l_seg))
                 depth(l_seg)=amax1(D_min(l_seg),depth(l_seg))
c
                 lat_flow(l_seg)=qout(l_seg)-qin(l_seg)
               end do
c
c  Set the value of the tributary flow due to the reach, NR
c
               q_trib(nr)=qout(l_seg)
	     end do
c
c     Main stem inflows and outflows for each reach first
c     Flows are cumulative and do not include tributaries if
c     tributaries are downstream of the inflow junction
c
c     Headwaters flow and temperature
c
 90            continue
c
c     Begin cycling through the reaches 
c
               do nr=1,nreach
                  nc_head=segment_cell(nr,1)
                  T_smth(nr)=b_smooth*T_smth(nr)+a_smooth*dbt(nc_head)
                  T_head(nr)=mu(nr)
     &            +(alf_Mu(nr)/(1.+exp(gmma(nr)*(beta(nr)-T_smth(nr)))))
c                  
                  temp(nr,0,n1)=T_head(nr)
                  temp(nr,-1,n1)=T_head(nr)
                  temp(nr,-2,n1)=T_head(nr)
                  temp(nr,no_celm(nr)+1,n1)=temp(nr,no_celm(nr),n1)
                  x_head=x_dist(nr,0)
                  x_bndry=x_head-1.0

c     First do the reverse particle tracking
c

                  do ns=no_celm(nr),1,-1
c
c     Segment is in cell SEGMENT_CELL(NC)
c

                     ncell=segment_cell(nr,ns)
                     nx_s=1
                     nx_part=ns
                     dt_part(ns)=dt(ncell)
                     dt_total=dt_part(ns)
                     x_part(ns)=x_dist(nr,ns)
 100                 continue
c
c     Determine if the total elapsed travel time is equal to the
c     computational interval
c

                     if(dt_total.lt.dt_comp) then
                        x_part(ns)=x_part(ns)
     .				         +dx(segment_cell(nr,nx_part))
c     If the particle has started upstream from the boundary point, give it
c     the value of the boundary
c

                        if(x_part(ns).ge.x_bndry) then
                           x_part(ns)=x_head
                           dt_part(ns)=dt(segment_cell(nr,nx_part))
                           dt_total=dt_total+dt_part(ns)
c                           nx_part=head_cell(nr)
                           go to 200
                        end if
c
c     Increment the segment counter if the total time is less than the
c     computational interval
c
                        nx_s=nx_s+1
                        nx_part=nx_part-1
                        dt_part(ns)=dt(segment_cell(nr,nx_part))
                        dt_total=dt_total+dt_part(ns)
                        go to 100
                     else
c
c     For the last segment of particle travel, adjust the particle location
c     such that the total particle travel time is equal to the computational
c     interval.
c

                        dt_before=dt_part(ns)
                        dt_part(ns)
     .                       =dt_comp-dt_total+dt_part(ns)
                        x_part(ns)=x_part(ns)
     .                            +u(segment_cell(nr,nx_part))
     .                            *dt_part(ns)
                        if(x_part(ns).ge.x_head) then
                           x_part(ns)=x_head
                           nx_s=nx_s-1
                           dt_part(ns)=dt(head_cell(nr))
                        end if
                     end if
 200                 continue
                     if(nx_part.lt.1) nx_part=1
                     nstrt_elm(ns)=nx_part
                     no_dt(ns)=nx_s
                  end do
                  DONE=.FALSE.
                  do ns=1,no_celm(nr)
                     ncell=segment_cell(nr,ns)
                     itest=no_celm(nr)
c
c     Net solar radiation (kcal/meter^2/second)
c
c     qns
c
c     Net atmospheric radiation (kcal/meter^2/second)
c
c     qna
c
c     Dry bulb temperature (deg C)
c
c     dbt
c
c     Wind speed (meters/second)
c
c     wind
c      time=day_fract+hr_fract
c     Factor for Bowen ratio ((deg C)^-1)
c
c     pf
c
c     Vapor pressure at given air temperature (mb)
c
c     ea
c
c     Photo period (fraction of a day.  Not used in the energy budget)
c
c     phper
c


 250                 continue
c
c     Now do the third-order interpolation to
c     establish the starting temperature values
c     for each parcel
c
                     nseg=nstrt_elm(ns)
                     npndx=1
c
c     If starting element is the first one, then set
c     the initial temperature to the boundary value
c
                     if (nseg.eq.1) then
                        t0=T_head(nr)
                        go to 350
                     end if
c
c     Perform polynomial interpolation
c
                     do ntrp=1,nterp(npndx)
                        npart=nseg+ntrp+ndltp(npndx)-1
                        nptest(ntrp)=npart
                        xa(ntrp)=x_dist(nr,npart)
                        ta(ntrp)=temp(nr,npart,n1)
                     end do
                     x=x_part(ns)
  280                continue
c
c     Call the interpolation function
c

                     t0=tntrp(xa,ta,x,nterp(npndx))
                     ttrp=t0
 300                 continue
 350                 continue
                     dt_calc=dt_part(ns)
                     nncell=segment_cell(nr,nstrt_elm(ns))
c
c    Set NCELL0 for purposes of tributary input
c
                     ncell0=nncell
                     dt_total=dt_calc
                     do nm=no_dt(ns),1,-1
                       u_river=u(nncell)/3.2808
                       z=depth(nncell)
                       call energy
     &                      (t0,QSURF,A,B,ncell)
                       t_eq=-B/A
                       qdot=qsurf/(z*rfac)
                       t0=t0+qdot*dt_calc  

                       if(t0.lt.0.0) t0=0.0
 400                   continue
c
c     Look for a tributary.
c 
                       q1=qin(nncell)

                       ntribs=no_tribs(nncell)
                       if (ntribs.gt.0.and..not.DONE) then
                         do ntrb=1,ntribs
                           nr_trib=trib(nncell,ntrb)
                           q2=q1+q_trib(nr_trib)
                           t0=(q1*t0+q_trib(nr_trib)*T_trib(nr_trib))/q2
                           q1=q1+q_trib(nr_trib)
c
  450 continue
                         end do
                         DONE=.TRUE.
                       end if
                       t00=t0
                       if (lat_flow(nncell).gt.0) then
                          q1=0.5*(qin(nncell)+qout(nncell))
                          q2=q1+lat_flow(nncell)
c
c  Modified nonpoint source temperature so as to be the same
c  as the instream simulated temperature for Connecticut River 7/2015
                          T_dist=t0
                          t0=(q1*t0+lat_flow(nncell)*T_dist)/q2
                          dtlat=t0-t00
                        end if
 500                    continue
                        nseg=nseg+1
                        nncell=segment_cell(nr,nseg)
c
c     Reset tributary flag is this is a new cell
c
                        if (ncell0.ne.nncell) then
                           ncell0=nncell
                           DONE=.FALSE.
                        end if
                        dt_calc=dt(nncell)
                        dt_total=dt_total+dt_calc
                      end do
                      if (t0.lt.0.5) t0=0.5
                     temp(nr,ns,n2)=t0
	             T_trib(nr)=t0
c
c   Write file 20 with all temperature output 11/19/2008
c
                     nsmod=mod(ns,2)
                     time=year+(day-1.+hour_inc*period)/xd_year
                     if(nsmod.eq.0) then
                       rmile_plot=x_dist(nr,ns)/5280.
                       write(20,'(f11.5,i5,1x,i4,1x,2i5,1x,5f7.2,f9.2)') 
     &                       time,nyear,nd,ncell,ns,t0
     &                      ,T_head(nr),dbt(ncell)
     &                      ,depth(ncell),u(ncell),qin(ncell)

                     end if
c
c     End of computational element loop
c

                  end do
c     End of reach loop
c

               end do
               ntmp=n1
               n1=n2
               n2=ntmp
c
c     End of weather period loop (NDD=1,NWPD)
c
            end do
c
c Reset daily loop counter
c
          nd_start=1
c 
C
c     End of main loop (ND=1,365/366)
c

         end do
c
c    Update initial time for new year
c
      year=year+1
c
c     End of year loop
c
      end do
c
c Finish
c
  900 Continue
c
c
c     ******************************************************
c                        return to rmain
c     ******************************************************
c

  950 return
      end
      SUBROUTINE ENERGY
     &           (TSURF,QSURF,A,B,ncell)
      REAL*4 Ksw,LVP
      real*4 q_fit(2),T_fit(2),evrate
      INCLUDE 'RBM.fi'
      data evrate/1.5e-9/,pf/0.640/
      parameter (pi=3.14159)
c     
      td=nd
      evap_rate=evrate
      T_fit(1)=tsurf-0.5
      T_fit(2)=tsurf+0.5
      do i=1,2
         E0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
         RB=PF*(DBT(ncell)-T_fit(i))
         LVP=597.0-0.57*T_fit(i)
         QEVAP=1000.*LVP*evap_rate*WIND(ncell)
         if(qevap.lt.0.0) qevap=0.0
         QCONV=RB*QEVAP
         QEVAP=QEVAP*(E0-EA(ncell))
         QWS=6.693E-2+1.471E-3*T_fit(i)
         q_fit(i)=QNS(ncell)+0.97*QNA(ncell)-QWS-QEVAP+QCONV
      end do

      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))
     .     /(T_fit(1)-T_fit(2))
      qsurf=0.5*(q_fit(1)+q_fit(2))
C
C     ******************************************************
C               Return to Subroutine RIVMOD
C     ******************************************************
C
      RETURN
      END
      function nodays(jtime,jy0)
      dimension ndmo(12)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
      jy=(jtime/10000)
      jrem=(jtime-jy*10000)
      jm=jrem/100
      jd=jrem-jm*100
      ny=jy-jy0
      nodays=365*ny+ndmo(jm)+jd
      return
      end
      function ndate(n,ny0)
      dimension ndmo(13)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334,365/
      ny=(n-1)/365
      njul=n-ny*365
      nm=1
 10   continue
      if(ndmo(nm+1).ge.njul) go to 50
      nm=nm+1
      go to 10
 50   continue
      nday=njul-ndmo(nm)
      nmon=100*nm
      nyear=10000*(ny0+ny)
      ndate=nyear+nmon+nday
      return
      end
c
c	Third-order polynomial interpolation using Lagrange
c     polynomials.  FUNCTION is SUBROUTINE POLINT from
c     Numerial Recipes
c
      FUNCTION tntrp(XA,YA,X,n)
c      PARAMETER (N=4)
c      DIMENSION XA(N),YA(N),C(N),D(N)
      DIMENSION XA(4),YA(4),C(4),D(4)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) DEN=0.001
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
	  tntrp=y
      RETURN
      END
c      function julian(iy,id)
c	  julian=10000*iy+id
c	  return
c	  end
c
c
      INTEGER FUNCTION Julian (YEAR,MONTH,DAY)
C
C---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
C   DATE (YEAR,MONTH,DAY).
C
      INTEGER YEAR,MONTH,DAY,I,J,K
C

      I= YEAR
      J= MONTH
      K= DAY
      write(*,*) 'Julian ',i,j,k
C

      Julian=
     1   K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)
     2  /12-3*((I+4900+(J-14)/12)/100)/4
C

      RETURN
      END









