Program Create_File
implicit none
!
! Integer variables
!
integer::delta_t,ierror,iostat,n,nc,nd_start,no_cycles
integer::narray,nf,nfile,n_head,nn,no_years,no_seg
integer::no_dt,no_days,nobs_start,nobs_end
integer::start_day,start_mon,start_yr,end_day,end_mon,end_yr,start_hour,end_hour
integer::Julian,start_jul,end_jul
integer,allocatable,dimension(:):: seg_no,seg_indx,seg_seq,seg_net
integer,allocatable,dimension(:):: dummy
!
! Real variables
real::press=1013.
real,allocatable,dimension(:)::depth,out_flow,in_flow,lat_flow
real,allocatable,dimension(:,:)::forcing
!
! Character variables
!
character (len=1)  :: colon=':'
character (len=19) :: start_date,end_date 
character (len=19) :: time_stamp0,time_stamp
character (len=4)  :: path
character (len=6)  :: fluff
character (len=8)  :: sequence
character (len=200):: InDirectry,Project
!
!
integer iargc
integer numarg

!
! Command line input
!
      numarg = iargc ( )
      if (numarg .lt. 2) then
        write (*,*) 'Too few arguments were given'
        write (*,*) ' '
        write (*,*) 'First:  Directory with forcing files *(*.Only)'
        write (*,*) 'Second:  Project Name'
        write (*,*) 'eg: $ ./Create_File <Input directory> <Project Name>'
        write (*,*) ' '
        stop
      end if
      call getarg ( 1, InDirectry )
      call getarg ( 2, Project )
!
!write(*,*) 'Name of Project.  Required files include:'
!write(*,*) 'ProjectName.map'
!write(*,*) 'ATP.Only     - Air temperature'
!write(*,*) 'NLW.Only     - Net longwave radiation'
!write(*,*) 'NSW.Only     - Net shortwave radiation'
!write(*,*) 'VP.Only      - Vapor pressure'
!write(*,*) 'WND.Only     - Wind speed'
!write(*,*) 'Inflow.Only  - Segment inflow'
!write(*,*) 'Outflow.Only - Segment outflow'
!
!read(*,*) Project
!
write(*,*) TRIM(Project)//'.map'
open(10,file=TRIM(Project)//'.segmap',status='old')
open(20,file=TRIM(InDirectry)//'/ATP.Only',status='old')
open(21,file=TRIM(InDirectry)//'/NLW.Only',status='old')
open(22,file=TRIM(InDirectry)//'/NSW.Only',status='old')
open(23,file=TRIM(InDirectry)//'/VP.Only',status='old')
open(24,file=TRIM(InDirectry)//'/WND.Only',status='old')
open(25,file=TRIM(InDirectry)//'/Inflow.Only',status='old')
open(26,file=TRIM(InDirectry)//'/Outflow.Only',status='old')
!
open(30,file=TRIM(Project)//'.forcing',status='unknown')
!
! Get some information about the network from the "<Project>.map" file
!
read(10,*) n_head,no_seg
!
! Allocate arrays after increasing expected size to account
! for differences in the number of segments and the size of
! the largest stream segment number in DHSVM.
!
narray=no_seg+no_seg/2
allocate (dummy(narray))
allocate (seg_no(narray))
allocate (seg_indx(narray))
allocate (seg_net(narray))
allocate (seg_seq(narray))
allocate (in_flow(narray))
allocate (out_flow(narray))
!allocate (lat_flow(narray))
!allocate (depth(narray))
allocate (forcing(5,narray))
!
do n=1,no_seg
  read(10,*) sequence,nn,path,seg_no(n)
end do
!
nfile=20
do nf=1,7
  read(nfile,'(A19,1x,A19,1x,i2)') start_date,end_date, no_dt
  write(*,*) 'start ',start_date
  write(*,*) 'end ',end_date
  read(start_date,'(i2,1x,i2,1x,i4,1x,i2,a6)') start_mon,start_day,start_yr  &
                                              ,start_hour,fluff
  start_jul=Julian(start_yr,start_mon,start_day)
  read(end_date,'(i2,1x,i2,1x,i4,1x,i2,a6)') end_mon,end_day,end_yr          &
                                              ,end_hour,fluff
  nfile=nfile+1
end do
!
!
delta_t=no_dt
!
! Determine number of time steps per day
write(*,*) 'delta_t ',no_dt,delta_t
no_dt=24/delta_t
!
nobs_start=no_dt
nobs_end=no_dt
!
if (no_dt .gt. 1) then
  nobs_start = (24 - start_hour)/delta_t
  nobs_end   = 1+(end_hour/delta_t)
end if 
nd_start=1+(start_hour/delta_t)
if (nobs_start .lt. no_dt) then
  write(*,*)
  write(*,*) '!!!===WARNING: There is (are) only',nobs_start,'value(s) for the first simulation day'
  write(*,*) 'rather than ',no_dt, 'values. Daily averaging scripts will exclude the first' 
  write(*,*) 'day in the daily temperature output. ===!!!'
  write(*,*)
end if
if (nobs_end .lt. no_dt) then
  write(*,*)
  write(*,*) '!!!===WARNING: There is (are) only',nobs_end,'value(s) for the last simulation day'
  write(*,*) 'rather than ',no_dt, 'values. Daily averaging scripts will exclude the last'
  write(*,*) 'day in the daily temperature output. ===!!!'
  write(*,*)
end if
!

!
! Julian date of start and end of forcing records
!
end_jul=Julian(end_yr,end_mon,end_day)
no_days=end_jul - start_jul
!
write(*,*) 'Julian',start_jul,end_jul
!
no_cycles=nobs_start+nobs_end+no_dt*(no_days-1)

!
write(*,*) 'no_cycles ',no_cycles
write(30,'(2(i4.4,i2.2,i2.2,a1,i2.2,1x),2i4)')          &
     start_yr,start_mon,start_day,colon,start_hour          &
    ,end_yr,end_mon,end_day,colon,end_hour                  &
    ,no_dt,nd_start
!
! Read segment mapping
!
nfile=19
!
! Read the segment sequencing from the header of the ATP.Only file
! and establish the relationship between DHSVM segment numbers and
! the indexed location in the forcing files
!
  nfile=nfile+1
  read(nfile,*) (seg_seq(n),n=1,no_seg)
  do nf=1,no_seg
    seg_net(seg_seq(nf))=nf
  end do
do nf=2,7
  nfile=nfile+1
  read(nfile,*) (dummy(n),n=1,no_seg)
!  read(nfile,*) time_stamp0
!  write(*,"('Initial Time Stamp - ',a19)"), time_stamp0
end do!
! Read the forcings from the DHSVM file
!
do nc=1,no_cycles
  nfile=19
  do nf=1,5
    nfile=nfile+1
    read(nfile,*) time_stamp,(forcing(nf,n),n=1,no_seg)
  end do
  do n=1,no_seg
    forcing(4,n)=0.01*forcing(4,n)
    forcing(2,n)=2.3884e-04*forcing(2,n)
    forcing(3,n)=2.3884e-04*forcing(3,n)
  end do
!
! Read the streamflow from the DHSVM file
!
  read(25,*) time_stamp,(in_flow(n),n=1,no_seg)
  read(26,*) time_stamp,(out_flow(n),n=1,no_seg)
  do n=1,no_seg
    ! Convert the unit from cubic meter per sec to cubic feet per sec
    in_flow(n) = in_flow(n) * 35.315;
    out_flow(n) = out_flow(n) * 35.315;
    if(in_flow(n) .lt. 0.01 .and.out_flow(n).lt. 0.01) then 
      in_flow(n)  = 0.01
      out_flow(n) = in_flow(n)
    end if
    if (in_flow(n) .lt. 0.01 .and. out_flow(n) .ge. 0.01) then
      in_flow(n) = 0.01
    end if
  end do
!
! Write the output file
!
  do n=1,no_seg
    nf=seg_no(n)
    nn=seg_net(nf)
    !write(*,*) n,nf,nn
    write(30,*) n,press,(forcing(nf,nn),nf=1,5)                  &
               ,in_flow(nn),out_flow(nn)
  end do
end do
end Program Create_File
