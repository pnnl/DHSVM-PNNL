      INTEGER FUNCTION Julian (year,month,day)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).
      implicit none
      integer,intent(in) :: year,month,day
      integer:: i,j,k
!
      i=year
      j=month
      k=day
!
      Julian=k-32075+1461*(i+4800+(j-14)/12)/4+367*(j-2-(j-14)/12*12) &
            /12-3*((i+4900+(j-14)/12)/100)/4
!
      return
      END function Julian
