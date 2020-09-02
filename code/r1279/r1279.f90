!*******************************************************************************
	! SUBROUTINE ir1279vector(random_vector, no_of_rand)
  !
	! USE	 r1279block
	! IMPLICIT NONE
	! INTEGER	 no_of_rand, i, random_vector(*)
	! REAL	 ran2
  !
	! do i = 1, no_of_rand
	!     ioffset = iand(ioffset + 1, 2047)
	!     irand(ioffset) = (irand(index1(ioffset))*irand(index2(ioffset)))
	!     random_vector(i) = lshift(irand(ioffset), -1)
	! end do
  !
	! END

!*******************************************************************************
FUNCTION r1279()

    IMPLICIT NONE
    INCLUDE "r1279block.h"
    REAL    r1279, inv_max
    REAL    INV_MAXINT
    PARAMETER (INV_MAXINT = 1.0/2147483647.)

    ioffset = iand(ioffset + 1, 2047)
    irand(ioffset) = (irand(index1(ioffset))*irand(index2(ioffset)))
    r1279 = ishft(irand(ioffset), -1) * INV_MAXINT

END

!*******************************************************************************
FUNCTION ir1279()

    IMPLICIT NONE
    INCLUDE  "r1279block.h"
    INTEGER  ir1279

    ioffset = iand(ioffset + 1, 2047)
    irand(ioffset) = (irand(index1(ioffset))*irand(index2(ioffset)))
    ir1279 = ishft(irand(ioffset), -1)

END

!*******************************************************************************
FUNCTION ir1279range(imin, imax)

    IMPLICIT NONE
    INCLUDE "r1279block.h"
    INTEGER ir1279range, imin, imax
    REAL INV_MAXINT
    PARAMETER (INV_MAXINT = 1.0/2147483647.)
    REAL    range

    range = float(imax - imin +1) * INV_MAXINT
    ioffset = iand(ioffset + 1, 2047)
    irand(ioffset) = (irand(index1(ioffset))*irand(index2(ioffset)))
    ir1279range = imin + int(ishft(irand(ioffset), -1) * range)
    if  (ir1279range >  imax) ir1279range = imax

END

!*******************************************************************************
SUBROUTINE setr1279(iseed)

    IMPLICIT	NONE
    INCLUDE "r1279block.h"
    INTEGER	ibit, ispoke, one_bit, iseed, localseed, NBITM1
    REAL	ran2
    PARAMETER (NBITM1 = 31)
!
!	Initialize ioffset. This will be increased by (1 mod 2048) for
!	each random number which is called.
!
    ioffset = 0
!
!	Set up the two arrays which give locations of the two random
!	numbers which will be multiplied to get the new random number
!
    do ispoke = 0, 2047
	index1(ispoke) = iand(ispoke - 1279, 2047)
	index2(ispoke) = iand(ispoke - 418, 2047)
    end do
!
!	set up the initial array of 2048 integer random numbers
!	Each bit is separately initialized using ran2 from numerical recipes
!
    localseed = -abs(iseed)


! Matteo:  I have substituted lshift with ishft


    do ispoke = 0, 2047

	irand(ispoke) = 0
	do ibit = 0, NBITM1
	    one_bit = 0
	    if (ran2(localseed) > 0.5) one_bit = 1
	    irand(ispoke) = ior(irand(ispoke), ishft(one_bit, ibit))
	end do
	irand(ispoke) = 2 * irand(ispoke) + 1

    end do

END

!*******************************************************************************
! SUBROUTINE setseed(iseed)
!     IMPLICIT NONE
!     integer iseed, idum , secs, mygetpid, iran1
!
!     iseed = abs(mygetpid() * secs())
!     idum = -iseed
!     iseed = abs (iseed * iran1(idum) )
!
! END
!*******************************************************************************

      FUNCTION iran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,    &
      NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      iran1=iy
      return
      END
