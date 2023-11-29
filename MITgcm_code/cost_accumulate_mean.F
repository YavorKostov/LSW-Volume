C $Header: YKK modified the routine $
C $Name: checkpoint65x $

#include "COST_OPTIONS.h"

      subroutine cost_accumulate_mean( myThid )
C     *==========================================================*
C     | subroutine cost_accumulate_mean                          |
C     | o accumulate mean state for cost evalualtion             |
C     *==========================================================*
C     |                                                          |
C     *==========================================================*
      IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"

#include "cost.h"

C     == Routine arguments ==
C     myThid - Thread number for this instance of the routine.
      integer bi, bj
      integer myThid

#ifdef ALLOW_COST
C     == Local variables
c      _RL pressRef
      _RL pressRefL
      _RL ygmskLoc
      _RL PVloc
      _RL TatL, SalL, RhoL
      _RL TatU, SalU, RhoU
      _RL SigmaL
C Pressure p = 1000 dbar = 100 bar = 1e7 Pa
c      parameter ( pressRef = 1.0E7 )
      integer i, j, k, tst
      integer itlo,ithi
      integer jtlo,jthi

c      integer max_len_mbufY
c      parameter ( max_len_mbufY = 512 )
c      character*(max_len_mbufY) msgbufY
c      character*(*) squeeze_rightY
c      parameter ( squeeze_rightY = 'R' )
c      integer iount2Y
c      parameter ( iount2Y = 6 )


      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C--   Accumulate the cost function on the tile of this instance
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do j=1,sNy
             do i=1,sNx


              do k = 2, Nr
             ygmskLoc = ygmsk(i,j,k,bi,bj)
             if (ygmskLoc > 0) then

      if ((salt(i,j,k,bi,bj)>0) .and. (salt(i,j,k-1,bi,bj)>0)) then
      pressRefL=(-1)*rC(k)*1.0E4
         TatL =
     &    theta(i,j,k,bi,bj)
         SalL =
     &    salt(i,j,k,bi,bj)
      CALL FIND_RHO_SCALAR(
     I     TatL, SalL, pressRefL,
     O     RhoL,
     I     myThid )
         TatU =
     &    theta(i,j,k-1,bi,bj)
         SalU =
     &    salt(i,j,k-1,bi,bj)
      CALL FIND_RHO_SCALAR(
     I     TatU, SalU, pressRefL,
     O     RhoU,
     I     myThid )

      CALL FIND_RHO_SCALAR(
     I     TatL, SalL, pressRefSig*1.0E4,
     O     SigmaL,
     I     myThid )

      PVloc =         fCori(i,j,bi,bj) 
     &                *maskC(i,j,k,bi,bj)*maskC(i,j,k-1,bi,bj)
     &                *(1/(rC(k)-rC(k-1)))*(-1)
     &                *((RhoL-RhoU)/RhoL)
     &                /1.0E-12 

                cMeanTheta(i,j,k,bi,bj) = cMeanTheta(i,j,k,bi,bj)
     &                +maskC(i,j,k,bi,bj)*maskC(i,j,k-1,bi,bj)
     &                *drF(k)*hfacC(i,j,k,bi,bj)*rA(i,j,bi,bj)
     &                *1.5E7
     &   *((1+tanh((SigmaL-SigmaLowBnd)*100000))/2)
     &   *(1-((1+tanh((SigmaL-SigmaUpBnd)*100000))/2))
     &   *((1+tanh((PVloc-PVlowBnd)*100000))/2)
     &   *(1-((1+tanh((PVloc-PVupBnd)*100000))/2))
     &               /lastinterval*deltaTClock 

      endif
              endif
              end do
            end do
          end do
        end do
      end do


#endif

      RETURN
      END
