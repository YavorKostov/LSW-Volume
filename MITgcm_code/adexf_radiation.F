C $Header ykk Exp $
C $Name:   $

#include "EXF_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: ADEXF_RADIATION
C     !INTERFACE:
      SUBROUTINE ADEXF_RADIATION( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE ADEXF_RADIATION
C     *==========================================================*
C     \ev
C
C

C     !USES:
       IMPLICIT NONE
C     === Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"

#include "EXF_PARAM.h"
#include "EXF_FIELDS.h"
#include "EXF_CONSTANTS.h"

#ifdef ALLOW_AUTODIFF_TAMC
#include "tamc.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     input:
C     myTime  :: Current time in simulation
C     myIter  :: Current iteration number in simulation
C     myThid  :: My Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
C     output:
CEOP


      RETURN
      END
