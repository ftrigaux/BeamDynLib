MODULE BeamDynLib_Types

   USE BeamDyn
   USE BeamDyn_Subs

   IMPLICIT NONE

   TYPE , PUBLIC :: BD_UsrDataType
      LOGICAL        :: DynamicSolve  ! flag for dynamic or static solve (static:false)
      LOGICAL        :: GlbRotBladeT0 ! Initial blade root orientation is also the GlbRot reference frame

      INTEGER(IntKi)                          :: nxD,nxL;        ! Number of nodes for Displacement and Loads
      !REAL(ReKi),DIMENSION(:,:), ALLOCATABLE  :: x;         ! Nodes coordinates

      REAL(DbKi)     :: t             ! time
      REAL(DbKi)     :: dt            ! time increment
      INTEGER(IntKi) :: nt            ! number of substep
      
      REAL(ReKi)     :: GlbPos(3)        ! Initial vector position 
      REAL(ReKi)     :: RootOri(3,3)     ! DCM of the initial root orientation
      REAL(ReKi)     :: GlbRot(3,3)      ! 
      REAL(ReKi)     :: RootRelInit(3,3)

      REAL(ReKi)     :: theta_rot          ! Angle for the rotation
      REAL(ReKi)     :: omega(3)           ! Angular velocity vector
      REAL(ReKi)     :: dOmega(3)          ! Angular acceleration vector
      REAL(ReKi)     :: PAngInp_rad        ! Pitch angle input command in radians

      REAL(ReKi),DIMENSION(:,:), ALLOCATABLE         :: loads    ! Forces and moment at the node positions; Size(NNodes,6)
      REAL(ReKi)                                     :: grav(3)  ! Gravity vector

      !REAL(DbKi),DIMENSION(:,:), ALLOCATABLE         :: u        ! Displacement variables (u,v,w,phi,th1,th2); Size(NNodes,6)
      !REAL(DbKi),DIMENSION(:,:), ALLOCATABLE         :: du       ! Velocity variables d(u,v,w,phi,th1,th2)/dt; Size(NNodes,6)

      CHARACTER(1024)  :: InputFile      !< Name of the input file; remove if there is no file [-]
      CHARACTER(1024)  :: OutputFile     !< Name of the output file; remove if there is no file [-]

      ! Contains the variables for the BeamDyn Analysis
      TYPE(BD_InitInputType)           :: BD_InitInput
      TYPE(BD_ParameterType)           :: BD_Parameter
      TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
      TYPE(BD_InitOutputType)          :: BD_InitOutput
      TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
      TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
      TYPE(BD_OtherStateType)          :: BD_OtherState
      TYPE(BD_MiscVarType)             :: BD_MiscVar
      TYPE(BD_InputType) ,ALLOCATABLE  :: BD_Input(:)
      REAL(DbKi),         ALLOCATABLE  :: BD_InputTimes(:)
      TYPE(BD_OutputType)              :: BD_Output
      !INTEGER(IntKi)                   :: DvrOut 

      ! Private variables used inside BeamDyn code
      TYPE(MeshType)                                :: RotationCenter
      TYPE(MeshMapType)                             :: Map_RotationCenter_to_RootMotion

   END TYPE

CONTAINS

! This routines packs all the BD_UsrDataType data into array, to be later written on files
 SUBROUTINE BD_PackBeamDyn_Data( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(BD_UsrDataType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i1
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'BD_PackBeamDyn_Data'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat     = ErrID_None
  ErrMsg      = ""
  Re_BufSz    = 0
  Db_BufSz    = 0
  Int_BufSz   = 0

  ! Continuous states
  Int_BufSz   = Int_BufSz   + 1     ! x allocated yes/no
  !IF ( ALLOCATED(InData%BD_ContinuousState) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*2  ! x upper/lower bounds for each dimension
   ! Allocate buffers for subtypes, if any (we'll get sizes from these) 
    !DO i2 = LBOUND(InData%x,2), UBOUND(InData%x,2)
    !DO i1 = LBOUND(InData%x,1), UBOUND(InData%x,1)
        Int_BufSz   = Int_BufSz + 3  ! x: size of buffers for each call to pack subtype
        CALL BD_PackContState( Re_Buf, Db_Buf, Int_Buf, InData%BD_ContinuousState, ErrStat2, ErrMsg2, .TRUE. ) ! x 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! x
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! x
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! x
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
    !END DO
  !END IF

  ! Discrete States
  Int_BufSz   = Int_BufSz   + 1     ! xd allocated yes/no
  !IF ( ALLOCATED(InData%BD_DiscreteState) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*2  ! xd upper/lower bounds for each dimension
    !DO i2 = LBOUND(InData%xd,2), UBOUND(InData%xd,2)
    !DO i1 = LBOUND(InData%xd,1), UBOUND(InData%xd,1)
      Int_BufSz   = Int_BufSz + 3  ! xd: size of buffers for each call to pack subtype
      CALL BD_PackDiscState( Re_Buf, Db_Buf, Int_Buf, InData%BD_DiscreteState, ErrStat2, ErrMsg2, .TRUE. ) ! xd 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! xd
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! xd
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! xd
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
    !END DO
  !END IF

  ! Constrain states
  Int_BufSz   = Int_BufSz   + 1     ! z allocated yes/no
  !IF ( ALLOCATED(InData%BD_ConstraintState) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*2  ! z upper/lower bounds for each dimension
    !DO i2 = LBOUND(InData%z,2), UBOUND(InData%z,2)
    !DO i1 = LBOUND(InData%z,1), UBOUND(InData%z,1)
      Int_BufSz   = Int_BufSz + 3  ! z: size of buffers for each call to pack subtype
      CALL BD_PackConstrState( Re_Buf, Db_Buf, Int_Buf, InData%BD_ConstraintState, ErrStat2, ErrMsg2, .TRUE. ) ! z 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! z
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! z
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! z
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
    !END DO
  !END IF

  ! Other States
  Int_BufSz   = Int_BufSz   + 1     ! OtherSt allocated yes/no
  !IF ( ALLOCATED(InData%BD_OtherState) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*2  ! OtherSt upper/lower bounds for each dimension
    !DO i2 = LBOUND(InData%OtherSt,2), UBOUND(InData%OtherSt,2)
    !DO i1 = LBOUND(InData%OtherSt,1), UBOUND(InData%OtherSt,1)
      Int_BufSz   = Int_BufSz + 3  ! OtherSt: size of buffers for each call to pack subtype
      CALL BD_PackOtherState( Re_Buf, Db_Buf, Int_Buf, InData%BD_OtherState, ErrStat2, ErrMsg2, .TRUE. ) ! OtherSt 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! OtherSt
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! OtherSt
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! OtherSt
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
    !END DO
  !END IF

  ! Parameters
  Int_BufSz   = Int_BufSz   + 1     ! p allocated yes/no
  !IF ( ALLOCATED(InData%BD_Parameter) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*1  ! p upper/lower bounds for each dimension
    !DO i1 = LBOUND(InData%p,1), UBOUND(InData%p,1)
      Int_BufSz   = Int_BufSz + 3  ! p: size of buffers for each call to pack subtype
      CALL BD_PackParam( Re_Buf, Db_Buf, Int_Buf, InData%BD_Parameter, ErrStat2, ErrMsg2, .TRUE. ) ! p 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! p
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! p
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! p
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
  !END IF

  ! Input Type
  Int_BufSz   = Int_BufSz   + 1     ! u allocated yes/no
  !IF ( ALLOCATED(InData%BD_Input) ) THEN
    Int_BufSz   = Int_BufSz   + 2*1  ! u upper/lower bounds for each dimension
    DO i1 = LBOUND(InData%BD_Input,1), UBOUND(InData%BD_Input,1)
      Int_BufSz   = Int_BufSz + 3  ! u: size of buffers for each call to pack subtype
      CALL BD_PackInput( Re_Buf, Db_Buf, Int_Buf, InData%BD_Input(i1), ErrStat2, ErrMsg2, .TRUE. ) ! u 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! u
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! u
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! u
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    END DO
  !END IF

  ! Output Type
  Int_BufSz   = Int_BufSz   + 1     ! y allocated yes/no
  !IF ( ALLOCATED(InData%BD_Output) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*1  ! y upper/lower bounds for each dimension
    !DO i1 = LBOUND(InData%y,1), UBOUND(InData%y,1)
      Int_BufSz   = Int_BufSz + 3  ! y: size of buffers for each call to pack subtype
      CALL BD_PackOutput( Re_Buf, Db_Buf, Int_Buf, InData%BD_Output, ErrStat2, ErrMsg2, .TRUE. ) ! y 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! y
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! y
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! y
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
  !END IF

  ! Misc vars
  Int_BufSz   = Int_BufSz   + 1     ! m allocated yes/no
  !IF ( ALLOCATED(InData%BD_MiscVar) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*1  ! m upper/lower bounds for each dimension
    !DO i1 = LBOUND(InData%m,1), UBOUND(InData%m,1)
      Int_BufSz   = Int_BufSz + 3  ! m: size of buffers for each call to pack subtype
      CALL BD_PackMisc( Re_Buf, Db_Buf, Int_Buf, InData%BD_MiscVar, ErrStat2, ErrMsg2, .TRUE. ) ! m 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! m
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! m
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! m
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    !END DO
  !END IF


  ! Skip Output and y_interp 

  ! Skip Input

  
  Int_BufSz   = Int_BufSz   + 1     ! InputTimes allocated yes/no
  IF ( ALLOCATED(InData%BD_InputTimes) ) THEN
    !Int_BufSz   = Int_BufSz   + 2*2  ! InputTimes upper/lower bounds for each dimension
    Int_BufSz   = Int_BufSz   + 2  ! InputTimes upper/lower bounds for each dimension
    Db_BufSz    = Db_BufSz   + SIZE(InData%BD_InputTimes)  ! InputTimes
  END IF
  
  
  ! Allocate the buffers

  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

  !IF ( .NOT. ALLOCATED(InData%BD_ContinuousState) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%x,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%x,1)
    !Int_Xferred = Int_Xferred + 2
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%x,2)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%x,2)
    !Int_Xferred = Int_Xferred + 2

    !DO i2 = LBOUND(InData%x,2), UBOUND(InData%x,2)
    !DO i1 = LBOUND(InData%x,1), UBOUND(InData%x,1)
      CALL BD_PackContState( Re_Buf, Db_Buf, Int_Buf, InData%BD_ContinuousState, ErrStat2, ErrMsg2, OnlySize ) ! x 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
    !END DO
  !END IF

  ! Discrete states
  !IF ( .NOT. ALLOCATED(InData%BD_DiscreteState) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%xd,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%xd,1)
    !Int_Xferred = Int_Xferred + 2
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%xd,2)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%xd,2)
    !Int_Xferred = Int_Xferred + 2
!
    !DO i2 = LBOUND(InData%xd,2), UBOUND(InData%xd,2)
    !DO i1 = LBOUND(InData%xd,1), UBOUND(InData%xd,1)
      CALL BD_PackDiscState( Re_Buf, Db_Buf, Int_Buf, InData%BD_DiscreteState, ErrStat2, ErrMsg2, OnlySize ) ! xd 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
    !END DO
  !END IF

  ! Constrain state
  !IF ( .NOT. ALLOCATED(InData%BD_ConstraintState) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%z,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%z,1)
    !Int_Xferred = Int_Xferred + 2
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%z,2)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%z,2)
    !Int_Xferred = Int_Xferred + 2

    !DO i2 = LBOUND(InData%z,2), UBOUND(InData%z,2)
    !DO i1 = LBOUND(InData%z,1), UBOUND(InData%z,1)
      CALL BD_PackConstrState( Re_Buf, Db_Buf, Int_Buf, InData%BD_ConstraintState, ErrStat2, ErrMsg2, OnlySize ) ! z 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
    !END DO
  !END IF

  ! Other states
  !IF ( .NOT. ALLOCATED(InData%BD_OtherState) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%OtherSt,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%OtherSt,1)
    !Int_Xferred = Int_Xferred + 2
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%OtherSt,2)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%OtherSt,2)
    !Int_Xferred = Int_Xferred + 2
!
    !DO i2 = LBOUND(InData%OtherSt,2), UBOUND(InData%OtherSt,2)
    !DO i1 = LBOUND(InData%OtherSt,1), UBOUND(InData%OtherSt,1)
      CALL BD_PackOtherState( Re_Buf, Db_Buf, Int_Buf, InData%BD_OtherState, ErrStat2, ErrMsg2, OnlySize ) ! OtherSt 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
    !END DO
  !END IF

  ! Parameters
  !IF ( .NOT. ALLOCATED(InData%BD_Parameter) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%p,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%p,1)
    !Int_Xferred = Int_Xferred + 2
!
    !DO i1 = LBOUND(InData%p,1), UBOUND(InData%p,1)
      CALL BD_PackParam( Re_Buf, Db_Buf, Int_Buf, InData%BD_Parameter, ErrStat2, ErrMsg2, OnlySize ) ! p 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
  !END IF

  ! Input
  !IF ( .NOT. ALLOCATED(InData%BD_Input) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%BD_Input,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%BD_Input,1)
    Int_Xferred = Int_Xferred + 2

    DO i1 = LBOUND(InData%BD_Input,1), UBOUND(InData%BD_Input,1)
      CALL BD_PackInput( Re_Buf, Db_Buf, Int_Buf, InData%BD_Input(i1), ErrStat2, ErrMsg2, OnlySize ) ! u 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    END DO
  !END IF

  ! Output
  !IF ( .NOT. ALLOCATED(InData%BD_Output) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%y,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%y,1)
    !Int_Xferred = Int_Xferred + 2
!
    !DO i1 = LBOUND(InData%y,1), UBOUND(InData%y,1)
      CALL BD_PackOutput( Re_Buf, Db_Buf, Int_Buf, InData%BD_Output, ErrStat2, ErrMsg2, OnlySize ) ! y 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
  !END IF

  ! Misc
  !IF ( .NOT. ALLOCATED(InData%BD_MiscVar) ) THEN
  !  IntKiBuf( Int_Xferred ) = 0
  !  Int_Xferred = Int_Xferred + 1
  !ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%m,1)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%m,1)
    !Int_Xferred = Int_Xferred + 2
!
    !DO i1 = LBOUND(InData%m,1), UBOUND(InData%m,1)
      CALL BD_PackMisc( Re_Buf, Db_Buf, Int_Buf, InData%BD_MiscVar, ErrStat2, ErrMsg2, OnlySize ) ! m 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    !END DO
  !END IF

  ! Skip output - y_interp - Input 

  ! Input times
  IF ( .NOT. ALLOCATED(InData%BD_InputTimes) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%BD_InputTimes,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%BD_InputTimes,1)
    Int_Xferred = Int_Xferred + 2
    !IntKiBuf( Int_Xferred    ) = LBOUND(InData%InputTimes,2)
    !IntKiBuf( Int_Xferred + 1) = UBOUND(InData%InputTimes,2)
    !Int_Xferred = Int_Xferred + 2

    DO i1 = LBOUND(InData%BD_InputTimes,1), UBOUND(InData%BD_InputTimes,1)
        DbKiBuf(Db_Xferred) = InData%BD_InputTimes(i1)
        Db_Xferred = Db_Xferred + 1
    END DO

  END IF
 END SUBROUTINE BD_PackBeamDyn_Data







! This routine unpacks the data and fills the Beamdyn user data for restart
! Inspired from FAST_UnPackBeamDyn_Data
 SUBROUTINE BD_UnPackBeamDyn_Data( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(BD_UsrDataType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'BD_UnPackBeamDyn_Data'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred   = 1
  Db_Xferred   = 1
  Int_Xferred  = 1

  ! Continous states
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! x not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !i2_l = IntKiBuf( Int_Xferred    )
    !i2_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%x)) DEALLOCATE(OutData%x)
    !ALLOCATE(OutData%x(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%x.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i2 = LBOUND(OutData%x,2), UBOUND(OutData%x,2)
    !DO i1 = LBOUND(OutData%x,1), UBOUND(OutData%x,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackContState( Re_Buf, Db_Buf, Int_Buf, OutData%BD_ContinuousState, ErrStat2, ErrMsg2 ) ! x 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO
    !END DO
  END IF

  ! Discrete states
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! xd not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !i2_l = IntKiBuf( Int_Xferred    )
    !i2_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%xd)) DEALLOCATE(OutData%xd)
    !ALLOCATE(OutData%xd(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%xd.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i2 = LBOUND(OutData%xd,2), UBOUND(OutData%xd,2)
    !DO i1 = LBOUND(OutData%xd,1), UBOUND(OutData%xd,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackDiscState( Re_Buf, Db_Buf, Int_Buf, OutData%BD_DiscreteState, ErrStat2, ErrMsg2 ) ! xd 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO
    !END DO
  END IF

  ! Constraint types
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! z not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !i2_l = IntKiBuf( Int_Xferred    )
    !i2_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%z)) DEALLOCATE(OutData%z)
    !ALLOCATE(OutData%z(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%z.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i2 = LBOUND(OutData%z,2), UBOUND(OutData%z,2)
    !DO i1 = LBOUND(OutData%z,1), UBOUND(OutData%z,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackConstrState( Re_Buf, Db_Buf, Int_Buf, OutData%BD_ConstraintState, ErrStat2, ErrMsg2 ) ! z 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO
    !END DO
  END IF

  ! Other states
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! OtherSt not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !i2_l = IntKiBuf( Int_Xferred    )
    !i2_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%OtherSt)) DEALLOCATE(OutData%OtherSt)
    !ALLOCATE(OutData%OtherSt(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%OtherSt.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i2 = LBOUND(OutData%OtherSt,2), UBOUND(OutData%OtherSt,2)
    !DO i1 = LBOUND(OutData%OtherSt,1), UBOUND(OutData%OtherSt,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackOtherState( Re_Buf, Db_Buf, Int_Buf, OutData%BD_OtherState, ErrStat2, ErrMsg2 ) ! OtherSt 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO
    !END DO
  END IF

  ! Paramters
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! p not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%p)) DEALLOCATE(OutData%p)
    !ALLOCATE(OutData%p(i1_l:i1_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%p.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i1 = LBOUND(OutData%p,1), UBOUND(OutData%p,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackParam( Re_Buf, Db_Buf, Int_Buf, OutData%BD_Parameter, ErrStat2, ErrMsg2 ) ! p 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO

  ! Input
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! u not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%u)) DEALLOCATE(OutData%u)
    !ALLOCATE(OutData%u(i1_l:i1_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%u.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    DO i1 = i1_l, i1_u
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackInput( Re_Buf, Db_Buf, Int_Buf, OutData%BD_Input(i1), ErrStat2, ErrMsg2 ) ! u 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    END DO
  END IF

  ! BD_Output
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! y not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%y)) DEALLOCATE(OutData%y)
    !ALLOCATE(OutData%y(i1_l:i1_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%y.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i1 = LBOUND(OutData%y,1), UBOUND(OutData%y,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackOutput( Re_Buf, Db_Buf, Int_Buf, OutData%BD_Output, ErrStat2, ErrMsg2 ) ! y 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO
  END IF

  ! Misc Vars
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! m not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    !i1_l = IntKiBuf( Int_Xferred    )
    !i1_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%m)) DEALLOCATE(OutData%m)
    !ALLOCATE(OutData%m(i1_l:i1_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%m.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !DO i1 = LBOUND(OutData%m,1), UBOUND(OutData%m,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL BD_UnpackMisc( Re_Buf, Db_Buf, Int_Buf, OutData%BD_MiscVar, ErrStat2, ErrMsg2 ) ! m 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    !END DO
  END IF

  ! Skip Output - y_interp - Input

  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! InputTimes not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    !i2_l = IntKiBuf( Int_Xferred    )
    !i2_u = IntKiBuf( Int_Xferred + 1)
    !Int_Xferred = Int_Xferred + 2
    !IF (ALLOCATED(OutData%InputTimes)) DEALLOCATE(OutData%InputTimes)
    !ALLOCATE(OutData%InputTimes(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    !IF (ErrStat2 /= 0) THEN 
    !   CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%InputTimes.', ErrStat, ErrMsg,RoutineName)
    !   RETURN
    !END IF
    !  DO i2 = LBOUND(OutData%InputTimes,2), UBOUND(OutData%InputTimes,2)
        DO i1 = LBOUND(OutData%BD_InputTimes,1), UBOUND(OutData%BD_InputTimes,1)
          OutData%BD_InputTimes(i1) = DbKiBuf(Db_Xferred)
          Db_Xferred = Db_Xferred + 1
        END DO
    !  END DO
  END IF

 END SUBROUTINE BD_UnPackBeamDyn_Data





 !----------------------------------------------------------------------------------------------------------------------------------
!> Routine that packs all of the data from one turbine instance into arrays and writes checkpoint files. If Unit is present and
!! greater than 0, it will append the data to an already open file. Otherwise, it opens a new file and writes header information
!! before writing the turbine data to the file.
!! Inspired from FAST_CreateCheckpoint_T
SUBROUTINE BD_CreateCheckpoint_T(t_initial, n_t_global, NumBeams, BD_Data, CheckpointRoot, ErrStat, ErrMsg, Unit )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   INTEGER(IntKi),           INTENT(IN   ) :: NumBeams         !< Number of turbines in this simulation
   TYPE(BD_UsrDataType),     INTENT(INOUT) :: BD_Data             !< all data for one instance of a turbine (INTENT(OUT) only because of hack for Bladed DLL)
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), OPTIONAL, INTENT(INOUT) :: Unit                !< unit number for output file

      ! local variables:
   REAL(ReKi),               ALLOCATABLE   :: ReKiBuf(:)
   REAL(DbKi),               ALLOCATABLE   :: DbKiBuf(:)
   INTEGER(IntKi),           ALLOCATABLE   :: IntKiBuf(:)

   INTEGER(B4Ki)                           :: ArraySizes(3)

   INTEGER(IntKi)                          :: iBeam
   INTEGER(IntKi)                          :: unOut               ! unit number for output file
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'BD_CreateCheckpoint_T'

   CHARACTER(1024)                         :: FileName            ! Name of the (output) checkpoint file

      ! init error status
   ErrStat = ErrID_None
   ErrMsg  = ""

   FileName    = TRIM(CheckpointRoot)//'.chkp'

   unOut=-1
   IF (PRESENT(Unit)) unOut = Unit

   IF ( unOut < 0 ) THEN

      CALL GetNewUnit( unOut, ErrStat2, ErrMsg2 )
      CALL OpenBOutFile ( unOut, FileName, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev ) then
            call cleanup()
            IF (.NOT. PRESENT(Unit)) THEN
               CLOSE(unOut)
               unOut = -1
            END IF

            RETURN
         end if

         ! checkpoint file header:
      WRITE (unOut, IOSTAT=ErrStat2)   INT(ReKi              ,B4Ki)     ! let's make sure we've got the correct number of bytes for reals on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   INT(DbKi              ,B4Ki)     ! let's make sure we've got the correct number of bytes for doubles on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   INT(IntKi             ,B4Ki)     ! let's make sure we've got the correct number of bytes for integers on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   AbortErrLev
      WRITE (unOut, IOSTAT=ErrStat2)   NumBeams                      ! Number of turbines
      WRITE (unOut, IOSTAT=ErrStat2)   t_initial                        ! initial time
      WRITE (unOut, IOSTAT=ErrStat2)   n_t_global                       ! current time step

   END IF

      ! Get the arrays of data to be stored in the output file
    CALL BD_PackBeamDyn_Data( ReKiBuf, DbKiBuf, IntKiBuf, BD_Data, ErrStat2, ErrMsg2 )
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        if (ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
        end if


   ArraySizes = 0
   IF ( ALLOCATED(ReKiBuf)  ) ArraySizes(1) = SIZE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) ArraySizes(2) = SIZE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) ArraySizes(3) = SIZE(IntKiBuf)


      ! data from current turbine at time step:
   WRITE (unOut, IOSTAT=ErrStat2)   ArraySizes                       ! Number of reals, doubles, and integers written to file
   WRITE (unOut, IOSTAT=ErrStat2)   ReKiBuf                          ! Packed reals
   WRITE (unOut, IOSTAT=ErrStat2)   DbKiBuf                          ! Packed doubles
   WRITE (unOut, IOSTAT=ErrStat2)   IntKiBuf                         ! Packed integers


   IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)

   WRITE (unOut, IOSTAT=ErrStat2)   BD_Data%theta_rot

   CLOSE(unOut)
   unOut = -1


   call cleanup()

contains
   subroutine cleanup()
      IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
      IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
      IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)
   end subroutine cleanup
END SUBROUTINE BD_CreateCheckpoint_T





!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is the inverse of FAST_CreateCheckpoint_T. It reads data from a checkpoint file and populates data structures for
!! the turbine instance.
!! Inspired from FAST_RestoreFromCheckpoint_T
SUBROUTINE BD_RestoreFromCheckpoint_T(t_initial, n_t_global, NumBeams, BD_Data, CheckpointRoot, ErrStat, ErrMsg, Unit )

   REAL(DbKi),               INTENT(INOUT) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(INOUT) :: n_t_global          !< loop counter
   INTEGER(IntKi),           INTENT(INOUT) :: NumBeams            !< Number of turbines in this simulation
   TYPE(BD_UsrDataType),     INTENT(INOUT) :: BD_Data          !< all data for one instance of a beam (bjj: note that is intent INOUT instead of OUT only because of a gfortran compiler memory issue)
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), OPTIONAL, INTENT(INOUT) :: Unit                !< unit number for output file

      ! local variables:
   REAL(ReKi),               ALLOCATABLE   :: ReKiBuf(:)
   REAL(DbKi),               ALLOCATABLE   :: DbKiBuf(:)
   INTEGER(IntKi),           ALLOCATABLE   :: IntKiBuf(:)

   INTEGER(B4Ki)                           :: ArraySizes(3)
   
   INTEGER(IntKi)                          :: iBeam
   INTEGER(IntKi)                          :: unIn                ! unit number for input file
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'BD_RestoreFromCheckpoint_T'

   CHARACTER(1024)                         :: FileName            ! Name of the (input) checkpoint file


   ErrStat=ErrID_None
   ErrMsg=""

   FileName    = TRIM(CheckpointRoot)//'.chkp'
   
   unIn=-1
   IF (PRESENT(Unit)) unIn = Unit

   IF ( unIn < 0 ) THEN

      CALL GetNewUnit( unIn, ErrStat2, ErrMsg2 )

      CALL OpenBInpFile ( unIn, FileName, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev ) RETURN

         ! checkpoint file header:
      READ (unIn, IOSTAT=ErrStat2)   ArraySizes     ! let's make sure we've got the correct number of bytes for reals, doubles, and integers on restart.

      IF ( ArraySizes(1) /= ReKi  ) CALL SetErrStat(ErrID_Fatal,"ReKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF ( ArraySizes(2) /= DbKi  ) CALL SetErrStat(ErrID_Fatal,"DbKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF ( ArraySizes(3) /= IntKi ) CALL SetErrStat(ErrID_Fatal,"IntKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(unIn)
         unIn = -1
         IF (PRESENT(Unit)) Unit = unIn
         RETURN
      END IF

      READ (unIn, IOSTAT=ErrStat2)   AbortErrLev
      READ (unIn, IOSTAT=ErrStat2)   NumBeams                      ! Number of turbines
      READ (unIn, IOSTAT=ErrStat2)   t_initial                        ! initial time
      READ (unIn, IOSTAT=ErrStat2)   n_t_global                       ! current time step

   END IF

      ! in case the Turbine data structure isn't empty on entry of this routine:
   !call FAST_DestroyTurbineType( Turbine, ErrStat2, ErrMsg2 )
    !DO iBeam=1,NumBeams
      ! data from current time step:
   READ (unIn, IOSTAT=ErrStat2)   ArraySizes                       ! Number of reals, doubles, and integers written to file

   ALLOCATE(ReKiBuf( ArraySizes(1)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate ReKiBuf", ErrStat, ErrMsg, RoutineName )
   ALLOCATE(DbKiBuf( ArraySizes(2)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate DbKiBuf", ErrStat, ErrMsg, RoutineName )
   ALLOCATE(IntKiBuf(ArraySizes(3)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate IntKiBuf", ErrStat, ErrMsg, RoutineName )

      ! Read the packed arrays
   IF (ErrStat < AbortErrLev) THEN

      READ (unIn, IOSTAT=ErrStat2)   ReKiBuf    ! Packed reals
         IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not read ReKiBuf", ErrStat, ErrMsg, RoutineName )
      READ (unIn, IOSTAT=ErrStat2)   DbKiBuf    ! Packed doubles
         IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not read DbKiBuf", ErrStat, ErrMsg, RoutineName )
      READ (unIn, IOSTAT=ErrStat2)   IntKiBuf   ! Packed integers
         IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not read IntKiBuf", ErrStat, ErrMsg, RoutineName )

   END IF

      ! Put the arrays back in the data types
   IF (ErrStat < AbortErrLev) THEN
      CALL BD_UnPackBeamDyn_Data( ReKiBuf, DbKiBuf, IntKiBuf, BD_Data, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

   READ (unIn, IOSTAT=ErrStat2)   BD_Data%theta_rot
   write(*,*) "Rotation angle = ",BD_Data%theta_rot

   !END DO


      ! close file if necessary (do this after unpacking turbine data, so that TurbID is set)
   !!!!IF (Turbine%TurbID == NumTurbines .OR. .NOT. PRESENT(Unit)) THEN
      CLOSE(unIn)
      unIn = -1
   !!!!END IF

   IF (PRESENT(Unit)) Unit = unIn


   IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)

   ! deal with files that were open:
   !IF (Turbine%p_FAST%WrTxtOutFile) THEN
   !   CALL OpenFunkFileAppend ( Turbine%y_FAST%UnOu, TRIM(Turbine%p_FAST%OutFileRoot)//'.out', ErrStat2, ErrMsg2)
   !   IF ( ErrStat2 >= AbortErrLev ) RETURN
   !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !   CALL WrFileNR ( Turbine%y_FAST%UnOu, '#Restarting here')
   !   WRITE(Turbine%y_FAST%UnOu, '()')
   !END IF
   ! (ignoring for now; will have fort.x files if any were open [though I printed a warning about not outputting binary files earlier])


END SUBROUTINE BD_RestoreFromCheckpoint_T


END MODULE BeamDynLib_Types
