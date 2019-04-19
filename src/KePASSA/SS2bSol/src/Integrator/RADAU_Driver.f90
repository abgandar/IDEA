Module RADAU_Driver

Use RADAU_Core

Implicit None

Private

Public:: RADAU
Public:: RADAU_Sol
Public:: RADAU_Mon

Public:: RADAU_RHS
Public:: RADAU_JAC
Public:: RADAU_MAS




Type RADAU

    Real    :: RelTol = 1e-13          ! Relative integration tolerance
    Real    :: AbsTol = 1e-7           ! Absolute integration tolerance
    Real    :: SafetyFactor = 0.9      ! Safety factor                                    (if 0, default value is used)
    Real    :: Hratio_Lower = 0.2      ! Lower bound to hnew/h ratio                      (if 0, default value is used))
    Real    :: Hratio_Upper = 8.       ! Upper bound to hnew/h ratio                      (if 0, default value is used))
    Real    :: Hratio_unchanged = 1.2  ! Don't update stepsize while hnew/h < this        (if 0, default value is used))
    Real    :: CF_Lower = 0.002        ! Contractivity Factor threshold to increase order (if 0, default value is used))
    Real    :: CF_Upper = 0.8          ! Contractivity Factor threshold to decrease order (if 0, default value is used))
    Real    :: Jac_Recomp = 0.001      ! Decides whether Jacobian should be recomputed    (if 0, default value is used; if negative, recomputes every accepted step)
    Real    :: Hmax  = 0.              ! Maximum stepsize                                 (if 0, the whole integration interval is used)
    Real    :: Hinit = 0.              ! Initial stepsize                                 (if 0, a suitable initial stepsize will be guessed)
    Logical :: FixedStepsize = .FALSE. ! If TRUE, fixed-stepsize integration is performed with stepsize 'Hinit'.
                                       !  Careful: with large stepsizes in stiff problems, Newton iterations may not converge and issues may arise yielding NaN; decreasing the stepsize should fix this

    Logical :: Hessenberg = .TRUE.     ! Transform Jacobian to Hessenberg form. Doesn't apply to banded jacobians nor implicit systems
    Logical :: Gustafsson = .TRUE.     ! Use Gustafsson's predictive stepsize control

    Integer :: MinStages = 3           ! Minimal number of stages. Possible values: 1, 3, 5, 7 (yielding orders 1, 5, 9, 13, respectively)
    Integer :: MaxStages = 7           ! Maximal number of stages. Possible values: 1, 3, 5, 7 (yielding orders 1, 5, 9, 13, respectively)

    Integer :: MaxSteps = 1000000      ! Max. number of allowed integr. steps (if 0, default value is used)

    Integer :: MLJAC = -1              ! If a Banded Jacobian structure is assumed, this represents the number of non-zero diagonals below the main diagonal. If negative, a Full Jacobian is assumed.
    Integer :: MUJAC = -1              ! If a Banded Jacobian structure is assumed, this represents the number of non-zero diagonals above the main diagonal. If negative, a Full Jacobian is assumed.

    Integer :: MLMAS = -1              ! If a Banded Mass-Matrix structure is assumed, this represents the number of non-zero diagonals below the main diagonal. If negative, a Full matrix is assumed.
    Integer :: MUMAS = -1              ! If a Banded Mass-Matrix structure is assumed, this represents the number of non-zero diagonals above the main diagonal. If negative, a Full matrix is assumed.

    Integer :: M1 = 0                  ! Take advantage of the ODE structure where: Y(i)' = Y(i+M2)  for  i=1,...,M1.
    Integer :: M2 = 0                  ! Very useful for second-order systems: P'=V, V'=G(P,V) where P and V are vectors of dim. N/2, one has to put M1=M2=N/2.

    Logical :: Verbose     = .FALSE.   ! If TRUE displays real-time integration progress
    Logical :: Statistics  = .FALSE.   ! If TRUE displays integration statistics
    Integer :: Console     = 6         ! Console unit (6 is default terminal output)
    Real    :: RefreshTime = 2.        ! Time interval (in seconds) to refresh integr. progress (if verbose)

    Real,    allocatable :: rpar(:)    ! Real    parameters array passed to ODE RHS
    Integer, allocatable :: ipar(:)    ! Integer parameters array passed to ODE RHS

    Integer :: Fcalls                  ! Integr. Stats: Total number of ODE RHS  evaluations (exluding calls for Jacobian calculation)
    Integer :: Jcalls                  ! Integr. Stats: Total number of Jacobian evaluations (either analytically or numerically)
    Integer :: GoodSteps               ! Integr. Stats: Total number of successful integr. steps
    Integer :: BadSteps                ! Integr. Stats: Total number of failed integr. steps
    Integer :: LUdecomp                ! Integr. Stats: Total number of LU decompositions of the matrices
    Integer :: Subst                   ! Integr. Stats: Number of Forward-Backward substitutions, of all systems that have to be solved per simplified Newton iter. (stepsize selection excluded)
    Real    :: RunTime                 ! Integr. Stats: Total integration runtime

Contains

    Procedure, private, pass :: Integrate_Straight
    Procedure, private, pass :: Integrate_Dense

    Generic :: Integrate => Integrate_Straight, Integrate_Dense

End Type RADAU




Type RADAU_Sol

    Real, Allocatable :: t(:)     ! Value of independent variable: ti(1:N)
    Real, Allocatable :: y(:,:)   ! Values of dependent variables: y(1:Dim, 1:N)

    Integer :: Dim
    Integer :: N

Contains

    Procedure, pass :: Initialize   =>  RADAU_Sol_Initialize
    Procedure, pass :: Finalize     =>  RADAU_Sol_Finalize
    Procedure, pass :: Append       =>  RADAU_Sol_Append
    Procedure, pass :: Export       =>  RADAU_Sol_Export

End Type RADAU_Sol




Type RADAU_Mon

    Real,    Allocatable :: t(:)          ! Value of independent variable: ti(1:N)
    Real,    Allocatable :: Stepsize(:)   ! Value of integr. stepsize taken  to reach the current integration step (assumed 0 for initial step)
    Integer, Allocatable :: Order(:)      ! Value of the integr. order  used to reach the current integration step (assumed 0 for initial step)
    Integer, Allocatable :: Steps(:)      ! Value of previous steps     used to reach the current integration step (assumed 0 for initial step)
    Integer, Allocatable :: Stages(:)     ! Value of Runge-Kutta stages used to reach the current integration step (assumed 0 for initial step)
    Integer, Allocatable :: NewtonIter(:) ! Number of Newton Iterations needed at     the current integration step (assumed 0 for initial step)

    Integer, Allocatable :: Fcalls(:)     ! Cumulative number of ODE RHS  evaluations (exluding calls for Jacobian calculation)
    Integer, Allocatable :: Jcalls(:)     ! Cumulative number of Jacobian evaluations (either analytically or numerically)
    Integer, Allocatable :: GoodSteps(:)  ! Cumulative number of successful (i.e. accepted) integration steps
    Integer, Allocatable :: BadSteps(:)   ! Cumulative number of failed (i.e. rejected) integration steps
    Integer, Allocatable :: LUdecomp(:)   ! Cumulative number of LU decompositions of the matrices
    Integer, Allocatable :: Subst(:)      ! Cumulative number of Forward-Backward substitutions, of all systems that have to be solved per simplified Newton iter. (stepsize selection excluded)

    Integer :: N

Contains

    Procedure, pass :: Initialize   =>  RADAU_Mon_Initialize
    Procedure, pass :: Finalize     =>  RADAU_Mon_Finalize
    Procedure, pass :: Append       =>  RADAU_Mon_Append
    Procedure, pass :: Export       =>  RADAU_Mon_Export

End Type RADAU_Mon




Abstract Interface

    Subroutine RADAU_RHS ( N, t, y, yp, rpar, ipar )
        Integer, intent(in)  ::  N
        Real,    intent(in)  ::  t
        Real,    intent(in)  ::  y(N)
        Real,    intent(out) ::  yp(N)
        Real,    intent(in)  ::  rpar(*)
        Integer, intent(in)  ::  ipar(*)
    End Subroutine

    Subroutine RADAU_JAC ( N, t, y, dfy, LDFY, rpar, ipar )
        Integer, intent(in)  ::  N
        Real,    intent(in)  ::  t
        Real,    intent(in)  ::  y(N)
        Real,    intent(out) ::  dfy(LDFY,N)
        Integer, intent(in)  ::  LDFY
        Real,    intent(in)  ::  rpar(*)
        Integer, intent(in)  ::  ipar(*)
    End Subroutine

    Subroutine RADAU_MAS (N, AM, LMAS, RPAR, IPAR)
        Integer, intent(in)  ::  N
        Real,    intent(out) ::  AM(LMAS,N)
        Integer, intent(in)  ::  LMAS
        Real,    intent(in)  ::  rpar(*)
        Integer, intent(in)  ::  ipar(*)
    End Subroutine

End Interface



Real,             Save, Pointer ::  Grid(:)
Real,             Save          ::  RefreshTime
Logical,          Save          ::  Verbose
Integer,          Save          ::  Console = 6
Type (RADAU_Sol), Save, Pointer ::  Solution_ptr
Type (RADAU_Mon), Save, Pointer ::  Monitor_ptr






Contains





!/////////////////////////////////////////////////////////////////
!//                     Integrator Drivers                      //
!/////////////////////////////////////////////////////////////////

Subroutine Integrate_Straight (This, ODE, T0, Tf, Y0, Solution, Monitor, JAC, MAS)

    Class (RADAU), intent(inout) :: This    ! · Instances of class 'RADAU'

    Procedure (RADAU_RHS)   :: ODE          ! · Procedure of the form ODE(t, y(:), F(:)), where F(:) contains the derivatives of the ODE system
    
    Real,    intent(inout)  :: T0           ! · As input is the initial time, at which the simulation begins
                                            !   As output, returns the time at which the integration stops, be it Tf or sooner if an event-driven stop happened
    Real,    intent(in)     :: Tf           ! · Final time, at which the simulation stops
    Real,    intent(inout)  :: Y0(:)        ! · As input is a vector containing the initial conditions at T0, i.e. Y(t=T0).
                                            !   As output, returns the solution at Tf, i.e. Y(t=Tf)
    
    Class (RADAU_Sol), target,           intent(inout) :: Solution  ! · Object containing the solution computed at every time-step of the integration
    Class (RADAU_Mon), target, optional, intent(out)   :: Monitor   ! · Object containing the integration process monitoring information at every integr. step

    Procedure (RADAU_JAC), optional :: JAC  ! · Procedure providing the Jacobian dF(t,y)/dy. If none is provided, the Jacobian is computed numerically.
    Procedure (RADAU_MAS), optional :: MAS  ! · Procedure providing the mass matrix M such that M*y' = F(t,y). If none is provided the identity matrix is assumed.


    ! Declare other variables
    Real,    allocatable :: RelTol(:)
    Real,    allocatable :: AbsTol(:)
    Real,    allocatable :: WORK(:)
    Integer, allocatable :: IWORK(:)

    Real(8) :: Time1, Time2

    Integer :: NEQ, NSMAX, LJAC, LMAS, LE
    Integer :: IJAC, MLJAC, MUJAC
    Integer :: IMAS, MLMAS, MUMAS
    Integer :: LWORK, LIWORK, IDID

    Procedure(RADAU_JAC), pointer :: Jacobian
    Procedure(RADAU_MAS), pointer :: MassMatrix

    

    ! Initialization and setup of the integrator
    NEQ    =  size (Y0)
    NSMAX  =  This % MaxStages
    
    If (Present(JAC)) Then
        IJAC = 1
        Jacobian => JAC
    Else
        IJAC = 0
    End If

    If (Present(MAS)) Then
        IMAS = 1
        MassMatrix => MAS
    Else
        IMAS = 0
    End If
    
    MLJAC = NEQ - This % M1
    MUJAC = NEQ - This % M1
    MLMAS = NEQ - This % M1
    MUMAS = NEQ - This % M1

    If (This%MLJAC >= 0 .AND. This%MLJAC < NEQ - This % M1)   MLJAC = This % MLJAC
    If (This%MUJAC >= 0 .AND. This%MUJAC < NEQ - This % M1)   MUJAC = This % MUJAC
    If (This%MLMAS >= 0 .AND. This%MLMAS < NEQ - This % M1)   MLMAS = This % MLMAS
    If (This%MUMAS >= 0 .AND. This%MUMAS < NEQ - This % M1)   MUMAS = This % MUMAS
    
    If (MLJAC == NEQ) Then
        LJAC  =  NEQ                 ! Full Jacobian
        LE    =  NEQ
    Else
        LJAC  =    MLJAC + MUJAC + 1 ! Banded Jacobian
        LE    =  2*MLJAC + MUJAC + 1
    End If
    
    If (IMAS == 0) Then
        LMAS  =  0
    Else
        If (MLMAS == NEQ) Then
            LMAS  =  NEQ               ! Full Mass-Matrix
        Else
            LMAS  =  MLMAS + MUMAS + 1 ! Banded Mass-Matrix
        End If
    End If
    
    LWORK  =  NEQ * (LJAC + LMAS + NSMAX*LE + 3*NSMAX + 3) + 20
    LIWORK =  (2 + (NSMAX-1) / 2) * NEQ + 20
    

    Allocate (RelTol (NEQ) )  ;  RelTol = This % RelTol
    Allocate (AbsTol (NEQ) )  ;  AbsTol = This % AbsTol
    Allocate (WORK  (LWORK))  ;  WORK   = 0.
    Allocate (IWORK (LIWORK)) ;  IWORK  = 0

    WORK(1)  = epsilon(0.)
    WORK(2)  = This % SafetyFactor
    WORK(3)  = This % Jac_Recomp
    WORK(6)  = This % Hratio_unchanged
    WORK(7)  = This % Hmax
    WORK(8)  = This % Hratio_Lower
    WORK(9)  = This % Hratio_Upper
    WORK(10) = This % CF_Lower
    WORK(11) = This % CF_Upper
    
    If (      This % Hessenberg) IWORK(1) = 1
    If (.NOT. This % Gustafsson) IWORK(8) = 2
    
    IWORK(5)  = NEQ
    IWORK(2)  = This % MaxSteps
    IWORK(9)  = This % M1
    IWORK(10) = This % M2
    IWORK(11) = This % MinStages
    IWORK(12) = This % MaxStages

    Verbose     =  This % Verbose
    Console     =  This % Console
    RefreshTime =  This % RefreshTime


    ! Initialize Solution
    Call Solution%Initialize (T0, Y0)
    Solution_ptr => Solution


    ! Initialize Monitor
    If (Present(Monitor)) Then
        Call Monitor%Initialize (T0)
        Monitor_ptr => Monitor
    Else
        Monitor_ptr => Null()
    End If


    ! Print Message on console
    If (This%Verbose) Write(This%Console, '(/,A25,/)') '> Integration Running...'
    If (This%Verbose) Write(This%Console, '(I15, F10.2, A3)') 0, 0., "%"


    ! Begin timer
    Call cpu_time (Time1)


    ! Use 'Grid' to remember initial and final values
    Allocate (Grid, Source = [T0, Tf])


    ! Perform integration
    Call RADAU_Integrator (NEQ, ODE, T0, Y0, Tf,                            &
                           This%Hinit, RelTol, AbsTol, 1,                   &
                           Jacobian,   IJAC, MLJAC, MUJAC,                  &
                           MassMatrix ,IMAS, MLMAS, MUMAS,                  &
                           Straight_Output, 1, WORK, LWORK, IWORK, LIWORK,  &
                           This%RPAR, This%IPAR, IDID,                      &
                           This%Console, This%FixedStepsize)


    ! Finalize Integration
    Call cpu_time (Time2)
    This % RunTime = real(Time2 - Time1)


    ! Update statistical counters
    This % Fcalls    = IWORK(14)
    This % Jcalls    = IWORK(15)
    This % GoodSteps = IWORK(17)
    This % BadSteps  = IWORK(18)
    This % LUdecomp  = IWORK(19)
    This % Subst     = IWORK(20)


    ! Finalize/Consolidate Solution
    Call Solution%Finalize ()


    ! Finalize/Consolidate Monitoring information
    If (Present(Monitor))  Call Monitor%Finalize ()


    ! Deallocate pointers
    Deallocate (RelTol)
    Deallocate (AbsTol)

    Nullify (Grid)


    ! Print Message on console
    If (This%Verbose .AND. IDID >= 0) Then
        Write(This%Console, '(I15, F10.2, A3, /)') This % GoodSteps, 100., "%"
    End If


    ! Display Integration Summary
    If (This % Statistics) Then
        Write(This%Console, '(A18,F10.3,A9)')       "Elapsed Time:",      This % RunTime, "seconds"
        Write(This%Console, '(A18,I10)')            "Function Calls:",    This % Fcalls
        Write(This%Console, '(A18,I10)')            "Jacobian Calls:",    This % Jcalls
        Write(This%Console, '(A18,I10)')            "LU Decompositions:", This % LUdecomp
        Write(This%Console, '(A18,I10)')            "Forw-Backw. Subst:", This % Subst
        Write(This%Console, '(A18,I10,F8.3,A2)')    "Correct Steps:",     This % GoodSteps, &
                         100. * This%GoodSteps / Max(1, (This%GoodSteps + This%BadSteps) ), "%"
        Write(This%Console, '(A18,I10,F8.3,A2)')    "Bad Steps:",         This % BadSteps,  &
                         100. * This%BadSteps  / Max(1, (This%GoodSteps + This%BadSteps) ), "%"
        Write(This%Console, '(A18,I10,/)')          "Solution Points:",   Solution % N
    End If


    ! Notify if any integration error occurred
    Select Case (IDID)
    Case (-1)
      Write(This%Console, *) "<RADAU ERROR>  INPUT IS NOT CONSISTENT"
    Case (-2)
      Write(This%Console, *) "<RADAU ERROR>  LARGER NMAX IS NEEDED"
    Case (-3)
      Write(This%Console, *) "<RADAU ERROR>  STEP SIZE BECOMES TOO SMALL"
    Case (-4)
      Write(This%Console, *) "<RADAU ERROR>  MATRIX IS REPEATEDLY SINGULAR"
    End Select
    

End Subroutine Integrate_Straight





Subroutine Integrate_Dense (This, ODE, T0, T, Y0, Solution, Monitor, JAC, MAS)

    Class (RADAU),  intent(inout) :: This  ! · Instances of class 'DOP853'

    Procedure (RADAU_RHS)       :: ODE     ! · Procedure of the form ODE(t, y(:), F(:)), where F(:) contains the derivatives of the ODE system
    
    Real,         intent(inout) :: T0      ! · As input is the initial time, at which the simulation begins
                                           !   As output, returns the time at which the integration stops, be it Tf or sooner if an event-driven stop happened
    Real, target, intent(in)    :: T(:)    ! · Time at which the simulation provides dense output
    Real,         intent(inout) :: Y0(:)   ! · As input is a vector containing the initial conditions at T0, i.e. Y(t=T0).
                                           !   As output, returns the solution at Tf, i.e. Y(t=Tf)
    
    Class (RADAU_Sol), target,           intent(inout) :: Solution  ! · Object containing the solution computed at every time-step of the integration
    Class (RADAU_Mon), target, optional, intent(out)   :: Monitor   ! · Object containing the integration process monitoring information at every integr. step

    Procedure (RADAU_JAC), optional :: JAC ! · Procedure providing the Jacobian dF(t,y)/dy. If none is provided, the Jacobian is computed numerically.
    Procedure (RADAU_MAS), optional :: MAS ! · Procedure providing the mass matrix M such that M*y' = F(t,y). If none is provided the identity matrix is assumed.


    ! Declare other variables
    Real,    allocatable :: RelTol(:)
    Real,    allocatable :: AbsTol(:)
    Real,    allocatable :: WORK(:)
    Integer, allocatable :: IWORK(:)

    Real(8) :: Time1, Time2

    Integer :: NEQ, NSMAX, LJAC, LMAS, LE
    Integer :: IJAC, MLJAC, MUJAC
    Integer :: IMAS, MLMAS, MUMAS
    Integer :: LWORK, LIWORK, IDID

    Procedure(RADAU_JAC), pointer :: Jacobian
    Procedure(RADAU_MAS), pointer :: MassMatrix

    

    ! Initialization and setup of the integrator
    NEQ    =  size (Y0)
    NSMAX  =  This % MaxStages
    
    If (Present(JAC)) Then
        IJAC = 1
        Jacobian => JAC
    Else
        IJAC = 0
    End If

    If (Present(MAS)) Then
        IMAS = 1
        MassMatrix => MAS
    Else
        IMAS = 0
    End If

    MLJAC = NEQ - This % M1
    MUJAC = NEQ - This % M1
    MLMAS = NEQ - This % M1
    MUMAS = NEQ - This % M1

    If (This%MLJAC >= 0 .AND. This%MLJAC < NEQ - This % M1)   MLJAC = This % MLJAC
    If (This%MUJAC >= 0 .AND. This%MUJAC < NEQ - This % M1)   MUJAC = This % MUJAC
    If (This%MLMAS >= 0 .AND. This%MLMAS < NEQ - This % M1)   MLMAS = This % MLMAS
    If (This%MUMAS >= 0 .AND. This%MUMAS < NEQ - This % M1)   MUMAS = This % MUMAS
    
    If (MLJAC == NEQ) Then
        LJAC  =  NEQ                 ! Full Jacobian
        LE    =  NEQ
    Else
        LJAC  =    MLJAC + MUJAC + 1 ! Banded Jacobian
        LE    =  2*MLJAC + MUJAC + 1
    End If
    
    If (IMAS == 0) Then
        LMAS  =  0
    Else
        If (MLMAS == NEQ) Then
            LMAS  =  NEQ               ! Full Mass-Matrix
        Else
            LMAS  =  MLMAS + MUMAS + 1 ! Banded Mass-Matrix
        End If
    End If
    
    LWORK  =  NEQ * (LJAC + LMAS + NSMAX*LE + 3*NSMAX + 3) + 20
    LIWORK =  (2 + (NSMAX-1) / 2) * NEQ + 20
    

    Allocate (RelTol (NEQ) )  ;  RelTol = This % RelTol
    Allocate (AbsTol (NEQ) )  ;  AbsTol = This % AbsTol
    Allocate (WORK  (LWORK))  ;  WORK   = 0.
    Allocate (IWORK (LIWORK)) ;  IWORK  = 0

    WORK(1)  = epsilon(0.)
    WORK(2)  = This % SafetyFactor
    WORK(3)  = This % Jac_Recomp
    WORK(6)  = This % Hratio_unchanged
    WORK(7)  = This % Hmax
    WORK(8)  = This % Hratio_Lower
    WORK(9)  = This % Hratio_Upper
    WORK(10) = This % CF_Lower
    WORK(11) = This % CF_Upper
    
    If (      This % Hessenberg) IWORK(1) = 1
    If (.NOT. This % Gustafsson) IWORK(8) = 2
    
    IWORK(5)  = NEQ
    IWORK(2)  = This % MaxSteps
    IWORK(9)  = This % M1
    IWORK(10) = This % M2
    IWORK(11) = This % MinStages
    IWORK(12) = This % MaxStages

    Verbose     =  This % Verbose
    Console     =  This % Console
    RefreshTime =  This % RefreshTime


    ! Initialize Solution
    Solution % N   = Size(T)
    Solution % Dim = NEQ

    If ( Allocated (Solution % t) )  Deallocate (Solution % t)
    If ( Allocated (Solution % y) )  Deallocate (Solution % y)
    
    Allocate ( Solution % t (SIZE(T)) )
    Allocate ( Solution % y (NEQ,SIZE(T)) )
    
    Solution % t (  1) = T0
    Solution % y (:,1) = Y0

    Solution_ptr => Solution


    ! Initialize Monitor
    If (Present(Monitor)) Then
        Call Monitor%Initialize (T0)
        Monitor_ptr => Monitor
    Else
        Monitor_ptr => Null()
    End If


    ! Print Message on console
    If (This%Verbose) Write(This%Console, '(/,A25,/)') '> Integration Running...'
    If (This%Verbose) Write(This%Console, '(I15, F10.2, A3)') 0, 0., "%"


    ! Begin timer
    Call cpu_time (Time1)


    ! Remember grid points
    Grid => T


    ! Perform integration
    Call RADAU_Integrator (NEQ, ODE, T0, Y0, T(Size(T)),                 &
                           This%Hinit, RelTol, AbsTol, 1,                &
                           Jacobian,   IJAC, MLJAC, MUJAC,               &
                           MassMatrix ,IMAS, MLMAS, MUMAS,               &
                           Dense_Output, 1, WORK, LWORK, IWORK, LIWORK,  &
                           This%RPAR, This%IPAR, IDID,                   &
                           This%Console, This%FixedStepsize)


    ! Finalize Integration
    Call cpu_time (Time2)
    This % RunTime = real(Time2 - Time1)


    ! Update statistical counters
    This % Fcalls    = IWORK(14)
    This % Jcalls    = IWORK(15)
    This % GoodSteps = IWORK(17)
    This % BadSteps  = IWORK(18)
    This % LUdecomp  = IWORK(19)
    This % Subst     = IWORK(20)


    ! Finalize/Consolidate Solution
    Solution % t (  Size(T)) = T0
    Solution % y (:,Size(T)) = Y0


    ! Finalize/Consolidate Monitoring information
    If (Present(Monitor))  Call Monitor%Finalize ()


    ! Deallocate pointers
    Deallocate (RelTol)
    Deallocate (AbsTol)
    
    Nullify (Grid)


    ! Print Message on console
    If (This%Verbose .AND. IDID >= 0) Then
        Write(This%Console, '(I15, F10.2, A3, /)') This % GoodSteps, 100., "%"
    End If


    ! Display Integration Summary
    If (This % Statistics) Then
        Write(This%Console, '(A18,F10.3,A9)')       "Elapsed Time:",      This % RunTime, "seconds"
        Write(This%Console, '(A18,I10)')            "Function Calls:",    This % Fcalls
        Write(This%Console, '(A18,I10)')            "Jacobian Calls:",    This % Jcalls
        Write(This%Console, '(A18,I10)')            "LU Decompositions:", This % LUdecomp
        Write(This%Console, '(A18,I10)')            "Forw-Backw. Subst:", This % Subst
        Write(This%Console, '(A18,I10,F8.3,A2)')    "Correct Steps:",     This % GoodSteps, &
                         100. * This%GoodSteps / Max(1, (This%GoodSteps + This%BadSteps) ), "%"
        Write(This%Console, '(A18,I10,F8.3,A2)')    "Bad Steps:",         This % BadSteps,  &
                         100. * This%BadSteps  / Max(1, (This%GoodSteps + This%BadSteps) ), "%"
        Write(This%Console, '(A18,I10,/)')          "Solution Points:",   Solution % N
    End If


    ! Notify if any integration error occurred
    Select Case (IDID)
    Case (-1)
      Write(This%Console, *) "<RADAU ERROR>  INPUT IS NOT CONSISTENT"
    Case (-2)
      Write(This%Console, *) "<RADAU ERROR>  LARGER NMAX IS NEEDED"
    Case (-3)
      Write(This%Console, *) "<RADAU ERROR>  STEP SIZE BECOMES TOO SMALL"
    Case (-4)
      Write(This%Console, *) "<RADAU ERROR>  MATRIX IS REPEATEDLY SINGULAR"
    End Select
    

End Subroutine Integrate_Dense





Subroutine Straight_Output (NR, XOLD, X, Y, CONT, LRC, N, RPAR, IPAR, IRTRN, CONF, STATS)
! Get integration output at successful steps

        Integer :: NR, LRC, N, IPAR(:), IRTRN
        Integer :: CONF(4), STATS(6)
        Real    :: XOLD, X, Y(N), CONT(LRC), RPAR(:)
        
        Real(8), Save :: Time1 = 0d0
        Real(8), Save :: Time2 = 0d0


        If (NR == 1) Return


        ! Append state vector to solution
        Call Solution_ptr%Append (X, Y)


        ! Append current integration stats to Monitor
        If (Associated(Monitor_ptr))  Call Monitor_ptr % Append (X, CONF, STATS)


        ! Display integration progress
        Call cpu_time(Time2)

        If ( Verbose .AND. Time2-Time1 > RefreshTime ) Then
            Write(Console, '(I15, F10.2, A3)') NR-1,  &
             (X - Grid(1)) / (Grid(size(Grid)) - Grid(1))*100., "%"
            Time1 = Time2
        End If


End Subroutine Straight_Output





Subroutine Dense_Output (NR, XOLD, X, Y, CONT, LRC, N, RPAR, IPAR, IRTRN, CONF, STATS)
! Get integration output at requested values of independent integration variable
! by using 'dense output' capabilitites of the RADAU integrator, made available
! by the continuous collocation solution (Subroutine "CONTRA")

        Integer :: NR, LRC, N, IPAR(:), IRTRN
        Integer :: CONF(4), STATS(6), i
        Real    :: XOLD, X, Y(N), CONT(LRC), RPAR(:), XOUT, YOUT(N)
        
        Integer, Save :: Index

        Real(8), Save :: Time1 = 0d0
        Real(8), Save :: Time2 = 0d0
        

        If (NR == 1) Index = 2
        
        XOUT = Grid(Index)


        ! Append state vector to solution
        Do While ( X > XOUT )
            
            Do i = 1, N
              YOUT(i) = CONTRA (i, XOUT, CONT, LRC)
            End Do
            
            Solution_ptr % t (  Index) = XOUT
            Solution_ptr % y (:,Index) = YOUT

            Index = Index + 1
            XOUT  = Grid(Index)

        End Do


        ! Append current integration stats to Monitor
        If (Associated(Monitor_ptr) .AND. NR > 1)  Call Monitor_ptr % Append (X, CONF, STATS)


        ! Display integration progress
        Call cpu_time(Time2)
        
        If ( Verbose .AND. Time2-Time1 > RefreshTime ) Then
            Write(Console, '(I15, F10.2, A3)') NR-1,  &
             (X - Grid(1)) / (Grid(size(Grid)) - Grid(1))*100., "%"
            Time1 = Time2
        End If

 
End Subroutine Dense_Output





Subroutine RADAU_Sol_Initialize (This, t, y)

    Class (RADAU_Sol), intent(inout) ::  This
    Real,              intent(in)    ::  t
    Real,              intent(in)    ::  y(:)


    ! Initialize parameters
    This % N   = 1
    This % Dim = Size(y)

    ! Allocate variables to keep the solution information
    If ( Allocated (This % t) )  Deallocate (This % t)
    If ( Allocated (This % y) )  Deallocate (This % y)
    
    Allocate ( This % t (1) )
    Allocate ( This % y (This%Dim,1) )
    
    This % t (1)   = t
    This % y (:,1) = y


End Subroutine RADAU_Sol_Initialize





Subroutine RADAU_Sol_Append (This, t, y)

    Class (RADAU_Sol), intent(inout) ::  This
    Real,              intent(in)    ::  t
    Real,              intent(in)    ::  y(:)

    Real, Allocatable :: Temp_t(:)
    Real, Allocatable :: Temp_y(:,:)


    ! Check if there is room for new data
    If ( This%N == Size(This%t) ) Then

        ! Allocate temp. variables and save backup
        Allocate ( Temp_t, Source = This % t )
        Allocate ( Temp_y, Source = This % y )

        ! Reallocate
        Deallocate (This % t);    Allocate ( This % t (2*This%N) )
        Deallocate (This % y);    Allocate ( This % y (This%Dim, 2*This%N) )

        ! Transfer backup data
        This % t (1:This%N)    = Temp_t
        This % y (:, 1:This%N) = Temp_y

        ! Free up memory
        Deallocate (Temp_t)
        Deallocate (Temp_y)

    End If


    ! Update value of current index, N
    This%N = This%N + 1

    ! Append new data
    This % t (This%N)    = t
    This % y (:, This%N) = y


End Subroutine RADAU_Sol_Append





Subroutine RADAU_Sol_Finalize (This)

    Class (RADAU_Sol), intent(inout) ::  This
    
    Real, Allocatable :: Temp_t(:)
    Real, Allocatable :: Temp_y(:,:)


    ! Allocate temp. variables and save backup
    Allocate ( Temp_t, Source = This % t (1:This%N) )
    Allocate ( Temp_y, Source = This % y (:, 1:This%N) )

    ! Reallocate and transfer backup data
    Deallocate (This % t);    Allocate ( This % t, Source = Temp_t )
    Deallocate (This % y);    Allocate ( This % y, Source = Temp_y )

    ! Free up memory
    Deallocate (Temp_t)
    Deallocate (Temp_y)


End Subroutine RADAU_Sol_Finalize





Subroutine RADAU_Sol_Export (This, FileName)

    Class (RADAU_Sol), intent(in) :: This
    Character (len=*), intent(in) :: FileName

    Logical :: Unit_OK
    Logical :: Unit_Opened
    Integer :: u, i

    Character(len=10) :: Columns


    ! Check if Integration Solution is available and can be exported
    If ( .NOT. Allocated (This % t) ) Then
        Write(Console, *) "<SOLUTION ERROR>  INTEGRATION SOLUTION HAS NOT BEEN COMPUTED YET AND THUS CANNOT BE EXPORTED TO A FILE"
        Return
    End If


    ! Choose a suitable Unit number to open file
    Do u = 10, 100
        Inquire (unit = u, exist = Unit_OK, opened = Unit_Opened)
        If (Unit_OK .AND. .NOT. Unit_Opened) Exit
    End Do


    ! Open File
    Open (Unit = u, File = FileName)


    ! Write Header
    Write(u,'(A17)') 'ZONE t="Solution"'
    Write(u,*)


    ! Convet ODE system dimension into character
    Write(Columns, '(I10)') This % Dim + 1


    ! Write Trajectory
    Do i = 1, This%N
        Write (u, '(' // trim(adjustl(Columns)) // 'ES25.16E3)') This % t(i), This % y(:,i)
    End Do


    ! Close File
    Close(u)


End Subroutine RADAU_Sol_Export




Subroutine RADAU_Mon_Initialize (This, t)

    Class (RADAU_Mon), intent(inout) ::  This
    Real,              intent(in)    ::  t


    ! Initialize parameters
    This % N   = 1

    ! Allocate variables to keep the integration information
    If ( Allocated (This % t         ) )  Deallocate (This % t)
    If ( Allocated (This % Stepsize  ) )  Deallocate (This % Stepsize)
    If ( Allocated (This % Order     ) )  Deallocate (This % Order)
    If ( Allocated (This % Steps     ) )  Deallocate (This % Steps)
    If ( Allocated (This % Stages    ) )  Deallocate (This % Stages)
    If ( Allocated (This % NewtonIter) )  Deallocate (This % NewtonIter)
    If ( Allocated (This % Fcalls    ) )  Deallocate (This % Fcalls)
    If ( Allocated (This % Jcalls    ) )  Deallocate (This % Jcalls)
    If ( Allocated (This % GoodSteps ) )  Deallocate (This % GoodSteps)
    If ( Allocated (This % BadSteps  ) )  Deallocate (This % BadSteps)
    If ( Allocated (This % LUdecomp  ) )  Deallocate (This % LUdecomp)
    If ( Allocated (This % Subst     ) )  Deallocate (This % Subst)
    
    Allocate ( This % t          (1) )
    Allocate ( This % Stepsize   (1) )
    Allocate ( This % Order      (1) )
    Allocate ( This % Steps      (1) )
    Allocate ( This % Stages     (1) )
    Allocate ( This % NewtonIter (1) )
    Allocate ( This % Fcalls     (1) )
    Allocate ( This % Jcalls     (1) )
    Allocate ( This % GoodSteps  (1) )
    Allocate ( This % BadSteps   (1) )
    Allocate ( This % LUdecomp   (1) )
    Allocate ( This % Subst      (1) )
    
    This % t          (1)  =  t
    This % Stepsize   (1)  =  0.
    This % Order      (1)  =  0
    This % Steps      (1)  =  0
    This % Stages     (1)  =  0
    This % NewtonIter (1)  =  0
    This % Fcalls     (1)  =  0
    This % Jcalls     (1)  =  0
    This % GoodSteps  (1)  =  0
    This % BadSteps   (1)  =  0
    This % LUdecomp   (1)  =  0
    This % Subst      (1)  =  0


End Subroutine RADAU_Mon_Initialize





Subroutine RADAU_Mon_Append (This, t, CONF, STATS)

    Class (RADAU_Mon), intent(inout) ::  This
    Real,              intent(in)    ::  t
    Integer,           intent(in)    ::  CONF (4)
    Integer,           intent(in)    ::  STATS(6)

    Integer :: NN

    Real,    Allocatable :: Temp_t(:)
    Real,    Allocatable :: Temp_Stepsize(:)
    Integer, Allocatable :: Temp_Order(:)
    Integer, Allocatable :: Temp_Steps(:)
    Integer, Allocatable :: Temp_Stages(:)
    Integer, Allocatable :: Temp_NewtonIter(:)
    Integer, Allocatable :: Temp_Fcalls(:)
    Integer, Allocatable :: Temp_Jcalls(:)
    Integer, Allocatable :: Temp_GoodSteps(:)
    Integer, Allocatable :: Temp_BadSteps(:)
    Integer, Allocatable :: Temp_LUdecomp(:)
    Integer, Allocatable :: Temp_Subst(:)


    ! Check if there is room for new data
    If ( This%N == Size(This%t) ) Then

        ! Allocate temp. variables and save backup
        Allocate ( Temp_t,          Source = This % t          )
        Allocate ( Temp_Stepsize,   Source = This % Stepsize   )
        Allocate ( Temp_Order,      Source = This % Order      )
        Allocate ( Temp_Steps,      Source = This % Steps      )
        Allocate ( Temp_Stages,     Source = This % Stages     )
        Allocate ( Temp_NewtonIter, Source = This % NewtonIter )
        Allocate ( Temp_Fcalls,     Source = This % Fcalls     )
        Allocate ( Temp_Jcalls,     Source = This % Jcalls     )
        Allocate ( Temp_GoodSteps,  Source = This % GoodSteps  )
        Allocate ( Temp_BadSteps,   Source = This % BadSteps   )
        Allocate ( Temp_LUdecomp,   Source = This % LUdecomp   )
        Allocate ( Temp_Subst,      Source = This % Subst      )

        ! Reallocate
        NN = 2 * This%N
        Deallocate (This % t);            Allocate ( This % t          (NN) )
        Deallocate (This % Stepsize);     Allocate ( This % Stepsize   (NN) )
        Deallocate (This % Order);        Allocate ( This % Order      (NN) )
        Deallocate (This % Steps);        Allocate ( This % Steps      (NN) )
        Deallocate (This % Stages);       Allocate ( This % Stages     (NN) )
        Deallocate (This % NewtonIter);   Allocate ( This % NewtonIter (NN) )
        Deallocate (This % Fcalls);       Allocate ( This % Fcalls     (NN) )
        Deallocate (This % Jcalls);       Allocate ( This % Jcalls     (NN) )
        Deallocate (This % GoodSteps);    Allocate ( This % GoodSteps  (NN) )
        Deallocate (This % BadSteps);     Allocate ( This % BadSteps   (NN) )
        Deallocate (This % LUdecomp);     Allocate ( This % LUdecomp   (NN) )
        Deallocate (This % Subst);        Allocate ( This % Subst      (NN) )

        ! Transfer backup data
        This % t          (1:This%N)  =  Temp_t
        This % Stepsize   (1:This%N)  =  Temp_Stepsize
        This % Order      (1:This%N)  =  Temp_Order
        This % Steps      (1:This%N)  =  Temp_Steps
        This % Stages     (1:This%N)  =  Temp_Stages
        This % NewtonIter (1:This%N)  =  Temp_NewtonIter
        This % Fcalls     (1:This%N)  =  Temp_Fcalls
        This % Jcalls     (1:This%N)  =  Temp_Jcalls
        This % GoodSteps  (1:This%N)  =  Temp_GoodSteps
        This % BadSteps   (1:This%N)  =  Temp_BadSteps
        This % LUdecomp   (1:This%N)  =  Temp_LUdecomp
        This % Subst      (1:This%N)  =  Temp_Subst

        ! Free up memory
        Deallocate (Temp_t         )
        Deallocate (Temp_Stepsize  )
        Deallocate (Temp_Order     )
        Deallocate (Temp_Steps     )
        Deallocate (Temp_Stages    )
        Deallocate (Temp_NewtonIter)
        Deallocate (Temp_Fcalls    )
        Deallocate (Temp_Jcalls    )
        Deallocate (Temp_GoodSteps )
        Deallocate (Temp_BadSteps  )
        Deallocate (Temp_LUdecomp  )
        Deallocate (Temp_Subst     )

    End If


    ! Update value of current index, N
    This%N = This%N + 1

    ! Append new data
    This % t          (This%N)  =  t
    This % Stepsize   (This%N)  =  abs( t - This % t (This%N-1) )
    This % Order      (This%N)  =  CONF (1)
    This % Steps      (This%N)  =  CONF (2)
    This % Stages     (This%N)  =  CONF (3)
    This % NewtonIter (This%N)  =  CONF (4)
    This % Fcalls     (This%N)  =  STATS(1)
    This % Jcalls     (This%N)  =  STATS(2)
    This % GoodSteps  (This%N)  =  STATS(3)
    This % BadSteps   (This%N)  =  STATS(4)
    This % LUdecomp   (This%N)  =  STATS(5)
    This % Subst      (This%N)  =  STATS(6)


End Subroutine RADAU_Mon_Append





Subroutine RADAU_Mon_Finalize (This)

    Class (RADAU_Mon), intent(inout) ::  This
    
    Real,    Allocatable :: Temp_t(:)
    Real,    Allocatable :: Temp_Stepsize(:)
    Integer, Allocatable :: Temp_Order(:)
    Integer, Allocatable :: Temp_Steps(:)
    Integer, Allocatable :: Temp_Stages(:)
    Integer, Allocatable :: Temp_NewtonIter(:)
    Integer, Allocatable :: Temp_Fcalls(:)
    Integer, Allocatable :: Temp_Jcalls(:)
    Integer, Allocatable :: Temp_GoodSteps(:)
    Integer, Allocatable :: Temp_BadSteps(:)
    Integer, Allocatable :: Temp_LUdecomp(:)
    Integer, Allocatable :: Temp_Subst(:)


    ! Allocate temp. variables and save backup
    Allocate ( Temp_t,          Source = This % t          (1:This%N) )
    Allocate ( Temp_Stepsize,   Source = This % Stepsize   (1:This%N) )
    Allocate ( Temp_Order,      Source = This % Order      (1:This%N) )
    Allocate ( Temp_Steps,      Source = This % Steps      (1:This%N) )
    Allocate ( Temp_Stages,     Source = This % Stages     (1:This%N) )
    Allocate ( Temp_NewtonIter, Source = This % NewtonIter (1:This%N) )
    Allocate ( Temp_Fcalls,     Source = This % Fcalls     (1:This%N) )
    Allocate ( Temp_Jcalls,     Source = This % Jcalls     (1:This%N) )
    Allocate ( Temp_GoodSteps,  Source = This % GoodSteps  (1:This%N) )
    Allocate ( Temp_BadSteps,   Source = This % BadSteps   (1:This%N) )
    Allocate ( Temp_LUdecomp,   Source = This % LUdecomp   (1:This%N) )
    Allocate ( Temp_Subst,      Source = This % Subst      (1:This%N) )

    ! Reallocate and transfer backup data
    Deallocate (This % t);            Allocate ( This % t,          Source = Temp_t         )
    Deallocate (This % Stepsize);     Allocate ( This % Stepsize,   Source = Temp_Stepsize  )
    Deallocate (This % Order);        Allocate ( This % Order,      Source = Temp_Order     )
    Deallocate (This % Steps);        Allocate ( This % Steps,      Source = Temp_Steps     )
    Deallocate (This % Stages);       Allocate ( This % Stages,     Source = Temp_Stages    )
    Deallocate (This % NewtonIter);   Allocate ( This % NewtonIter, Source = Temp_NewtonIter)
    Deallocate (This % Fcalls);       Allocate ( This % Fcalls,     Source = Temp_Fcalls    )
    Deallocate (This % Jcalls);       Allocate ( This % Jcalls,     Source = Temp_Jcalls    )
    Deallocate (This % GoodSteps);    Allocate ( This % GoodSteps,  Source = Temp_GoodSteps )
    Deallocate (This % BadSteps);     Allocate ( This % BadSteps,   Source = Temp_BadSteps  )
    Deallocate (This % LUdecomp);     Allocate ( This % LUdecomp,   Source = Temp_LUdecomp  )
    Deallocate (This % Subst);        Allocate ( This % Subst,      Source = Temp_Subst     )

    ! Free up memory
    Deallocate (Temp_t         )
    Deallocate (Temp_Stepsize  )
    Deallocate (Temp_Order     )
    Deallocate (Temp_Steps     )
    Deallocate (Temp_Stages    )
    Deallocate (Temp_NewtonIter)
    Deallocate (Temp_Fcalls    )
    Deallocate (Temp_Jcalls    )
    Deallocate (Temp_GoodSteps )
    Deallocate (Temp_BadSteps  )
    Deallocate (Temp_LUdecomp  )
    Deallocate (Temp_Subst     )


End Subroutine RADAU_Mon_Finalize





Subroutine RADAU_Mon_Export (This, FileName)

    Class (RADAU_Mon), intent(in) :: This
    Character (len=*), intent(in) :: FileName

    Logical :: Unit_OK
    Logical :: Unit_Opened
    Integer :: u, i

    Character(len=10) :: Columns


    ! Check if Monitoring information is available and can be exported
    If ( .NOT. Allocated (This % t) ) Then
        Write(Console, *) "<MONITORING ERROR>  MONITORING INFO IS NOT AVAILABLE AND THUS CANNOT BE EXPORTED TO A FILE"
        Return
    End If


    ! Choose a suitable Unit number to open file
    Do u = 10, 100
        Inquire (unit = u, exist = Unit_OK, opened = Unit_Opened)
        If (Unit_OK .AND. .NOT. Unit_Opened) Exit
    End Do


    ! Open File
    Open (Unit = u, File = FileName)


    ! Write Header
    Write(u,'(A135)') 'VARIABLES = '   //  &
                      '"t", '          //  &
                      '"Stepsize", '   //  &
                      '"Order", '      //  &
                      '"Steps", '      //  &
                      '"Stages", '     //  &
                      '"NewtonIter", ' //  &
                      '"Fcalls", '     //  &
                      '"Jcalls", '     //  &
                      '"GoodSteps", '  //  &
                      '"BadSteps", '   //  &
                      '"LUdecomp", '   //  &
                      '"Subst"'
    
    Write(u,'(A16)') 'ZONE t="Monitor"'
    Write(u,*)


    ! Write Trajectory
    Do i = 1, This%N
        Write (u, '(2ES25.16E3, 4I4, 6I13)')  &
                    This % t(i),              &
                    This % Stepsize(i),       &
                    This % Order(i),          &
                    This % Steps(i),          &
                    This % Stages(i),         &
                    This % NewtonIter(i),     &
                    This % Fcalls(i),         &
                    This % Jcalls(i),         &
                    This % GoodSteps(i),      &
                    This % BadSteps(i),       &
                    This % LUdecomp(i),       &
                    This % Subst(i)
    End Do


    ! Close File
    Close(u)


End Subroutine RADAU_Mon_Export




End Module RADAU_Driver