Module SS2b

Use RADAU_Driver

Implicit None

Private
Public Problem_SS2b



Contains


  Subroutine Problem_SS2b (Pos_Offset, Vel_Offset, RelTol, AbsTol, Final_State, Final_Time, Verbosity)

    Real, intent(in)  ::  Pos_Offset(3)     ! Initial Cartesian position offset vector [m]
    Real, intent(in)  ::  Vel_Offset(3)     ! Initial Cartesian velocity offset vector [m/s]
    Real, intent(in)  ::  RelTol            ! Integration relative tolerance
    Real, intent(in)  ::  AbsTol            ! Integration absolute tolerance
    Real, intent(out) ::  Final_State(6)    ! Final Cartesian state vector [m, m/s]
    
    Real,    optional, intent(in)  ::  Final_Time        ! Final Time when to stop integration
    Logical, optional, intent(in)  ::  Verbosity         ! Indicates whether to display or not information on terminal


  ! Declaration of variables
    
    Type (RADAU)      :: Integrator      ! Object containing the RADAU integrator
    Type (RADAU_Sol)  :: Solution        ! Object containing the integration solution/output

    Real              :: Y(6)            ! Integrator input/output: Initial/final state-vector
    Real              :: T0              ! Initial Time
    Real              :: Tf              ! Final Time
    Real              :: Error           ! Final integration error


  ! Set Integration parameters

      Integrator % RelTol   = RelTol          ! Relative Allowed Error or Tolerance
      Integrator % AbsTol   = AbsTol          ! Absolute Allowed Error or Tolerance
      Integrator % MaxSteps = 0               ! Maximum number of allowed integration steps. (unlimited if 0)

      Integrator % Hessenberg = .TRUE.        ! Transform Jacobian to Hessenberg form (doesn't apply to banded jacobians nor implicit systems)
      Integrator % Gustafsson = .TRUE.        ! Use Gustafsson's predictive stepsize control

      Integrator % MinStages = 1              ! Minimal number of stages. Possible values: 1, 3, 5, 7 (yielding orders 1, 5, 9, 13, respectively)
      Integrator % MaxStages = 7              ! Maximal number of stages. Possible values: 1, 3, 5, 7 (yielding orders 1, 5, 9, 13, respectively)
      
      If (Present(Verbosity) .AND. Verbosity) Then
        Integrator % Verbose     = Verbosity  ! To display real-time integration progress on console
        Integrator % Statistics  = Verbosity  ! To display integration statistics on console
        Integrator % RefreshTime = 1.         ! Time interval (in seconds) to refresh integr. progress (if verbose)
      Else
        Integrator % Verbose     = .FALSE.    ! To display real-time integration progress on console
        Integrator % Statistics  = .FALSE.    ! To display integration statistics on console
      End If



  ! Set Time Span for simulation
      T0 = 0.

      If (Present(Final_Time)) Then
        Tf = Final_Time
      Else
        Tf = 288.12768941 * 3600. * 24.
      End If


  ! Set initial conditions for simulation

      Y(1) =     0.0         +  Pos_Offset(1)
      Y(2) = -5888.9727e3    +  Pos_Offset(2)
      Y(3) = -3400.0e3       +  Pos_Offset(3)
      Y(4) =    10.691338e3  +  Vel_Offset(1)
      Y(5) =     0.0         +  Vel_Offset(2)
      Y(6) =     0.0         +  Vel_Offset(3)

      If (Present(Verbosity) .AND. Verbosity) Then
        Write(*,*)
        Write(*,*) "Initial State Vector"
        Write(*,*) "X: ", Y(1), "m"
        Write(*,*) "Y: ", Y(2), "m"
        Write(*,*) "Z: ", Y(3), "m"
        Write(*,*) "Vx:", Y(4), "m/s"
        Write(*,*) "Vy:", Y(5), "m/s"
        Write(*,*) "Vz:", Y(6), "m/s"
      End If


  ! Perform integration of Stiefel & Scheifele's 2b Problem with the RADAU integrator
      Call Integrator % Integrate ( RHS_SS2b, T0, [T0, Tf], Y, Solution )


  ! Display results and compare the error with respect to the Stiefel & Scheifele's exact solution

      If (Present(Verbosity) .AND. Verbosity) Then
        Write(*,*)
        Write(*,*) "Final State Vector"
        Write(*,*) "X: ", Y(1), "m"
        Write(*,*) "Y: ", Y(2), "m"
        Write(*,*) "Z: ", Y(3), "m"
        Write(*,*) "Vx:", Y(4), "m/s"
        Write(*,*) "Vy:", Y(5), "m/s"
        Write(*,*) "Vz:", Y(6), "m/s"
        Write(*,*)
      End If


  ! Return output: final state vector
      Final_State = Y


  End Subroutine




  Subroutine RHS_SS2b ( N, t, y, yp, rpar, ipar )
      
      Integer, intent(in)  ::  N
      Real,    intent(in)  ::  t
      Real,    intent(in)  ::  y(N)
      Real,    intent(out) ::  yp(N)
      Real,    intent(in)  ::  rpar(*)
      Integer, intent(in)  ::  ipar(*)

      Real :: Rho, Rho2, Rho3, z2
      Real :: mu, R_p, J2, cost_J2, mu_m, d_m, w_m
      Real :: aX_J2, aY_J2, aZ_J2, aX_m, aY_m, aZ_m
      Real :: XM, YM, ZM, XMrel, YMrel, ZMrel, Rho_M


      ! Derived Quantities
      Rho   =  norm2( y(1:3) )
      Rho2  =  Rho**2
      Rho3  =  Rho**3
      z2    =  y(3)**2


      ! Parameters
      mu    =  3.98601e14
      R_p   =  6371220.0
      J2    =  0.00108265
      mu_m  =  4.90266e12
      d_m   =  384400000.
      w_m   =  2.665315780887e-6
      

      ! J2 Perturbation
      cost_J2 = - 1.5 * mu * J2 * R_p**2 / Rho**7   

      aX_J2 = cost_J2 * y(1) * (  Rho2 - 5.*z2)
      aY_J2 = cost_J2 * y(2) * (  Rho2 - 5.*z2)
      aZ_J2 = cost_J2 * y(3) * (3.*Rho2 - 5.*z2)
      

      ! Moon Perturbation           
      XM =   d_m * sin(w_m*t)
      YM = - d_m * cos(w_m*t) / 2. * sqrt(3.)
      ZM = - d_m * cos(w_m*t) / 2.

      XMrel = y(1) - XM;
      YMrel = y(2) - YM;
      ZMrel = y(3) - ZM;
      Rho_M = norm2( [XMrel, YMrel, ZMrel] )
         
      aX_m = - mu_m * (XMrel/Rho_M**3 + XM/d_m**3);
      aY_m = - mu_m * (YMrel/Rho_M**3 + YM/d_m**3);
      aZ_m = - mu_m * (ZMrel/Rho_M**3 + ZM/d_m**3);
      

      ! RHS of ODE system
      yp(1) =   y(4)
      yp(2) =   y(5)
      yp(3) =   y(6)
      yp(4) = - y(1) / Rho3 * ( mu  ) + (aX_J2 + aX_m)
      yp(5) = - y(2) / Rho3 * ( mu  ) + (aY_J2 + aY_m)
      yp(6) = - y(3) / Rho3 * ( mu  ) + (aZ_J2 + aZ_m)

  End Subroutine


End Module
