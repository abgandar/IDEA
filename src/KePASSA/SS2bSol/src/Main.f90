Program Example

Use SS2b

Implicit None

Real :: Pos_Offset(3)
Real :: Vel_Offset(3)
Real :: RelTol
Real :: AbsTol
Real :: Final_State(6)
Real :: Final_Time
Real :: scaleX
Real :: scaleV

Character(len=30):: inputvar


    ! If the number of arguments is not correct, complain and stop
        If (COMMAND_ARGUMENT_COUNT() < 10 .OR. COMMAND_ARGUMENT_COUNT() > 11) Then
          Write(*,*) 'ERROR: Incorrect number of arguments'
          Stop
        End if


    ! Read the input arguments and store in the corresponding variables
    
        ! Initial Position offset along X axis [m]
        Call GET_COMMAND_ARGUMENT (1, inputvar)
        Read (inputvar,*) Pos_Offset(1)
        
        ! Initial Position offset along Y axis [m]
        Call GET_COMMAND_ARGUMENT (2, inputvar)
        Read (inputvar,*) Pos_Offset(2)
        
        ! Initial Position offset along Z axis [m]
        Call GET_COMMAND_ARGUMENT (3, inputvar)
        Read (inputvar,*) Pos_Offset(3)
        
        ! Initial Velocity offset along X axis [m/s]
        Call GET_COMMAND_ARGUMENT (4, inputvar)
        Read (inputvar,*) Vel_Offset(1)
        
        ! Initial Velocity offset along Y axis [m/s]
        Call GET_COMMAND_ARGUMENT (5, inputvar)
        Read (inputvar,*) Vel_Offset(2)

        ! Initial Velocity offset along Z axis [m/s]
        Call GET_COMMAND_ARGUMENT (6, inputvar)
        Read (inputvar,*) Vel_Offset(3)

        ! position scaling factor [to m/s]
        Call GET_COMMAND_ARGUMENT (7, inputvar)
        Read (inputvar,*) scaleX

        ! Initial Velocity offset along Z axis [m/s]
        Call GET_COMMAND_ARGUMENT (8, inputvar)
        Read (inputvar,*) scaleV

        ! Integration relative tolerance
        Call GET_COMMAND_ARGUMENT (9, inputvar)
        Read (inputvar,*) RelTol

        ! Integration absolute tolerance
        Call GET_COMMAND_ARGUMENT (10, inputvar)
        Read (inputvar,*) AbsTol

        ! Integration Timespan (optional argument)
        If (COMMAND_ARGUMENT_COUNT() == 11 ) Then
            Call GET_COMMAND_ARGUMENT (11, inputvar)
            Read (inputvar,*) Final_Time
        Else
            Final_Time = 288.12768941 * 3600. * 24.
        End If


    ! Run SS2b problem with prescribed parameters
    Call Problem_SS2b (scaleX*Pos_Offset, scaleV*Vel_Offset,  &
                       RelTol, AbsTol, Final_State, Final_Time = Final_Time, Verbosity = .FALSE.)


    ! Print solution on console
    Write(*,*) Final_State(1:6)
    


End Program
