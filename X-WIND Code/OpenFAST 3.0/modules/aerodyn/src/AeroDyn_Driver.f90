!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
program AeroDyn_Driver

   use AeroDyn_Driver_Subs   
    
   implicit none   
   
      ! Program variables

   real(DbKi)                                     :: time                 !< Variable for storing time, in seconds 
   real(DbKi)                                     :: dT_Dvr               !< copy of DT, to make sure AD didn't change it
                                                    
   type(Dvr_SimData)                              :: DvrData              ! The data required for running the AD driver
   type(AeroDyn_Data)                             :: AD                   ! AeroDyn data 
                                                  
   integer(IntKi)                                 :: iCase                ! loop counter (for driver case)
   integer(IntKi)                                 :: nt                   ! loop counter (for time step)
   integer(IntKi)                                 :: j                    ! loop counter (for array of inputs)
   integer(IntKi)                                 :: errStat              ! Status of error message
   character(ErrMsgLen)                           :: errMsg               ! Error message if ErrStat /= ErrID_None

   !integer                                        :: StrtTime (8)                            ! Start time of simulation (including intialization)
   !integer                                        :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
   !real(ReKi)                                     :: PrevClockTime                           ! Clock time at start of simulation in seconds
   !real                                           :: UsrTime1                                ! User CPU time for simulation initialization
   !real                                           :: UsrTime2                                ! User CPU time for simulation (without intialization)
   !real                                           :: UsrTimeDiff                             ! Difference in CPU time from start to finish of program execution
   !real(DbKi)                                     :: TiLstPrn                                ! The simulation time of the last print
   !real(DbKi)                                     :: SttsTime                                ! Amount of time between screen status messages (sec)
   !integer                                        :: n_SttsTime                              ! Number of time steps between screen status messages (-)
   logical                                        :: AD_Initialized

   real(ReKi)                                      :: RotAzimuth                             ! Rotor Azimuth (aligned with blade 1)
   !real(ReKi)                                      :: TeetAng                                ! Teeter angle
   !real(ReKi)                                      :: TeetAngVel                             ! Teeter angular velocity

                            

   errStat     = ErrID_None
   errMsg      = ''
   AD_Initialized = .false.
   
   time        = 0.0 ! seconds
      
            
      ! Get the current time
   !call date_and_time ( Values=StrtTime )                               ! Let's time the whole simulation
   !call cpu_time ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   
   
      ! initialize this driver:
   call Dvr_Init( DvrData, ErrStat, ErrMsg)
      call CheckError()
   
   
   do iCase = 1, DvrData%NumCases
   
!      call WrScr( NewLine//'Running case '//trim(num2lstr(iCase))//' of '//trim(num2lstr(DvrData%NumCases))//'.' )
   
         ! Set the Initialization input data for AeroDyn based on the Driver input file data, and initialize AD
         ! (this also initializes inputs to AD for first time step)
      dT_Dvr   = DvrData%Cases(iCase)%dT
      call Init_AeroDyn(iCase, DvrData, AD, dT_Dvr, errStat, errMsg)
         call CheckError()
         AD_Initialized = .true.
         
         if (.not. EqualRealNos( dT_Dvr, DvrData%Cases(iCase)%dT ) ) then
            ErrStat = ErrID_Fatal
            ErrMsg = 'AeroDyn changed the time step for case '//trim(num2lstr(iCase))//'. Change DTAero to "default".'
            call CheckError()
         end if
                                    
      if (iCase.eq.1) then
         call Dvr_InitializeOutputFile(DvrData%numBlades, iCase, DvrData%Cases(iCase), DvrData%OutFileData, errStat, errMsg)
            call CheckError()
      endif
      
      RotAzimuth = 0.0_ReKi
      do nt = 1, DvrData%Cases(iCase)%numSteps
         
         !...............................
         ! set AD inputs for nt (and keep values at nt-1 as well)
         !...............................
         
         call Set_AD_Inputs(iCase,nt,RotAzimuth,DvrData,AD,errStat,errMsg) ! u(1) is at nt+1, u(2) is at nt
            call CheckError()
   
         time = AD%InputTime(2)
            
            ! Calculate outputs at nt - 1

         call AD_CalcOutput( time, AD%u(2), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat, errMsg )
            call CheckError()
  



         call Dvr_WriteOutputLine(DvrData%OutFileData, nt, RotAzimuth, AD%y%rotors(1)%WriteOutput, DvrData%Cases(iCase), iCase, errStat, errMsg)
            call CheckError()


            
            
            ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt

         call AD_UpdateStates( time, nt-1, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat, errMsg )
            call CheckError()
      
                  print *, nt,"step complete----------------------"
      end do !nt=1,numSteps
      
   end do !iCase = 1, DvrData%NumCases
   
   call Dvr_End()
   
contains
!................................   
   subroutine CheckError()
   
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            call Dvr_End()
         end if
      end if
         
   end subroutine CheckError
!................................   
   subroutine Dvr_End()
   
         ! Local variables
      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      
      character(*), parameter                       :: RoutineName = 'Dvr_End'
         ! Close the output file
      if (DvrData%OutFileData%unOutFile > 0) close(DvrData%OutFileData%unOutFile)
            
      if ( AD_Initialized ) then
         call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end if
           
      call AD_Dvr_DestroyDvr_SimData( DvrData, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      call AD_Dvr_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
      
      
   end subroutine Dvr_End
!................................   
end program AeroDyn_Driver
   
