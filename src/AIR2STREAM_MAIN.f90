PROGRAM air2stream

USE commondata
IMPLICIT NONE
INTEGER :: i,j,k,status
REAL(KIND=8):: T1,T2

WRITE(*,*) '       .__       ________            __                                  '
WRITE(*,*) '_____  |__|______\_____  \   _______/  |________   ____ _____    _____   '
WRITE(*,*) '\__  \ |  \_  __ \/  ____/  /  ___/\   __\_  __ \_/ __ \\__  \  /     \  '
WRITE(*,*) ' / __ \|  ||  | \/       \  \___ \  |  |  |  | \/\  ___/ / __ \|  Y Y  \ '
WRITE(*,*) '(____  /__||__|  \_______ \/____  > |__|  |__|    \___  >____  /__|_|  / '
WRITE(*,*) '     \/                  \/     \/                    \/     \/      \/  '
WRITE(*,*) 'Version 1.0.0 - October 2015'

!-------------------------------------------------------------------------------------
!
! Provided by Sebastiano Piccolroaz and Marco Toffolon
!
! Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
! email contacts: s.piccolroaz@unitn.it, marco.toffolon@unitn.it
!
! How to cite:
! Toffolon M. and Piccolroaz S. (2015), A hybrid model for river water temperature 
! as a function of air temperature and discharge, Environmental Research Letters  
! (under review)
!
! Related works (air2water):
! Piccolroaz S., M. Toffolon, and B. Majone (2013), A simple lumped model to convert
! air temperature into surface water temperature in lakes, Hydrol. Earth Syst. Sci., 
! 17, 3323-3338, doi:10.5194/hess-17-3323-2013
!
! Toffolon M., S. Piccolroaz, B. Majone, A.M. Soja, F. Peeters, M. Schmid and A. Wüest 
! (2014), Prediction of surface water temperature from air temperature in lakes with 
! different morphology, Limnology and Oceanography, 59(6), 2185-2202, doi: 10.4319/lo.2014.59.6.2185
!
! Toffolon M., S. Piccolroaz, and B. Majone (2015), The role of stratification on lakes' thermal 
! response: The case of Lake Superior, Water Resources Research, doi: 10.1002/2014WR016555, in press
!-------------------------------------------------------------------------------------

CALL CPU_TIME(T1)

! allocation of parameter matrices
ALLOCATE(parmin(n_par),stat=status) 
ALLOCATE(parmax(n_par),stat=status)
ALLOCATE(flag_par(n_par),stat=status) 
ALLOCATE(par(n_par),stat=status) 
ALLOCATE(par_best(n_par),stat=status) 

! read input data
CALL read_calibration

! aggregate calibration data on the basis of time resolution
CALL aggregation

! calculae mean, TSS and standard deviation
CALL statis
WRITE(*,*) 'mean, TSS and standard deviation (calibration)'
WRITE(*,*)  SNGL(mean_obs),SNGL(TSS_obs),SNGL(std_obs)

OPEN(UNIT=11,FILE=TRIM(folder)//'/1_'//TRIM(runmode)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')

IF (runmode .eq. 'FORWARD') THEN
    CALL forward_mode
ELSE IF (runmode .eq. 'PSO') THEN
    CALL PSO_mode
ELSE IF (runmode .eq. 'LATHYP') THEN
    CALL LH_mode
END IF

! forward with the best set of parameters   
CALL forward

CALL CPU_TIME(T2)
print *, 'Computation time was ', T2-T1, 'seconds.'

STOP
END PROGRAM air2stream
