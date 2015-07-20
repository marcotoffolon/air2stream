PROGRAM PSO   ! Particle Swarm Optimization

USE commondata
IMPLICIT NONE
INTEGER :: i,j,k,status
INTEGER :: count
REAL(KIND=8):: eff_index
REAL(KIND=8):: T1,T2
REAL(KIND=8):: norm
CHARACTER(LEN=100) :: fmt
! PSO parameters
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)::x,v,pbest
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::fit,r,gbest,fitbest
REAL(KIND=8)::w,dw,dxmax,dvmax,foptim


WRITE(*,*) '       .__       ________            __                                  '
WRITE(*,*) '_____  |__|______\_____  \   _______/  |________   ____ _____    _____   '
WRITE(*,*) '\__  \ |  \_  __ \/  ____/  /  ___/\   __\_  __ \_/ __ \\__  \  /     \  '
WRITE(*,*) ' / __ \|  ||  | \/       \  \___ \  |  |  |  | \/\  ___/ / __ \|  Y Y  \ '
WRITE(*,*) '(____  /__||__|  \_______ \/____  > |__|  |__|    \___  >____  /__|_|  / '
WRITE(*,*) '     \/                  \/     \/                    \/     \/      \/  '
WRITE(*,*) 'Version 1.0 - July 2015'

!-------------------------------------------------------------------------------------
!
! Provided by Marco Toffolon and Sebastiano Piccolroaz
!
! Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
! email contacts: marco.toffolon@unitn.it, s.piccolroaz@unitn.it
!
! How to cite:
! Toffolon M. and Piccolroaz S. (2015), A hybrid model for river water temperature 
! as a function of air temperature and discharge., Environmental Research Letters  
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
fmt='(8(f10.6,1x),f11.8)'

! aggregate calibration data on the basis of time resolution
CALL aggregation

! call of subroutine to calculae mean, TSS and standard deviation
CALL statis
WRITE(*,*) 'mean, TSS and standard deviation (calibration)'
WRITE(*,*)  SNGL(mean_obs),SNGL(TSS_obs),SNGL(std_obs)

WRITE(*,*) 'N. particles = ',n_particles,', N. run = ', n_run

ALLOCATE(gbest(n_par))

IF (runmode=='RANSAM') THEN
    foptim=-9999
   
    DO i=1,n_run                    ! number of iterations
        CALL random_seed()
        DO j=1,n_par
            CALL random_number(rRANSAM)
            par(j)=parmin(j) +(parmax(j) - parmin(j))*rRANSAM
        END DO
        CALL sub_1(eff_index) 
	    
        IF (eff_index .ge. mineff_index) THEN
            WRITE(10)(par(j),j=1,n_par),eff_index
        END IF
        
        IF (eff_index .gt. foptim) THEN
            foptim=eff_index
            gbest=par
        END IF
	    IF (i>=10) THEN
		    IF (MOD(i,INT(REAL(n_run)/10))==0 ) THEN
			    WRITE(*,1003) REAL(i)/REAL(n_run)*100.
		    END IF
	    END IF
    END DO
ELSEIF (runmode=='PSO') THEN
    
    ALLOCATE (x(n_par,n_particles),v(n_par,n_particles),pbest(n_par,n_particles))  ! initialization of PSO
    ALLOCATE (r(2*n_par),fit(n_particles),fitbest(n_particles))
    x=0; v=0;   
    r=0
    pbest=0; gbest=0; 
    fit=0; fitbest=0
    dw=(wmax-wmin)/n_run
    w=wmax

    ! random set of parameters during the first step
    CALL random_seed()
    CALL random_number(x)
    CALL random_number(v)
    DO j=1,n_par
        DO k=1,n_particles
            dxmax=(parmax(j)-parmin(j))     ! range of each parameter
	        dvmax=1.0*dxmax                 ! mmax velocity for each particle
	        x(j,k)=x(j,k)*dxmax+parmin(j)   ! random parameter value of j-parameter and k-particle
	        v(j,k)=v(j,k)*dvmax             ! random velocity value of j-parameter and k-particle
	        pbest(j,k)=x(j,k)               ! inizialization of partial best
        END DO
    END DO
    DO k=1,n_particles
        par=x(:,k)                          ! set of parameter of k-particle
        CALL sub_1(eff_index)        
        fitbest(k)=eff_index                ! best fit (max efficiency index)
    END DO
    CALL best(fit,k,foptim)                 
    DO j=1,n_par
	    gbest(j)=x(j,k)                     ! global best (overall best position)
    END DO

    DO i=1,n_run                            !number of iterations
        CALL random_seed()      
	    DO k=1,n_particles
		    CALL random_number(r)
		    status=0
            !update the velocity and the position of the particles
	        DO j=1,n_par
		        v(j,k)=w*v(j,k)+c1*r(j)*(pbest(j,k)-x(j,k))+c2*r(j+n_par)*(gbest(j)-x(j,k))
                x(j,k)=x(j,k)+v(j,k)
    		
		        ! absorbing wall boundary condition
		        IF (x(j,k).gt.parmax(j)) THEN
		            x(j,k)= parmax(j)
		            v(j,k)=0.0          
                    status=1            
		        END IF
                IF (x(j,k).lt.parmin(j)) THEN
		            x(j,k)= parmin(j)
		            v(j,k)=0.0          
                    status=1            
	            END IF
	        END DO

           ! new performances
            IF (status.eq.0) THEN	
                par=x(:,k)                          ! set of k-particles
	            CALL sub_1(eff_index) 
			    fit(k)=eff_index
			    ! write on file if efficiency index is greater than mineff_index
                IF (eff_index .ge. mineff_index) THEN
	                WRITE(10)(x(j,k),j=1,n_par),eff_index
                END IF
            ELSE
		        fit(k)=-1e30
            ENDIF
            
            ! Evaluation if the particle is improving its efficency
	        IF (fit(k).gt.fitbest(k)) THEN
			    fitbest(k)=fit(k)
			    DO j=1,n_par
				    pbest(j,k)=x(j,k)
			    END DO
		    END IF
	    END DO     
        
        ! Evaluation which is the best particle
        CALL best(fitbest,k,foptim)
	    DO j=1,n_par
		    gbest(j)=pbest(j,k)
	    END DO

        w=w-dw

	    IF (i>=10) THEN
		    IF (MOD(i,INT(REAL(n_run)/10))==0 ) THEN
			    WRITE(*,1003) REAL(i)/REAL(n_run)*100.
		    END IF
	    END IF
    	
	   ! If the norm between pbest and gbest is less then tol for the perc percentage of the particles --> exit the cycle
        count=0
        DO k=1,n_particles
            norm=0.
            DO j=1,n_par
                IF (flag_par(j)) THEN
                    norm=norm+( (pbest(j,k)-gbest(j))/(parmax(j)-parmin(j)) )**2
                END IF
            END DO
            norm=SQRT(norm)
            IF (norm .lt. 0.0) THEN                         ! 0.01 --> 1 °/oo
                count=count+1
            END IF
        END DO
        IF (count .ge. (1.0*n_particles)) THEN              ! 1 --> 100 °/o of n_particles
            WRITE(*,*)  '- Warning:  PSO has been stopped'
            EXIT
        END IF
    	    
    END DO
ENDIF

CLOSE(10)

par_best=gbest
finalfit=foptim

WRITE(*,*) 'Indice efficienza calibrazione', finalfit

! forward with the best set of parameters   
CALL forward


1003 FORMAT('Calcolo al', 1x, f5.1 ,1x,'%')
	
CALL CPU_TIME(T2)
print *, 'Computation time was ', T2-T1, 'seconds.'

STOP
END PROGRAM PSO
