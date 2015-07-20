!-------------------------------------------------------------------------------
!				sub_1
!-------------------------------------------------------------------------------
SUBROUTINE sub_1(ei)
USE commondata

! Statement of variables
implicit none
REAL(KIND=8) :: ei

! call model air2stream
CALL call_model
	
! calculation of objective function
CALL funcobj(ei)

RETURN
END

!-------------------------------------------------------------------------------
!				MODEL AIR2STREAM
!-------------------------------------------------------------------------------

SUBROUTINE call_model

USE commondata

! Statement of variables
IMPLICIT NONE

INTEGER :: j
REAL(KIND=8) :: K1, K2 , K3, K4

IF (Twat_obs(1)==-999) THEN	            ! initial condition
	Twat_mod(1)=4.d0
ELSE
	Twat_mod(1)=Twat_obs(1)
END IF

DO j=1,n_tot-1
	
	    IF (mod_num .eq. 'RK2') THEN 
		    CALL RK4_air2stream(Tair(j), Q(j), Twat_mod(j), tt(j), K1)
            CALL RK4_air2stream(Tair(j+1), Q(j+1), Twat_mod(j) + K1, tt(j) + ttt, K2)   
		ENDIF
		
		IF (mod_num .eq. 'RK4') THEN
		    CALL RK4_air2stream(Tair(j), Q(j), Twat_mod(j), tt(j), K1)
	        CALL RK4_air2stream(0.5d0*(Tair(j) + Tair(j+1)), 0.5d0*(Q(j) + Q(j+1)), Twat_mod(j) + 0.5d0*K1, tt(j) + 0.5*ttt, K2)  !------------RK4
	        CALL RK4_air2stream(0.5d0*(Tair(j) + Tair(j+1)), 0.5d0*(Q(j) + Q(j+1)), Twat_mod(j) + 0.5d0*K2, tt(j) + 0.5*ttt, K3)  !------------RK4
	        CALL RK4_air2stream(Tair(j+1), Q(j+1), Twat_mod(j) + K3, tt(j) + ttt, K4)                                             !------------RK4
	    
	    ELSEIF (mod_num .eq. 'EUL') THEN 
		    CALL RK4_air2stream(Tair(j+1), Q(j+1), Twat_mod(j), tt(j+1), K1)
		
		ELSEIF (mod_num .eq. 'CRN') THEN  
		    IF ( version==8 .or. version==7 .or. version==4) THEN
                theta_j = Q(j)/(Qmedia) ;
                theta_j1 = Q(j+1)/(Qmedia) ;
                DD_j = theta_j**par(4);
                DD_j1 = theta_j1**par(4);
            
		        Twat_mod(j+1) = ( Twat_mod(j) + 0.5d0/DD_j*(par(1)+par(2)*Tair(j)-par(3)*Twat_mod(j)+ theta_j*(par(5)+par(6)*COS(2.d0*pi*((tt(j))-par(7)))-par(8)*Twat_mod(j))) + &
		        0.5d0/DD_j1*(par(1)+par(2)*Tair(j+1)+theta_j1*(par(5)+par(6)*COS(2.d0*pi*((tt(j+1))-par(7))) ) )) / (1.0d0+0.5d0*par(8)*theta_j1/DD_j1+0.5d0*par(3)/DD_j1);
		    
		    ELSEIF ( version==5 .or. version==3 ) THEN
		    
		        Twat_mod(j+1) = ( Twat_mod(j)*(1.0d0-0.5d0*par(3))+par(1)+0.5d0*par(2)*(Tair(j) + &
                Tair(j+1))+0.5d0*par(6)*COS(2.d0*pi*((tt(j))-par(7)))+0.5d0*par(6)*COS(2.d0*pi*((tt(j+1))-par(7))) )/(1.0d0+0.5d0*par(3));
		    ENDIF
		ENDIF
    
    
    IF (mod_num .eq. 'RK2') THEN 
		Twat_mod(j+1)=Twat_mod(j) + 0.5d0*(K1 + K2);               
	
	ELSEIF (mod_num .eq. 'RK4') THEN 
        Twat_mod(j+1)=Twat_mod(j) + 1.0d0/6.0d0*(K1 + 2.0d0*K2 + 2.0d0*K3 + K4 );	   
    
    ELSEIF (mod_num .eq. 'EUL') THEN    
	    Twat_mod(j+1)=Twat_mod(j) + K1;
	           
	ENDIF
   
    Twat_mod(j+1)=MAX(Twat_mod(j+1),Tice_cover);

END DO

RETURN
END



!-------------------------------------------------------------------------------
!				RUNGE KUTTA 4 - AIR2STREAM
!-------------------------------------------------------------------------------
SUBROUTINE RK4_air2stream(Ta, QQ, Tw, time, K)

USE commondata

! Statement of variables
IMPLICIT NONE

REAL(KIND=8), INTENT(OUT) :: K
REAL(KIND=8), INTENT(IN) :: Ta, Tw, time, QQ
REAL(KIND=8) :: DD

IF ( version==8 .or. version==4 ) THEN
	DD = (QQ / Qmedia )**par(4) ;	
ELSEIF ( version==5 .or. version==3 ) THEN
	DD = 0.0d0	
ELSE
	DD=1.0d0
END IF

      
IF ( version==3 ) THEN
  K = ( par(1) + par(2)*Ta -  par(3)*Tw )
END IF
      
IF ( version==5 ) THEN
  K = par(1) + par(2)*Ta -  par(3)*Tw + par(6)*COS(2.d0*pi*(time-par(7)))
END IF

IF ( version==8 .or. version==7 ) THEN
    K = par(1) + par(2)*Ta - par(3)*Tw  +  QQ/Qmedia * (  par(5) + par(6)*COS(2.d0*pi*(time-par(7))) - par(8)*Tw )
    K = K/DD
END IF

IF ( version==4 ) THEN
    K =( par(1) + par(2)*Ta - par(3)*Tw ) / DD;	
END IF

RETURN
END

!-------------------------------------------------------------------------------
!				CALCULATION OF THE OBJECTIVE FUNCTION
!-------------------------------------------------------------------------------
SUBROUTINE funcobj(ind)
! Statement of variables
USE commondata
IMPLICIT NONE

INTEGER :: i, j
REAL(KIND=8),INTENT(OUT) :: ind
REAL(KIND=8) :: TSS, mean_mod, TSS_mod, std_mod, covar_mod
REAL(KIND=8) :: tmp

Twat_mod_agg=-999
DO i=1,n_dat                            ! 1st year = warm-up period
    tmp=0.0d0
    DO j=I_inf(i,1),I_inf(i,2)
        tmp=tmp+Twat_mod(I_pos(j)) 
    END DO
    Twat_mod_agg(I_inf(i,3))=tmp/REAL(I_inf(i,2)-I_inf(i,1)+1)    
END DO
    
IF (fun_obj=='NSE') THEN
    TSS=0.d0
    DO i=1,n_dat		
	    TSS=TSS+(Twat_mod_agg(I_inf(i,3))-Twat_obs_agg(I_inf(i,3)))**2
    END DO
    ind=1.d0-TSS/TSS_obs
ELSEIF  (fun_obj=='KGE') THEN
    mean_mod=0.0d0
    DO i=1,n_dat		
	    mean_mod=mean_mod+Twat_mod_agg(I_inf(i,3))
    END DO
    mean_mod=mean_mod/REAL(n_dat)
    
    covar_mod=0.0d0
    TSS_mod=0.0d0
    DO i=1,n_dat
	    TSS_mod=TSS_mod+(Twat_mod_agg(I_inf(i,3))-mean_mod)**2
	    covar_mod=covar_mod+(Twat_mod_agg(I_inf(i,3))-mean_mod)*(Twat_obs_agg(I_inf(i,3))-mean_obs)
    END DO
    std_mod=DSQRT(TSS_mod/REAL(n_dat-1))
    covar_mod=covar_mod/REAL(n_dat-1)
    ind=1.0d0-DSQRT((std_mod/std_obs-1.0d0)**2 + (mean_mod/mean_obs-1.0d0)**2 + (covar_mod/(std_mod*std_obs)-1.0d0)**2)
ELSEIF  (fun_obj=='RMS') THEN
    TSS=0.d0
    DO i=1,n_dat
	    TSS=TSS+(Twat_mod_agg(I_inf(i,3))-Twat_obs_agg(I_inf(i,3)))**2
    END DO
    ind=-DSQRT(TSS/REAL(n_dat))
ELSE
	WRITE(*,*) 'Errore nella scelta della f. obiettivo'
END IF

RETURN
END
!------------------------------------------------------------------------------
!AGGREGATION (to calibrate the model with different time scale: daily, weekly, monthly, n-weeks )
!-------------------------------------------------------------------------------
SUBROUTINE aggregation

! Statement of variables
USE commondata
IMPLICIT NONE

INTEGER :: i, j, k, status, pp
INTEGER :: n_inf, n_pos, n_units, n_days, count, pos_tmp
INTEGER :: month, month_curr
INTEGER, DIMENSION(:,:), ALLOCATABLE :: A
INTEGER, DIMENSION(:), ALLOCATABLE :: B
REAL(KIND=8) :: tmp

pp=LEN_TRIM(time_res)
IF (pp==2) THEN
    unit=time_res(2:)
    READ(time_res,'(i1)') qty
ELSEIF (pp==3) THEN
    unit=time_res(3:)
    READ(time_res,'(i2)') qty
END IF

ALLOCATE(I_pos(n_tot),stat=status)
I_pos=-999
Twat_obs_agg=-999
n_inf=1   
n_pos=1
IF (time_res=='1d') THEN                    ! daily resolution
    n_units=n_tot-365
    ALLOCATE(I_inf(n_units,3))
    I_inf=-999
    DO i=366,n_tot                          ! 1st year = warm-up period
        IF (Twat_obs(i) .ne. -999) THEN
            I_inf(n_inf,2)=n_pos
            I_inf(n_inf,3)=i
            I_pos(n_pos)=i
            Twat_obs_agg(I_inf(n_inf,3))=Twat_obs(i)  
            n_inf=n_inf+1
            n_pos=n_pos+1
        END IF  
    END DO
ELSEIF (unit=='w') THEN                     ! weekly resolution
    n_days=qty*7      
    n_units=CEILING(REAL(n_tot-365)/REAL(n_days))
    ALLOCATE(I_inf(n_units,3))
    I_inf=-999
    DO i=366,n_tot,n_days
        tmp=0.0d0
        count=0
        pos_tmp=i+CEILING(0.5*n_days)-1
        DO j=0,n_days-1
            k=i+j
            IF (k .gt. n_tot) THEN
                EXIT
            END IF
            IF (Twat_obs(k) .ne. -999) THEN
                tmp=tmp+Twat_obs(k)
                I_pos(n_pos)=k
                n_pos=n_pos+1
                count=count+1
            END IF
        END DO
        IF (count .ge. n_days*prc) THEN           
            I_inf(n_inf,2)=n_pos-1
            I_inf(n_inf,3)=pos_tmp  
            Twat_obs_agg(I_inf(n_inf,3))=tmp/count
            n_inf=n_inf+1
        ELSE
            I_pos(n_pos-count:n_pos)=-999
            n_pos=n_pos-count
        END IF
    END DO
ELSEIF (unit=='m') THEN                     ! monthly resolution
    n_units=CEILING(n_tot/(30.5))
    ALLOCATE(I_inf(n_units,3))
    I_inf=-999
    n_days=0
    month_curr=-999    
    count=0
    DO i=366,n_tot
        month=date(i,2)
        IF (month .ne. month_curr) THEN
            IF (count .ge. n_days*prc .and. i .ne. 366) THEN  
                I_inf(n_inf,2)=n_pos-1
                I_inf(n_inf,3)=i-FLOOR(0.5*n_days)-1
                Twat_obs_agg(I_inf(n_inf,3))=tmp/count
                n_inf=n_inf+1
            ELSE
                I_pos(n_pos-count:n_pos)=-999
                n_pos=n_pos-count
            END IF        
            month_curr=month
            count=0           
            n_days=1
            tmp=0.0d0            
        ELSE
            n_days=n_days+1
        END IF  
        IF (Twat_obs(i) .ne. -999) THEN
            tmp=tmp+Twat_obs(i)
            I_pos(n_pos)=i
            n_pos=n_pos+1
            count=count+1
        END IF
    END DO
    ! Last month
    IF (count .ge. n_days*prc) THEN  
        I_inf(n_inf,2)=n_pos-1
        I_inf(n_inf,3)=i-FLOOR(0.5*n_days)-1
        Twat_obs_agg(I_inf(n_inf,3))=tmp/count
        n_inf=n_inf+1
    ELSE
        I_pos(n_pos-count:n_pos)=-999
        n_pos=n_pos-count
    END IF                    
ELSE
    WRITE(*,*) 'Error: variable time_res'
END IF

n_dat=n_inf-1
n_pos=n_pos-1

I_inf(1,1)=1
I_inf(2:n_dat,1)=I_inf(1:n_dat-1,2)+1
ALLOCATE(A(n_units,3),B(n_tot))
A=I_inf;    
B=I_pos
DEALLOCATE(I_inf,I_pos)
ALLOCATE(I_inf(n_dat,3),I_pos(n_pos))
I_inf=A(1:n_dat,:)
I_pos=B(1:n_pos)
DEALLOCATE(A,B)

RETURN
END
!-------------------------------------------------------------------------------
!				STATIS (to calculate errors) 
!-------------------------------------------------------------------------------
SUBROUTINE statis

! Statement of variables
USE commondata
IMPLICIT NONE

INTEGER :: i, k, d, status

mean_obs=0.d0
TSS_obs=0.d0
DO i=1,n_dat		    
	mean_obs=mean_obs+Twat_obs_agg(I_inf(i,3))
END DO
mean_obs=mean_obs/REAL(n_dat)

DO i=1,n_dat		
	TSS_obs=TSS_obs+(Twat_obs_agg(I_inf(i,3))-mean_obs)**2
END DO

std_obs=DSQRT(TSS_obs/REAL(n_dat-1))

RETURN
END

!-------------------------------------------------------------------------------
!				FORWARD
!-------------------------------------------------------------------------------
SUBROUTINE forward
USE commondata
IMPLICIT NONE

INTEGER :: i, j
REAL(KIND=8) :: ei_check,ei

par=par_best		                ! best set of parameters

CALL call_model

! Control: efficency index
CALL funcobj(ei_check)
IF (ABS(ei_check - finalfit) .gt. 0.0001) THEN
	WRITE(*,*) 'Errore efficienza in forward'
	WRITE(*,*) ei_check, finalfit
	PAUSE
ELSE
	WRITE(*,*) 'Controllo superato'
END IF
	
OPEN(UNIT=12,FILE=TRIM(folder)//'/2_'//TRIM(runmode)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'c_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')

DO i=1,n_tot
WRITE(12,1005) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),Twat_obs_agg(i),Twat_mod_agg(i),Q(i)
END DO

CLOSE(12)

CALL read_validation

IF (n_tot .lt. 365) THEN
    ei=-999
    GO TO 200
END IF

! aggregate calibration data on the basis of time_resolution
CALL aggregation

CALL statis
WRITE(*,*) 'mean, TSS and standard deviation (validation)'
WRITE(*,*)  SNGL(mean_obs),SNGL(TSS_obs),SNGL(std_obs)

CALL call_model

CALL funcobj(ei)

OPEN(UNIT=13,FILE=TRIM(folder)//'/3_'//TRIM(runmode)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'v_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
DO i=1,n_tot
WRITE(13,1005) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),Twat_obs_agg(i),Twat_mod_agg(i),Q(i)
END DO
CLOSE(13)

200 OPEN(UNIT=11,FILE=TRIM(folder)//'/1_'//TRIM(runmode)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
DO i=1,n_par
	WRITE(11,*) par_best(i)
END DO
WRITE(11,*) finalfit
WRITE(11,*) ei
CLOSE(11)


1004 FORMAT(i4,1x,i4,1x,i4,1x,5(1x,f10.5))
1005 FORMAT(i4,1x,i4,1x,i4,1x,6(1x,f10.5))

RETURN
END

!-------------------------------------------------------------------------------
!				BEST
!-------------------------------------------------------------------------------
SUBROUTINE best(fit,part,foptim)
USE commondata
IMPLICIT NONE

INTEGER,INTENT(OUT)::part
INTEGER:: k
REAL(KIND=8),INTENT(IN),DIMENSION(n_particles):: fit
REAL(KIND=8),INTENT(OUT):: foptim

foptim=-1e30				                ! really small value
DO k=1,n_particles
    IF(fit(k).gt.foptim) then
	    foptim=fit(k)
	    part=k
    END IF
END DO

RETURN
END

!-------------------------------------------------------------------------------
!				LEAP YEAR
!-------------------------------------------------------------------------------
SUBROUTINE leap_year(Y,I)
USE commondata
IMPLICIT NONE

INTEGER,INTENT(IN) :: Y
INTEGER,INTENT(OUT) :: I

IF(MOD(Y,100).NE.0.AND.MOD(Y,4).EQ.0) THEN 
    I=1
ELSEIF(MOD(Y,400).EQ.0) THEN
    I=1
ELSE 
    I=0
END IF

RETURN
END