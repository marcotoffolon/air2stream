SUBROUTINE read_calibration

USE commondata
USE ifport                          ! necessary for makedirqq

IMPLICIT NONE
INTEGER:: i, j, status, pp
CHARACTER(LEN=1) :: string
LOGICAL result                      ! necessary for makedirqq

! read input information
OPEN(unit=1,file='input.txt',status='old',action='read')
READ(1,*)		                    ! header
READ(1,*) name                      ! name of folder
READ(1,*) air_station               ! name/ID of the air station
READ(1,*) water_station             ! name/ID of the water station
READ(1,*) series                    ! type of series: c=continuous, m=mean year
READ(1,*) time_res                  ! time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
READ(1,*) version                   ! version: 3,4,5,7,8 parameters
READ(1,*) Tice_cover                ! Threshold temperature for ice formation
READ(1,*) fun_obj                   ! objective function: KGE, NSE, RMS
READ(1,*) prc                       ! minimum percentage of data in input: 0...1
READ(1,*) mod_num                   ! mod_num :   RK4 , EUL , RK2 , CRN
READ(1,*) runmode                   ! optimization algorithm: PSO or RANSAM
CLOSE(1)

model = 'air2stream'
station=TRIM(air_station)//'_'//TRIM(water_station)

WRITE(*,*) 'Objective function ',fun_obj

! read PSO parameters
OPEN(unit=1,file='PSO.txt',status='old',action='read')
READ(1,*)		                    ! header
READ(1,*) n_run                     ! number of iterations
READ(1,*) n_particles               ! number of particles
READ(1,*) c1,c2                     ! constant for the equation of motion of the particles
READ(1,*) wmax,wmin                 ! inertia min and max in the equation of motion of particles
READ(1,*) mineff_index              ! index for the minimum efficiency that the code memorizes
CLOSE(1)

! read model parameters
    OPEN(unit=1,file=TRIM(name)//'/parameters_air2stream.txt',status='old',action='read')

    READ(1,*) (parmin(i),i=1,n_par);	
    READ(1,*) (parmax(i),i=1,n_par);
    
! parameters that are not used are zeroed    
    flag_par=.true.
    IF (version == 3) THEN                                      !air2stream with 3 parameters
        parmin(4)=0;    parmax(4)=0;    flag_par(4)=.false.;      
	    parmin(5)=0;	parmax(5)=0;    flag_par(5)=.false.;
	    parmin(6)=0;	parmax(6)=0;    flag_par(6)=.false.;
	    parmin(7)=0;	parmax(7)=0;    flag_par(7)=.false.;
	    parmin(8)=0;	parmax(8)=0;    flag_par(8)=.false.; 
	END IF 
	IF (version == 4) THEN                                      !air2stream with 4 parameters
	     parmin(5)=0;    parmax(5)=0;    flag_par(5)=.false.;    
	     parmin(6)=0;    parmax(6)=0;    flag_par(6)=.false.; 
	     parmin(7)=0;    parmax(7)=0;    flag_par(7)=.false.; 
	     parmin(8)=0;    parmax(8)=0;    flag_par(8)=.false.; 
	END IF
    IF (version == 5) THEN                                      !air2stream with 5 parameters
  	    parmin(4)=0;	parmax(4)=0;    flag_par(4)=.false.;
	    parmin(5)=0;	parmax(5)=0;    flag_par(5)=.false.;
	    parmin(8)=0;	parmax(8)=0;    flag_par(8)=.false.;
	ENDIF
	IF (version == 7) THEN                                      !air2stream with 7 parameters
         parmin(4)=0;    parmax(4)=0;    flag_par(4)=.false.;      
    END IF 
    IF (version == 4) THEN                                      !air2stream with 8 parameters
	     parmin(5)=0;    parmax(5)=0;    flag_par(5)=.false.;    
	     parmin(6)=0;    parmax(6)=0;    flag_par(6)=.false.; 
	     parmin(7)=0;    parmax(7)=0;    flag_par(7)=.false.; 
	     parmin(8)=0;    parmax(8)=0;    flag_par(8)=.false.; 
	END IF
	
    CLOSE(1)


! write parameters
WRITE(string,'(i1)' ) version
folder = TRIM(name)//'/'//TRIM(model)
result=makedirqq(folder)
folder = TRIM(name)//'/'//TRIM(model)//'/output_'//mod_num//'_'//string//'/'
result=makedirqq(folder)

OPEN(unit=2,file=TRIM(folder)//'/parameters_air2stream.txt',status='unknown',action='write')

WRITE(2,'(I2,A)') n_par, '   !numero parametri'
WRITE(2,'(<n_par>(F10.5,1x))') (parmin(i),i=1,n_par)
WRITE(2,'(<n_par>(F10.5,1x))') (parmax(i),i=1,n_par)
CLOSE(2)

! read T series (calibration)
CALL read_Tseries('c')

! open file for the writing of all parameter set + efficiency index
OPEN(unit=10,file=TRIM(folder)//'/0_'//TRIM(runmode)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',status='unknown',action='write',form='binary')

RETURN
END

!-------------------------------------------------------------------------------
!				READ VALIDATION PERIOD
!-------------------------------------------------------------------------------
SUBROUTINE read_validation

USE commondata

IMPLICIT NONE

DEALLOCATE(date, tt, Tair, Twat_obs, Twat_obs_agg, Twat_mod, Twat_mod_agg, Q)
DEALLOCATE(I_pos, I_inf)

! read T series (validation)
CALL read_Tseries('v')

RETURN
END

!-------------------------------------------------------------------------------
!				READ TEMPERATURE FILE
!-------------------------------------------------------------------------------
SUBROUTINE read_Tseries(p)

USE commondata

IMPLICIT NONE

INTEGER :: i, j, k, status
INTEGER :: n_year, leap, year_ini
CHARACTER(LEN=1),INTENT(IN) :: p
CHARACTER(LEN=10) :: period

n_tot=0;

IF (p=='c') THEN
    period='calibration'
ELSE
    period='validation'
END IF

OPEN(unit=3,file=TRIM(name)//'/'//TRIM(station)//'_'//series//p//'.txt',status='unknown',action='read', iostat=status)
openif3: IF (status==0) THEN
	readloop3: DO
		READ(3,*,iostat=status)
		IF (status/=0) EXIT
		n_tot=n_tot+1
	END DO readloop3
	readif3: IF(status>0) THEN
	END IF readif3	
END IF openif3
REWIND(3)

! allocation + replication of the 1st year
WRITE(*,1001)  n_tot/365.25,TRIM(period)
1001 FORMAT('There are ',f4.1,' years for ', a12)

IF (p=='v' .and. n_tot .lt. 365) THEN
    WRITE(*,*) 'Validation period < 1 year --> validation is skipped'
    GO TO 100
END IF
n_year=CEILING(n_tot/365.25)
n_tot=n_tot+365             ! the 1st year is replicated. The 1st year is always considered 365 days long
ALLOCATE(date(n_tot,3),stat=status)
ALLOCATE(Tair(n_tot),stat=status)
ALLOCATE(Twat_obs(n_tot),stat=status) 
ALLOCATE(Twat_obs_agg(n_tot),stat=status) 
ALLOCATE(Twat_mod(n_tot),stat=status) 
ALLOCATE(Twat_mod_agg(n_tot),stat=status) 
ALLOCATE(tt(n_tot),stat=status)
ALLOCATE(Q(n_tot),stat=status) 

        Qmedia = 0.0d0
        n_Q = 0
        
        DO i=366,n_tot
	        READ(3,*) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Q(i)
	        IF ( Q(i) .ne. -999 ) THEN
	        n_Q = n_Q + 1
	        Qmedia = Qmedia + Q(i) 
	        END IF
        END DO
            Qmedia = Qmedia / REAL(n_Q)
           
 
year_ini=date(366,1)
date(1:365,:)=-999
Tair(1:365)=Tair(366:730)
Twat_obs(1:365)=Twat_obs(366:730)
Q(1:365)=Q(366:730)


CLOSE(3)

! check leap years + define tt
k=0
DO j=1,365
    tt(k+j)=REAL(j)/365.0d0
END DO
k=365
DO i=1,n_year
    CALL leap_year(year_ini+i-1,leap)
    IF(leap==0) THEN
        DO j=1,365
            IF (k+j .gt. n_tot) THEN
                EXIT
            END IF
            tt(k+j)=REAL(j)/365.0d0
        END DO
        k=k+365
    ELSE
        DO j=1,366
            IF (k+j .gt. n_tot) THEN
                EXIT
            END IF
            tt(k+j)=REAL(j)/366.0d0
        END DO
        k=k+366
    END IF
END DO

100 RETURN 
END