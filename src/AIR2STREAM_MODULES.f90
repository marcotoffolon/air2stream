MODULE commondata

IMPLICIT NONE
SAVE 
! Statement of variables
INTEGER, PARAMETER :: n_par = 8    
INTEGER :: n_Q      
REAL(KIND=8) :: Qmedia, theta_j, theta_j1, DD_j, DD_j1    
REAL(KIND=8), PARAMETER :: pi = ACOS(0.d0)*2.d0
REAL(KIND=8), PARAMETER :: ttt = 1.0d0/365.0d0

INTEGER :: n_tot, n_dat
INTEGER :: version, qty
INTEGER, ALLOCATABLE, DIMENSION(:) :: I_pos
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_inf
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: date
REAL(KIND=8) :: Tice_cover, prc
REAL(KIND=8) :: mean_obs, TSS_obs, std_obs
REAL(KIND=8) :: mineff_index,finalfit
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: tt, Tair, Twat_obs_agg, Twat_obs, Q
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Twat_mod, Twat_mod_agg

REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: parmin
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: parmax
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: par, par_best
CHARACTER(len=100) :: folder
CHARACTER(LEN=30) :: name, air_station, water_station, station, model, runmode
CHARACTER(LEN=1) :: series, unit
CHARACTER(LEN=3) :: time_res
CHARACTER(LEN=3) :: fun_obj, mod_num
LOGICAL,ALLOCATABLE,DIMENSION(:) :: flag_par
INTEGER :: n_run,n_particles                    ! number of run and number of particles for PSO (Particle Swarming Optimization)
REAL(KIND=8) ::c1,c2,wmin,wmax                  ! parameters for the equation of motion of particles in PSO

END MODULE commondata

