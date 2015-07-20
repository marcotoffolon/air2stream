air2stream needs the following files to run:


1-----> 'input.txt'   to be located in the same folder of the executable file of air2stream. This file contains the input information:

! Main input
Switzerland		! name of the river/basin/region
MAH      		! name/ID of the air temperature station
2369			! name/ID of the water temperature station
c				! type of series: c=continuous series, m=mean year
1d        		! time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
8           	! version: 3,4,5,7,8 parameters
0				! Threshold temperature for ice formation
RMS				! objective function: RMS (Root Mean Square Error), NSE (Nash-Sutcliffe Efficiency Index), KGE (Kling-Gupta Efficiency Index)
0.0001			! minimum percentage of data: 0...1. E.g., when using 1m time resolution, the monthly average is set to NaN when available data during one month are less than this percentage  
RK4          	! mod_num : numerical method used to solve the model equation EUL (Euler Explicit), RK2 (Runge-Kutta 2), RK4 (Runge-Kutta 4), CRN (Crank Nicolson)
PSO            	! optimization algorithm: PSO (Particle Swarm Optimization) or RANSAM (Random Sampling)


2-----> 'PSO.txt'   to be located in the same folder of the executable file of air2stream. This file contains the parameters of PSO:
		see e.g.:	Kennedy, J., and R. Eberhart (1995), Particle swarm optimization, in Proceedings of IEEE International Conference on Neural Networks, pp. 1942-1948, doi:10.1109/ICNN.1995.488968.
					Robinson, J., and Y. Rahmat-Samii (2004), Particle swarm optimization in electromagnetics, IEEE T. Antenn. Propag., 52, 397-407, doi:10.1109/TAP.2004.823969.

! PSO parameters
500			! number of iterations
500			! number of particles
2 2	   		! constant for the equation of motion of the particles
0.9  0.4    ! inertia min and max in the equation of motion of particles
-999		! minimum efficiency index (i.e., RMS, NSE, KGE). All model realization with efficiency index greater than this threshold are saved in file 0_...


3-----> A folder named according to the first line in file input.txt (e.g., Switzerland), located in the same folder of the executable file of air2stream, and containing:

	a-----> input files for calibration 'IDair_IDwater_typeofseriesc.txt' (e.g., MAH_2369_cc.txt) and validation 'IDair_IDwat_typeofseriesv.txt' (e.g., MAH_2369_cv.txt)
		The two files have 6 columns (year, month, day, air temperature, water temperature, discharge):

		2001		1		1		-0.125		1.000		0.480
		2001		1		2		 0.158		1.155		0.588
		2001		1		3		 1.155		-999		0.567
		2001		1		4		 3.458		-999		0.898
		2001		1		5		10.125		3.254		1.025
		...
		
		NOTE: 	the series of observed air temperature must be complete. It cannot have gaps or no data. The series is at daily time scale.
				the series of observed water temperature can contain no-data (-999). The series is at daily time scale.
			    the series of observed must be completeIt cannot have gaps or no data. The series is at daily time scale.

	b-----> file 'parameters_air2stream.txt'
		The file contains the range of variation of each parameter of the model. The structure of the file is as follows:

		-5	-5	 -5	 -1	 0	 0 	 0	-1	-1
		 15	 1.5  5	  1	 20	 10	 1	 5	 2

		The first line contains the minimum value of each parameters: a1, a2, a3, a4, a5, a6, a7, a8.
		The first line contains the maximum value of each parameters: a1, a2, a3, a4, a5, a6, a7, a8.


-----------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------ OUTPUTs ---------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------

The model writes the output in a folder called 'air2stream' located in the same folder of the input files.
The folder 'air2stream' contains:

0-----> file '0_PSO_..........   .out'
	Binary file that contains a matrix with 9 columns. Each row contains the set of parameters and the efficiency of each iteration
	of the optimization algorithm.
	

1-----> file '1_PSO_..........   .out'
	ASCII file that contains the final best set of parameters founded and the values of error (indicated as the objective function chosen. Es:RMSE)

2-----> file '2_PSO_..........   .out'
	ASCII file from calibration period in which the columns are:
	year
	month
	day
	air temperature
	observed water temperature
	simulated water temperature at daily time scale
	observed water temperature aggregated at the time scale that was chosen in file 'cosa.txt' (Es:weekly)
	simulated water temperature aggregated at the time scale that was chosen in file 'cosa.txt' (Es:weekly)
	discharge at daily time scale
	
	
3-----> file '3_PSO_..........   .out'
     ASCII file with the same structure as file '2_PSO... .out' but it contains the outputs data of validation period.



-----------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------ post_processing.m ----------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------

To run the post_processing you should write the main informations in the code (from line 11 to line 20).

The post_processing produces these outputs:

---> 	dotty-plots of parameters

---> 	file errors in The folder errors that contains two columns. 
	The first column is about the errors in calibration periods and the second column contains the same errors in validation periods.
	The lines of these two columns are:
	
	RMS calculated at DAILY time scale
	MAE calculated at DAILY time scale
	KGE calculated at DAILY time scale
	NSE calculated at DAILY time scale
	alpha calculated at DAILY time scale
	beta calculated at DAILY time scale
	gamma calculated at DAILY time scale
	RMS calculated at WEEKLY time scale
	MAE calculated at WEEKLY time scale
	KGE calculated at WEEKLY time scale
	NSE calculated at WEEKLY time scale
	alpha calculated at WEEKLY time scale
	beta calculated at WEEKLY time scale
	gamma calculated at WEEKLY time scale
	RMS calculated at MONTHLY time scale
	MAE calculated at MONTHLY time scale
	KGE calculated at MONTHLY time scale
	NSE calculated at MONTHLY time scale
	alpha calculated at MONTHLY time scale
	beta calculated at MONTHLY time scale
	gamma calculated at MONTHLY time scale
	RMS calculated at SEASONAL time scale
	MAE calculated at SEASONAL time scale
	KGE calculated at SEASONAL time scale
	NSE calculated at SEASONAL time scale
	alpha calculated at SEASONAL time scale
	beta calculated at SEASONAL time scale
	gamma calculated at SEASONAL time scale

	In particular:

	RMSE = root mean squared error
	MAE = mean absolute error
	KGE = 1-sqrt(alpha^2+beta^2+gamma^2)
	alpha = std_mod/std_obs-1
	beta = mean_mod/mean_obs-1
	gamma = cov_mod/(std_mod*std_obs)-1
	NSE = Nash-Sutcliffe efficiency

---> plots of the hysteresys loop that contains the data of mean year smoothed using a 30 days moving average filter.

---> plot 'Ta-TWoss-TWsim_am_c': the mean year in calibration period.

---> plot 'Ta-TWoss-TWsim_am_c': the mean year in validation period.

---> plot 'Ta-TWoss-TWsim_d_calib': the daily series in calibration period.

---> plot 'Ta-TWoss-TWsim_d_valid': the daily series in validation period.

---> plot 'Ta-TWoss-TWsim_d_calib_smooth': the daily series in calibration period smoothed using a 30 days moving average filter.

---> plot 'Ta-TWoss-TWsim_d_valid_smooth': the daily series in validation period smoothed using a 30 days moving average filter.





