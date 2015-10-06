       .__       ________            __                                 
_____  |__|______\_____  \   _______/  |________   ____ _____    _____  
\__  \ |  \_  __ \/  ____/  /  ___/\   __\_  __ \_/ __ \\__  \  /     \ 
 / __ \|  ||  | \/       \  \___ \  |  |  |  | \/\  ___/ / __ \|  Y Y  \
(____  /__||__|  \_______ \/____  > |__|  |__|    \___  >____  /__|_|  /
     \/                  \/     \/                    \/     \/      \/       
A model to predict River Water Temperature (RWT) using air temperature and discharge.
Version 1.0.0 - October 2015

Provided by Sebastiano Piccolroaz and Marco Toffolon

Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
email contacts: s.piccolroaz@unitn.it, marco.toffolon@unitn.it

How to cite:

- Toffolon M., and S. Piccolroaz (2015), A hybrid model for river water temperature as a function of air temperature and discharge, Environmental Research Letters, under review

Related articles (air2water: A model to predict Lake Surface Temperature using air temperature):

- Piccolroaz S., M. Toffolon, and B. Majone (2013), A simple lumped model to convert air temperature into surface water temperature in lakes, Hydrol. Earth Syst. Sci., 17, 3323-3338, doi:10.5194/hess-17-3323-2013

- Toffolon M., S. Piccolroaz, B. Majone, A.M. Soja, F. Peeters, M. Schmid and A. Wüest (2014), Prediction of surface water temperature from air temperature in lakes with different morphology, Limnology and Oceanography, 59(6), 2185-2202, doi: 10.4319/lo.2014.59.6.2185

- Toffolon M., S. Piccolroaz, and B. Majone (2015), The role of stratification on lakes' thermal response: The case of Lake Superior, Water Resources Research, doi: 10.1002/2014WR016555, in press
-----------------------------------------------------------------------------------------------------------------------------------------
 
How to use air2stream 

0. Preamble
The model is presented and discussed in Toffolon and Piccolroaz, 2015.

Precompiled executable files are available for Linux (air2stream_1.0.0.out) and Windows (air2stream_1.0.0.exe) systems.

The code is provided with a few example test cases. The study sites are three Swiss rivers characterized by different hydrological conditions:
1. The river Mentue flows at low altitude on the Swiss plateau through a sparsely inhabited area with predominant agricultural land use, unaffected by strong anthropic thermal alterations. This case is considered as a `natural low-land' river. (river station ID: 2369; meteo station ID: MAH)
2. The station of the river Rhone at Sion lies at the bottom of a populated alpine valley. Starting from the beginning of the 20th century its hydrological regime has been altered by the construction of a large high-head hydropower storage system, and the river is now affected by strong hydro- and thermo-peaking. This case is considered as representative of a `regulated' river. river station ID: 2011; meteo station ID: SIO)
3. The river Dischmabach is located at high altitude in a steep glacial valley that is uninhabited and used for mountain pastures, with a significant influence of snow melting. This case is taken as a `snowmelt-fed' river. river station ID: 2327; meteo station ID: DAV)
Data have been downloaded from (downloaded data have been pre-processed):
- Federal Office for the Environment (FOEN) (webpage: http://www.bafu.admin.ch/wasser/13465/13483/14087/index.html?lang=en) --> daily water temperature and discharge;
- Meteoswiss (webpage: http://www.meteoswiss.admin.ch/home/measurement-and-forecasting-systems/land-based-stations/automatisches-messnetz.html) --> daily air temperature.


1. Input files
-----> 'input.txt' to be located in the same folder of the executable file of air2stream. This file contains the input information:
! Main input
Switzerland		! name of the river/basin/region
MAH      		! name/ID of the air temperature station
2369			! name/ID of the water temperature station
c				! type of series: c=continuous series, m=mean year
1d        		! time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
8           	! version: 3,4,5,7,8 parameters
0				! Threshold temperature for ice formation
RMS				! objective function: RMS (Root Mean Square Error), NSE (Nash-Sutcliffe Efficiency Index), KGE (Kling-Gupta Efficiency Index)
RK4          	! mod_num : numerical method used to solve the model equation EUL (Euler Explicit), RK2 (Runge-Kutta 2), RK4 (Runge-Kutta 4), CRN (Crank-Nicolson)
PSO            	! run mode: PSO (Particle Swarm Optimization), LATHYP (Latin Hypercube), or FORWARD
0.60			! minimum percentage of data: 0...1. E.g., when using 1m time resolution, the monthly average is set to NaN when available data during one month are less than this percentage (60% in this case)   
500				! n_run: number of iterations
-999			! minimum efficiency index (i.e., RMS, NSE, KGE). All model realization with efficiency index greater than this threshold are saved in file 0_...

-----> 'PSO.txt'   to be located in the same folder of the executable file of air2stream. This file contains the parameters of PSO:
see e.g.:	
- Kennedy, J., and R. Eberhart (1995), Particle swarm optimization, in Proceedings of IEEE International Conference on Neural Networks, pp. 1942-1948, doi:10.1109/ICNN.1995.488968.
- Robinson, J., and Y. Rahmat-Samii (2004), Particle swarm optimization in electromagnetics, IEEE T. Antenn. Propag., 52, 397-407, doi:10.1109/TAP.2004.823969.

! PSO parameters
500			! number of particles
2 2	   		! c1, c2: constants known as cognitive and social learning factors (used in the equation of motion of particles)
0.9  0.4    ! inertia max and min (used in the equation of motion of particles)

-----> A folder named according to the first line in file input.txt (e.g., Switzerland), located in the same folder of the executable file of air2stream and containing:

-----> input files for calibration 'IDair_IDwater_typeofseriesc.txt' (e.g., MAH_2369_cc.txt) and validation 'IDair_IDwat_typeofseriesv.txt' (e.g., MAH_2369_cv.txt)
	The two files have 6 columns (year, month, day, air temperature, water temperature, discharge):

	2001		1		1		-0.125		1.000		0.480
	2001		1		2		 0.158		1.155		0.588
	2001		1		3		 1.155		-999		0.567
	2001		1		4		 3.458		-999		0.898
	2001		1		5		10.125		3.254		1.025	
	...
		
	NOTE: 	The series of observed air temperature and discharge must be complete. They cannot have gaps or no data. 
		The series of observed water temperature can contain no-data (-999). 
		Both series are always at daily time scale, as the equation of the model is solved with daily time step. The model automatically evaluates weekly, multi-weekly, or monthly averages (of water temperature) when using different time scales for model calibration. 
	
-----> file 'parameters.txt'
	The file contains the range of variation of each parameter of the model. 
	The range of parameters should be carefully checked during the post-processing step through the analysis of the dotty plots (the cloud of points obtained by performing latin hypercube should be centred within the searching domain, and in any case the best set of parameters should be sufficiently far from the upper and lower bounds.
	The structure of the file 'parameters.txt' is as follows:

	-5   -5   -5   -1   0    0    0   -1	 
	15   1.5   5    1   20   10   1    5	

	The first line contains the minimum value of each parameters, while the second contains the maximum values. 
		
-----> file 'parameters_forward.txt'
	The file contains a set of parameters to be used to run the model in forward mode (required only if the model is run in forward), see 10th line of file input.txt . The structure of the file is as follows:

	 0.000284   0.006668   0.006719   2.888324   0.017322   0.212055 148.709752   0.499891  -1.303755

	
2. Output files
The model writes the output in a folder called 'output_numberofparameters' (e.g., output_8) located in the same folder of the input files (e.g., Switzerland).
The folder contains:

-----> file '0_PSO_..........   .out'
	Binary file that contains a matrix with 9 columns. Each row contains the set of parameters (columns 1-8) and the associated efficiency index (column 9) of each iteration performed by the optimization algorithm. Values are saved in double precision.
	This file allows: a) drawing the dotty plots of parameters to evaluate whether the a priori range of variations has been reasonably defined (i.e., not too narrow, not too large), b) evaluating whether the optimization (searching) algorithm converged towards the best solution, c) perform uncertainty analyses (when using RANSAM or LATHYP).
	
-----> file '1_PSO_..........   .out'
	ASCII file that contains the best set of parameters (1st line), the value of the efficiency index during the calibration period (2nd line), and the value of the efficiency index during the validation period (if any, 3nd line)

-----> file '2_PSO_..........   .out'
	ASCII file containing the results of the calibration period in which the columns are:
	year
	month
	day
	observed air temperature
	observed water temperature
	simulated water temperature
	observed water temperature aggregated at the time scale chosen in file 'input.txt' (e.g., 1m) 
	simulated water temperature aggregated at the time scale chosen in file 'input.txt' (e.g., 1m) 
	discharge at daily time scale
	
	Note that the first year is replicated and is used as warm up year to mitigate the influence of initial conditions. 
	During the warm up year, year, month and day are equal to -999.
	
-----> file '3_PSO_..........   .out'
	ASCII file with the same structure of file '2_PSO... .out' but referred to the validation period. 


4. Post-processing
A matlab script is available for post-processing. The script ('post_processing.m') is located in the folder 'post_processing'. The output of the script is saved in the same folder where the output of air2stream is located. 
The figures produced by the post-processing script are:
- dotty plots of parameters;
- comparison between observed air temperature, observed water temperature, and simulated water temperature for the calibration period;
- comparison between observed air temperature, observed water temperature, and simulated water temperature for the validation period.
