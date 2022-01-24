README
------

These scripts are designed to
      - read in the simulation 2-dimensional maps of total matter, stellar matter and x-ray maps,
      - find the peaks in each map
      - associate and match the peaks to create configrurations of three mass components
      - create a vector between the gas and the stellar matter and work out where the dark matter is closest along that vector.
      - stack each of these vectors and find beta

Key terms
---------


Stellar		Inception	  Xray 
Matter Peak       Point          Gas Peak
    'S' ------------'I'---------'G'
		     |         /
		     |        /
		     |	     /
		     |      /
		     |     /
		     |    /
		     |   /
		     |  /
		     | /
		     |/
		    'D'
        	Dark Matter Peak



dist_sd / SD : SD always referes to the distance between "stars" and "dark matter"
dist_di / DI : T


Simulations
----------
CDM+baryons :  sigma/dm = 0, with the fiducial baryonic prescipriton
SIDM0.1+baryons : sigma/dm=0.1, with the fiducial baryonic prescipriton
SIDM0.3+baryons : sigma/dm=0.3, with the fiducial baryonic prescipriton
SIDM1+baryons : sigma/dm=1, with the fiducial baryonic prescipriton
CDM_low :  sigma/dm = 0, with the lowest level of AGN reheating that is compatible with obseravtions
CDM_hi :  sigma/dm = 0, with the highest level of AGN reheating that is compatible with obseravtions

Halos
---------
Each simulations have 300 of the most massive haloes from four redshift slice: 0., 0.125, 0.250, 0.375

Maps
--------
Each TOTAL MATTER and STELLAR MATTER maps are 2 Mpc x 2Mpc with a resolution of 1kpc
Each XRAY map is 10Mpc x 10Mpc with 5kpc resolution

Scripts
-------

HSTconstraints.py :

execution
>> python HSTconstraints.py [DM model]

This script is the first one to make sure that all the files are created and found properly. It will go through the different classes, create them, make sure data is there, determine the vectors and plot BETA as a function of cross-section,

REQUIREMENTS
-------------
pyraf
astropy
cosmolopy
numpy
scipy
matplotlib ( % PyQt5)
