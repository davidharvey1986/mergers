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



dist_XY : refers to the magnitude of the vector between component X and Y
S : Stellar Matter
D : Dark Matter
G : Gas
I : Inception point (closest point dark matter is along the SG vector)

Beta = dist_SI / dist_SG (the signal)
Beta_perp = dist_DI / dist_SG (should be consistent with zero)


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

INSTALL
-------

virtualenv mergers -p /usr/local/bin/python3.7
cd mergers
source bin/activate
unset PYTHONPATH
pip install pyraf
pip install astropy
pip install cosmolopy
pip install numpy
pip install matplotlib
pip install SciencePlots
pip install PyQt5
git clone git@github.com:davidharvey1986/mergers.git
cd mergers

python HSTconstraints.py




Scripts
-------

HSTconstraints.py :

execution
>> python HSTconstraints.py [DM model]

This script is the first one to make sure that all the files are created and found properly. It will go through the different classes, create them, make sure data is there, determine the vectors and plot BETA as a function of cross-section. It doesn't require any inout and if none is given will default to all the simulations.

>> python pl

REQUIREMENTS
-------------
pyraf
astropy
cosmolopy
numpy
scipy
matplotlib ( % PyQt5)
