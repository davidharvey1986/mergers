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
git clone git@github.com:davidharvey1986/mergers.git
cd mergers

python HSTconstraints.py




Two Main Scripts
-------

HSTconstraints.py
plotBetaAsFunctionMass.py

Class Files
-----------

clusterClass : a class for an individual cluster from a simulation at a redshift. This gets the component positions for an individual cluster and matches them. Does the majority of the work.

ClusterSample (inside getMergerSample.py): a class that compiles the ensemble of all clusterClasses and carries out calculations on clusterclass and compiles in to vectors.

Other files
----------

component_extractor.py : uses pysex.py and SEXtractor to get the positions of all peaks in the input fits file.

getAndPlotTrend.py : some simple linear regression scripts

MockedClusterClass.py : creates a fake cluster. Hasnt been tested in a while.

matchMassComponents.py : match two lists of mass components anbd return a matched catalogue

REQUIREMENTS
-------------
pyraf
astropy
cosmolopy
numpy
scipy
matplotlib ( % PyQt5)
