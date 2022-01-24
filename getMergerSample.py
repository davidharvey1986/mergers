
from clusterClass import *
import getSimData as gSD
import os as os
import glob as glob
import pickle as pkl
import xray_sim_concentration as xsc
from MockedClusterClass import *

class ClusterSample:

    def __init__( self, SimName, **kwargs ):
        self.SimName = SimName
        self.ClusterSample = []
        if 'Mocked' in SimName:
            self.getMockedSample( **kwargs )
        
        else:
            self.Sims = gSD.getSimList( SimName )
            if 'CDM_' in SimName:
                self.dataDir = \
                    '/data2/harvey/simulations/BAHAMAS'
            else:
                self.dataDir = \
                    '/data1/harvey/simulations/BAHAMAS'

            self.pickleFile = \
              'Pickles/'+SimName+'_MassPositions.pkl'
            
            if os.path.isfile( self.pickleFile):
                self.ClusterSample = \
                  pkl.load( open( self.pickleFile, 'rb'),
                                      encoding='latin1')
            else:
                self.getClusterSample()
        print("Number of ints", len([ i.clusterInt for i in self.ClusterSample]))

    def getClusterSample( self ):
        
        for iSim in self.Sims:
          
                self.SimRoot = self.dataDir+'/'+iSim
                self.getRedshiftList( iSim )

                for iRedshift in self.redshiftList:
                    simPickleFile = \
                      "Pickles/"+iSim+"_"+str(iRedshift)+".pkl"

                    if os.path.isfile( simPickleFile ):
                        simClusterClasses = \
                          pkl.load( open(simPickleFile, 'rb'),
                                      encoding='latin1')
                        [ self.ClusterSample.append( i ) for i in \
                        simClusterClasses ]
                        self.getClusterInts( iRedshift )
                    else:
                        self.getClusterInts( iRedshift )

                        self.getGasConcentration( iSim, iRedshift)
                
                        for iCluster, iClusterInt in \
                          enumerate(self.clusterInts):

                            iClusterClass = \
                            clusterClass( iClusterInt, iSim, iRedshift,\
                                    self.GasConcentration[iCluster])
                        
                            self.ClusterSample.append( iClusterClass )
                        
                        pkl.dump( self.ClusterSample, open(simPickleFile,'wb'))

        pkl.dump( self.ClusterSample, open(self.pickleFile, 'wb') )
                    
                                                 
            

    def getGasConcentration( self, sim, redshift ):
        redshiftStr = "z_%0.3f" % redshift
        
        self.GasConcentration = np.array([0.2 for i in self.clusterInts]) #\
        
          #xsc.get_all_xray_concentrations( self.clusterInts, \
           #                             sim=sim,\
           #                             clusterRedshiftStr=redshiftStr )
        

        
    def getRedshiftList( self, iSim ):

        redshiftListComplete = glob.glob(self.SimRoot+'/z_*')
        print(redshiftListComplete)

        self.redshiftList = [ np.float(i.split('/')[-1].split('_')[1])\
                              for i in redshiftListComplete if 'tar' not in i]


    def getClusterInts( self, redshift ):
        redshiftStr = "%0.3f" % redshift
        self.GalaxyCatalogs = \
          glob.glob(self.SimRoot+'/z_'+redshiftStr+'/HIRES_MAPS/*stellar*sph*')

        self.clusterInts = [ i.split('/')[-1].split('_')[1] \
                                 for i in self.GalaxyCatalogs ]

        
    def extractMergerHalos( self ):

        mergerHalos = []
        pickleFile = "Pickles/"+self.SimName+"_mergerHalos.pkl"

        if os.path.isfile( pickleFile ):
            self.mergerHalos = pkl.load(open(pickleFile, 'rb'),
                                      encoding='latin1')
        else:
            
            for iCluster in self.ClusterSample:

                

                iClusterPickle = 'Pickles/individual/'+iCluster.simulation+\
                    str(iCluster.redshift)+str(iCluster.clusterInt)+'.pkl'
                if os.path.isfile( iClusterPickle ):
                    iCombined = pkl.load(open(iClusterPickle,'rb'),
                                      encoding='latin1')
                else:
                    try:
                        flag = iCluster.combineMassComponents()
                    except:
                        continue
                    iCluster.flag = flag
                    pkl.dump(iCluster, open(iClusterPickle,'wb'))
                    iCombined = pkl.load(open(iClusterPickle,'rb'),
                                      encoding='latin1')
                if iCombined.flag > 0:
                    iCombined.calculateBulletVectors()
                    mergerHalos.append( iCombined )

                
                    
            self.mergerHalos = mergerHalos
            
            pkl.dump( mergerHalos, open(pickleFile, 'wb'))
            
    def CalculateOffsetVectors( self, nClusters=3 ):

        self.dist_si = np.array([])
        self.dist_sg = np.array([])
        self.dist_di = np.array([])
        self.beta = np.array([])
        self.betaPerp = np.array([])
        self.ClusterMass = np.array([])
        self.dist_sd = np.array([])
        self.componentMass = np.array([])
        
        for iMerger in self.mergerHalos:
            if  (len(iMerger.dist_si)<=1):
                continue


            #hange this to the first nClusteres not get rid if roo many
            if nClusters != -1:
                if len(iMerger.dist_si) > nClusters:
                    continue
            
            BinaryHalos = iMerger.BinaryCluster()
            

            self.dist_si = np.append( self.dist_si, \
                        iMerger.dist_si)
            self.dist_sg = np.append( self.dist_sg, \
                        iMerger.dist_sg)
            self.dist_di = np.append( self.dist_di, \
                        iMerger.dist_di)
            sd = np.sqrt(np.sum(iMerger.vector_sd**2, axis=0))
            self.dist_sd = np.append( self.dist_sd, \
                        sd)

            self.beta = np.append( self.beta, \
                        iMerger.beta)
            self.betaPerp = np.append( self.betaPerp, \
                        iMerger.betaPerp)
            self.ClusterMass = \
              np.append( self.ClusterMass, \
                np.zeros(len(iMerger.beta))+\
                             iMerger.mass)
                             
            iMerger.getComponentMass()
            
            self.componentMass = \
              np.append( self.componentMass, \
                np.log10(iMerger.componentMass))
            


    def getMockedSample( self, **kwargs ):

        for iCluster in xrange(kwargs['nClusters']):
            iClusterSample = \
            MockedCluster( kwargs['DarkMatterStellarOffset'],  \
                    kwargs['GasStellarOffset'],
                    kwargs['DarkMatterError'], \
                    kwargs['GasError'],\
                    kwargs['nSpuriousHalos'])
            self.ClusterSample.append( iClusterSample )
        

    
            
    def getComponentMasses( self ):

        for iCluster in self.mergerHalos:
            iCluster.getComponentMass()
