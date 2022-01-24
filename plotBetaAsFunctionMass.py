
from getMergerSample import ClusterSample
import getSimData as GSD
from matplotlib import pyplot as plt
from HSTconstraints import getMean, getMeanCauchy
from scipy.stats import norm
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import cauchy
from scipy.stats import rayleigh
from GetAndPlotTrend import *
import matplotlib as mpl
import pickle as pkl
plt.style.use('science')
from astro_tools import getColourFromRange
from lensing_parameters import virial_radius


def massMatchSample():

    SimNames = ['CDM_low','CDM_hi','SIDM0.1','SIDM0.3','SIDM1','vdSIDM'] #GSD.getSimNameList()
    ax = plt.gca()

    bins=40
    colors = ['k','b','g','r','y']
    mpl.rcParams['lines.linewidth'] = 2

    fig, axarr = plt.subplots(len(SimNames))
    
    for i, iSimName in enumerate(SimNames):

        iClusterSample = ClusterSample(iSimName)

        iClusterSample.extractMergerHalos()

        iClusterSample.CalculateOffsetVectors(nClusters=-1)

        axarr[i].hist(iClusterSample.ClusterMass, bins=40, alpha=1, label=iSimName)
        axarr[i].set_xlim(14, 15.2)
    plt.legend()

    plt.show()
    

   
def main():

    '''
    plot and devise the constraints for the HST proposal
    '''
    SimNames = ['SIDM0.1','SIDM0.3','SIDM1','vdSIDM'] #GSD.getSimNameList()
    fig, axarr = plt.subplots( 3, 1,  figsize=(5,7))
    fig.subplots_adjust(hspace=0)
    
    bins=40
    colors = ['b','g','r','y']
    mpl.rcParams['lines.linewidth'] = 2
    
    crossSections = [0.1, 0.3, 1.0, 2.0]
    plotCDM( axarr )
    allParams = []
    allParamsErr = []
    for i, iSimName in enumerate(SimNames):
            
        massBinCentres, totDistSI, totDistSG, totBeta = \
           getVectors( iSimName )


        label = iSimName.split('_')[0]
            
        axarr[0].errorbar(massBinCentres,  totDistSI[0,:], \
                           yerr=totDistSI[1:,:],capsize=2, \
                           color=colors[i], fmt='o', label=label)

        if iSimName != 'vdSIDM':
            func = fitfuncLin
            pinit = [1., -1.]
        else:
            func = fitfuncExp
            pinit = [4., 13., -10.]
        
        pfinal, pErr = getTrend(massBinCentres, totDistSI[0,:], \
                                    np.mean(totDistSI[1:,:],axis=0), \
                      func=func, pinit = pinit)
        print(pfinal)
        massPlot = np.linspace(12,14,100)
        axarr[0].plot(  massPlot, \
                      func(massPlot, *pfinal),\
                      '-', color=colors[i])

                        
        axarr[1].errorbar(massBinCentres,  totDistSG[0,:], \
                           yerr=totDistSG[1:,:], capsize=2,\
                           color=colors[i], fmt='o')
                           
        pfinal, pErr = getTrend(massBinCentres, totDistSG[0,:],\
                            np.mean(totDistSG[1:,:], axis=0), \
                      func=func,pinit = pinit)

        axarr[1].plot(  massPlot, \
                      func(massPlot, *pfinal),\
                      '-', color=colors[i])
        func = fitfuncExp
        pinit = [4., 13., -10.]
            
        axarr[2].errorbar(massBinCentres,  totBeta[0,:], \
                           yerr=totBeta[1:,:], capsize=2, \
                           color=colors[i], fmt='o')
                           
        pfinal, pErr = getTrend(massBinCentres, totBeta[0,:], \
                                    np.mean(totBeta[1:,:],axis=0), \
                      func=func, pinit = pinit)

        axarr[2].plot(  massPlot, \
                      func(massPlot, *pfinal),\
                      '-', color=colors[i])

        #if 'vd' not in iSimName:
        #    allParams.append( pfinal )
        #    allParamsErr.append( pErr)

    axarr[0].legend(loc=2)

    
    for i in range(3):
        axarr[i].set_xlim(12.25, 13.75)
        axarr[i].yaxis.set_label_coords(-0.07, 0.5)

        if i != 2:
            axarr[i].set_xticklabels([])
            axarr[i].set_xticklabels([])

    axarr[0].set_ylim(-1., 10)
    axarr[2].set_ylim(-0.05, 0.28)
    
    axarr[0].set_ylabel(r'$\bar{SD}$ / kpc')
    axarr[1].set_ylabel(r'$\bar{SG}$ / kpc')
    axarr[2].set_ylabel(r'$\beta=\bar{SD}/\bar{SG}$')
    axarr[2].set_xlabel(r'$\log[M(<100$kpc$)/M_\odot]$')
    plt.savefig('plots/constantSNR.pdf')
    plt.show()

    
    allParams = np.array(allParams)
    allParamsErr = np.array( allParamsErr)
    crossSections = np.array(crossSections)
    pfinalA, pErrA = getTrend(crossSections[:-1], allParams[:,0],\
                            allParamsErr[:,0], \
                      func=fitfuncLin)
    print("A = %0.3f + %0.3f  sigma" % tuple(pfinalA))

    pfinal, pErr = getTrend(crossSections[:-1], allParams[:,1],\
                            allParamsErr[:,1], \
                      func=fitfuncLin)
    print("B = %0.3f pm %0.3f s" % tuple(pfinal))

    fig, axarr = plt.subplots(2, 1)
    axarr[0].errorbar(crossSections[:-1], allParams[:,0], allParamsErr[:,0])
    axarr[0].plot( crossSections[:-1], fitfuncLin( crossSections[:-1], *pfinalA))
    axarr[1].errorbar(crossSections[:-1], allParams[:,1], allParamsErr[:,0])
    axarr[1].plot( crossSections[:-1], fitfuncLin( crossSections[:-1], *pfinal))
    plt.show()
    
def getVectors(iSimName):

    iClusterSample = ClusterSample(iSimName)

    iClusterSample.extractMergerHalos()

    iClusterSample.CalculateOffsetVectors(nClusters=-1)

    nClustersInSim = len(iClusterSample.componentMass)

    nMassBins = 9
    nMassPerBins = nClustersInSim // nMassBins
        
    massSortedIndex = np.argsort( iClusterSample.componentMass )

    

    massBinCentres = []
    totBeta = np.zeros( (3,  nMassBins))
    totDistSI = np.zeros( (3,   nMassBins))
    totDistSG = np.zeros( (3,  nMassBins))   
    
    stdTotDistSG = np.zeros( (2,nMassBins))
    crossSection = np.zeros(  nMassBins)
                      
    for iMassBin in range(nMassBins):
        if iMassBin == nMassBins -1:
            endInd = nClustersInSim
        else:
            endInd = np.min([(iMassBin+1)*nMassPerBins, nClustersInSim])
        
        index = massSortedIndex[ iMassBin*nMassPerBins:endInd ]
        #index = allIndexes[ iClusterSample.dist_sg[allIndexes] < 250]
        massBinCentres.append( np.mean(iClusterSample.componentMass[index]))
        
        nClusters = len(index)

        if nClusters == 0:
            print(iMassBin*nMassPerBins,endInd)
            continue
        print("selected %f/%f clusters" % \
                  (nClusters,len(iClusterSample.componentMass)))

        SG = getMean(iClusterSample.dist_sg[index] )
        SI = getMean(iClusterSample.dist_si[index])
        DI = getMean(iClusterSample.dist_di[index])
        beta = getMean(iClusterSample.beta[index])
        

        totDistSI[ :, iMassBin] = SI
        totDistSG[:, iMassBin] = SG
        totBeta[ :, iMassBin] = beta
        
    massBinCentres = np.array(massBinCentres)

    return massBinCentres, totDistSI, totDistSG, totBeta


def plotCDM( axarr ):

    simNames = ['CDM']
    massPlot = np.linspace(12,14,100)
    allValues = np.zeros( (len(simNames), len(axarr), len(massPlot), 2))
    for iSim, iSimName in enumerate(simNames):
        allVectors = getVectors(iSimName)
        
        for iAxis, iVector in enumerate(allVectors[1:]):
                           
            pfinal, pErr = getTrend(allVectors[0], iVector[0,:],\
                            np.mean(iVector[1:,:], axis=0), \
                        func=fitfuncLin,pinit =  [1., -1.])

            allValues[ iSim, iAxis, :, 0] = \
              fitfuncLin( massPlot, *[pfinal[0]+pErr[0], pfinal[1]-pErr[1]])

            allValues[ iSim, iAxis, :, 1] = \
              fitfuncLin( massPlot, *[pfinal[0]-pErr[0], pfinal[1]+pErr[1]])
              
    for iA, iAxis in enumerate(axarr):

        lower = np.min(np.min( allValues[ :, iA, :, :], axis=0), axis=1)
        upper = np.max(np.max( allValues[ :, iA, :, :], axis=0),axis=1)
        iAxis.fill_between( massPlot, lower, upper, label='CDM', \
                             facecolor='grey', hatch='//', alpha=0.2)
        iAxis.plot(  massPlot, lower, 'k-')
        iAxis.plot(  massPlot, upper, 'k-')
        

    

        
        
if __name__ == '__main__':
    main()
