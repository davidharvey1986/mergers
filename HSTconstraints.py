from getMergerSample import *
import getSimData as GSD
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from scipy.stats import norm
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import cauchy
from scipy.stats import rayleigh
from GetAndPlotTrend import *
import matplotlib as mpl

def HSTconstraints():

    '''
    plot and devise the constraints for the HST proposal
    '''

    SimNames = GSD.getSimNameList()
    fig, axarr = plt.subplots( 1 )

    bins=40
    colors = ['k','r','b']
    mpl.rcParams['lines.linewidth'] = 2
    nCross = len(SimNames)
    totBeta = np.zeros( (3, nCross))
    totDistSI = np.zeros( (3, nCross))
    totDistSG = np.zeros( (3, nCross))

    crossSection = np.array([ GSD.getCrossSection(i) for i in SimNames])

    for i, iSimName in enumerate(SimNames):

        iClusterSample = ClusterSample(iSimName)

        iClusterSample.extractMergerHalos()
        
        iClusterSample.CalculateOffsetVectors(nClusters=-1)

        #I use this to cut in different options, for example
        # - unphysical SI
        # - mass of cluter
        # - mass of component
        
        index = iClusterSample.dist_si < 250
        nClusters = len(index[index])

        print("selected %f/%f clusters" % \
                  (nClusters,len(iClusterSample.ClusterMass)))

        
        #plotBulletVectors( iClusterSample )

        
        SG = getMean(iClusterSample.dist_sg[index] )
        SI = getMean(iClusterSample.dist_si[index])
        DI = getMean(iClusterSample.dist_di[index])
        beta = getMean(iClusterSample.beta[index])
        
        totDistSI[:, i] = SI
        totDistSG[:, i] = SG
        totBeta[:, i] = beta        

        
        print("%s SG mean is %f +/- %0.2f" %  \
              ( iSimName, SG[0], np.mean(SG[1:])))
        print("%s SI mean is %f +/- %0.2f" %
              ( iSimName, SI[0], np.mean(SI[1:])))
        print("%s DI mean is %f +/- %0.2f" %  \
              ( iSimName, DI[0], np.mean(DI[1:])))
        print("%s Beta mean is %f +/- %0.2f" %  \
              ( iSimName, beta[0], np.mean(beta[1:])))


    #Dump this info in a pickle so i can use it for other scripts
    pkl.dump([crossSection,totBeta[0,:], totBeta[1,:], totBeta[2,:],\
              totDistSI[0,:], totDistSG[0,:]],\
                 open('Pickles/SimulationBetaToCross.pkl','wb'))

    plt.sca(axarr)
    ax = plt.gca()

    BetaCrossSectionTrend, pEr = \
       getTrend(crossSection,  totBeta[0, :], \
                                np.min(totBeta[:-1, :],axis=0), \
                                func=fitfunc)

    plotCross = np.linspace(0, 2, 100)
    ax.plot( plotCross, fitfunc( plotCross, *BetaCrossSectionTrend), '--')
    
    ax.errorbar( crossSection, totBeta[0,:], yerr=totBeta[1:,:], color='black', fmt='o', \
                           label='BAHAMAS-SIDM hydro simulations' )
    
    ax.errorbar( plotCross, 1.-np.exp(-1.*plotCross/6.5), ls='-', \
                       label='Approximate H14 analytical model', color='k')

    print("$A = %0.3f \pm %0.3f, B = %0.3f \pm %0.3f$" % \
              (BetaCrossSectionTrend[0], pEr[0], BetaCrossSectionTrend[1], pEr[1]))
              
    HarveyMeasurement = np.zeros(2)+0.07
    HarveyCrossSection =   betaToCross( HarveyMeasurement, 1., 6.5)
    BAHAMASCrossSection =  betaToCross( HarveyMeasurement, \
                                            *BetaCrossSectionTrend)



    #plotJauzacMeasurement( BetaCrossSectionTrend )
    
    #ax.errorbar( [1e-3, BAHAMASCrossSection[0]], \
    #                       HarveyMeasurement, color='red')
    #ax.text( 0.012, HarveyMeasurement[0]*1.1, \
    #    'Existing HST+Chandra archive of 30 clusters', color='red')
    #ax.errorbar(  HarveyCrossSection, \
    #            [-0.1,HarveyMeasurement[0]], color='red', ls='--', \
    #                        label='Published particle constraints')
    #ax.errorbar(  BAHAMASCrossSection, \
    #            [-0.1,HarveyMeasurement[0]], color='red',ls=':',\
    #                    label='Revised interpretation (cf.W17)')


    #plotSensitivity(SI,meanSG,BetaCrossSectionTrend)
    

    #ax.set_xscale('log')
    #ax.set_yscale('log')

    #axarr[0].legend(loc=0)
    ax.legend(loc=0)

    #axarr[0].set_xlabel('Star - Gas Distance / kpc')
    #axarr[0].set_ylabel('Star - Intersection / kpc')
    ax.set_xlabel(r'Dark Matter Cross-Section / cm$^2$/g', fontsize=12)
    ax.set_ylabel(r'Observable $\beta=\delta_{\rm SI}/\delta_{\rm SG}$', fontsize=12)
    ax.set_ylim(-1e-2,0.2)
    ax.set_xlim(-5e-2,1.75)
    plt.savefig('plots/HSTconstraints.pdf')
    plt.show()

def getMaxLike( y, bins=30 ):

    #mean = np.median(x)
    #std = np.std(x)/np.sqrt(len(x))

    ypdf, x = np.histogram( y, bins=bins, density=True)
    xc = (x[:-1]+x[1:])/2.
    maxLike = xc[ np.argmax(ypdf) ]
    deltax = xc[1] - xc[0]
    cumsum = np.cumsum(ypdf*deltax)
    std = x[np.argmin(np.abs(cumsum - 0.84))] - maxLike
    
    return maxLike, std/np.sqrt(len(y))


def getMean( x, ax=None ):

    x = np.sort(x)
    mean, std = norm.fit( x )

    low, med, hi = np.percentile(x, [16, 50, 84])
    std = (hi-low)/2.

    return med, (med-low)/np.sqrt(len(x)),  (hi-med)/np.sqrt(len(x))
    

def getMeanCauchy( x ):
    #Because a gaussian divided by a gaussian is a cauchy
    #therefor beta should be a cauchy distribution

    mean, std = cauchy.fit( x )
    
    
    return mean, std/np.sqrt(len(x))


def plot( x, y, label, c):
    nbins = 7
    bins = np.linspace(0,100,nbins+1)

    meany = np.zeros( nbins )
    std = np.zeros( nbins )

    for i in xrange(nbins):

        inbin = ( x > bins[i] ) & (x <= bins[i+1])

        meany[i] = np.mean( y[inbin] )
        std[i] = np.std( y[inbin] ) / np.sqrt(len(y[inbin]))

    bincentres = (bins[:-1] + bins[1:])/2.
    bincentres = bincentres[np.isfinite(std)]
    meany = meany[np.isfinite(std)]
    ax = plt.gca()
    std[:] = np.mean(std)
    plotTrend( bincentres, meany, None, label, color=c )
    ax.errorbar( bincentres, meany, yerr=std, color=c, label=label, fmt='o')
        






def plotSensitivity(si,sg,BetaCrossSectionTrend):
    '''
    Plot the sensitivity in beta
    '''
    error =60.
    sigmaBeta = 0.07
    nClusterBenchMark = 72
    nNewClusters = 72
    totalClusters = nClusterBenchMark + nNewClusters

    ax = plt.gca()
    plotCrossSections = np.linspace(0.001,2.0,2.)
    SensitivityBetaSingleBand = np.zeros(len(plotCrossSections))+sigmaBeta
    SensitivityBetaDualBand = np.zeros(len(plotCrossSections))+sigmaBeta/2.

    SensitivityBetaSingleBand *= np.sqrt(nClusterBenchMark)
    SensitivityBetaDualBand *= np.sqrt(nClusterBenchMark)
    yerr = np.zeros(len(plotCrossSections))+0.01



    BetaLimit = SensitivityBetaSingleBand[0]/np.sqrt(totalClusters)
    CrossSectionLimit =  betaToCross( BetaLimit, *BetaCrossSectionTrend)

    
    #ax.errorbar( [plotCrossSections[0],CrossSectionLimit], \
    #             [BetaLimit,BetaLimit], ls='-', \
    #             color='orange')
    #ax.text( 0.006, BetaLimit, \
    #            '85 Clusters + single band', color='orange')

    #ax.errorbar( [CrossSectionLimit,CrossSectionLimit],\
     #            [-1,BetaLimit], ls='-', \
     #                color='orange')
    
    
    BetaLimit = SensitivityBetaDualBand[0]/np.sqrt(totalClusters)
    CrossSectionLimit =  betaToCross( BetaLimit, *BetaCrossSectionTrend)
      
    ax.errorbar( [plotCrossSections[0],CrossSectionLimit], \
                 [BetaLimit,BetaLimit], ls='-', \
                color='green')
    ax.errorbar( [CrossSectionLimit,CrossSectionLimit], \
                 [1e-4,BetaLimit], ls='-', \
                color='green')
                
    ax.text( 0.012, BetaLimit*1.1, \
        '$\sim10^2$ Clusters Euclid', color='green')

    print('Cross section limit is ',CrossSectionLimit)
    ax.errorbar( [CrossSectionLimit,CrossSectionLimit],\
                 [-1,BetaLimit], ls='-', \
                     color='green')


def plotJauzacMeasurement(BetaCrossSectionTrend):

    dm = 0.5
    dmerr = 0.5
    gas = 1.7
    gaserr = 0.5


    JauzacMeasurement = dm/gas
    JauzacMeasurementError = np.sqrt((dmerr/dm)**2+(gaserr/gas)**2)*JauzacMeasurement
    JauzacCrossSection = (JauzacMeasurement - \
      BetaCrossSectionTrend[0])/BetaCrossSectionTrend[1]
    JauzacCrossSectionError = \
      JauzacMeasurementError/JauzacMeasurement*JauzacCrossSection

    ax = plt.gca()

    ax.plot( [1e-3, JauzacCrossSection], \
                 [JauzacMeasurement,JauzacMeasurement], \
                 '-', color='cyan', label='J18')

    ax.fill_between( [1e-3, 10], \
                    [JauzacMeasurement-JauzacMeasurementError,\
                         JauzacMeasurement-JauzacMeasurementError], \
                    [JauzacMeasurement+JauzacMeasurementError,\
                         JauzacMeasurement+JauzacMeasurementError], \
                         color='cyan', alpha=0.5)
                         
    ax.fill_between( [JauzacCrossSection-JauzacCrossSectionError, \
                          JauzacCrossSection+JauzacCrossSectionError], \
                    [JauzacMeasurement,JauzacMeasurement], \
                    [-1,-1], \
                         color='cyan', alpha=0.5)
    print("JAUZAC CONSTRAINTS ARE %0.2f +/- %0.2f" \
              %(JauzacCrossSection,JauzacCrossSectionError))

    x=np.linspace(-10,30,1000)
    gaussCDF = np.cumsum(norm.pdf(x, JauzacCrossSection, JauzacCrossSectionError))
    gaussCDF/=np.max(gaussCDF)
    JauzacCrossSectionLimit = x[np.argmin(np.abs(gaussCDF-0.68))]
    print("JAUZAC LIMITS IS <  %0.2f " \
              %(JauzacCrossSectionLimit))

    
    
def betaToCross( beta, Aparam, sigmaStar):
    return np.log(1.-beta/Aparam)*-sigmaStar


def plotBulletVectors( iClusterSample ):
     plt.sca( axarr[0] )
     plot( iClusterSample.dist_sg[index], iClusterSample.dist_si[index], \
                  iSimName, colors[i])

     axarr.scatter( iClusterSample.dist_sg[index], iClusterSample.dist_si[index], \
                          color=colors[i])

     axarr.hist(  iClusterSample.dist_sd[index], bins=30, alpha=1-i*0.3, \
                           density=True, color=colors[i] )
     axarr[1].hist(  iClusterSample.dist_sg[index], bins=30, alpha=1-i*0.3,\
                            density=True)
     axarr[2].hist(  iClusterSample.beta[index], bins=30, alpha=1-i*0.3,\
                            density=True)
if __name__ == '__main__':
    HSTconstraints()
