'''

A single script to get the sim data

its uses BCG_check_centre, but it has become all very
convoluted so this is a wrapper

'''
def getCrossSection( simName ):

    #cover the three CDM names'
    if 'CDM' in simName:
        simName = 'CDM'
        
    crossSections = \
        {'CDM':0, 'SIDM0.1':0.1, 'SIDM0.3':0.3, 'SIDM1':1, 'vdSIDM':1}

                     

def getSimNameList( baryonic='fiducial'):
    if baryonic == 'high':
       Sims = ['CDM_low','SIDM0.1','SIDM0.3','SIDM1','vdSIDM']
    elif baryonic == 'high':
       Sims = ['CDM_hi','SIDM0.1','SIDM0.3','SIDM1','vdSIDM']
    else:
       Sims = ['CDM','SIDM0.1','SIDM0.3','SIDM1','vdSIDM']
     
    
    return Sims


      
def getSimList( SimName ):
    
    if SimName == 'CDM_low':
        Sims = ['CDM_low']
    elif SimName == 'CDM_hi':
        Sims = ['CDM_hi']
    else:
        Sims = ['%s+baryons' % SimName]
    return Sims


