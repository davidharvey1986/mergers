'''

A single script to get the sim data

its uses BCG_check_centre, but it has become all very
convoluted so this is a wrapper

'''


def getSimNameList():
    Sims = ['CDM','SIDM0.1','SIDM1', 'CDMhires','SIDMhires','obs']
    
    return Sims


      
def getSimList( SimName ):
    
    if SimName == 'CDM_low':
        Sims = ['CDM_low']
    elif SimName == 'CDM_hi':
        Sims = ['CDM_hi']
    else:
        Sims = ['%s+baryons' % SimName]
    return Sims


