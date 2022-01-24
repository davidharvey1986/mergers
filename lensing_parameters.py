import numpy as np
import cosmolopy.distance as dist

def hubble_z( z ):
    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = dist.set_omega_k_0(cosmo)

    return dist.hubble_z( z, **cosmo)
def comoving_volume( z ):

    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = dist.set_omega_k_0(cosmo)

    return dist.comoving_volume( z, **cosmo )

def diff_comoving_volume( z ):

    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = dist.set_omega_k_0(cosmo)

    return dist.diff_comoving_volume( z, **cosmo )

def comoving_distance( z ):

    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = dist.set_omega_k_0(cosmo)

    return dist.comoving_distance( z, **cosmo )

def lum_distance( z, z0=0.):
    '''
    THis is just a wrapper around angular diameter distance so
    I ndont need to constantly set the cosmo

    Returns the luminoisty diameter distance in Mpc
    '''
    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = dist.set_omega_k_0(cosmo)    
    
    return dist.luminosity_distance(z, z0=z0, **cosmo)
    
def ang_distance( z, z0=0.):
    '''
    THis is just a wrapper around angular diameter distance so
    I ndont need to constantly set the cosmo

    Returns the angular diameter distance in Mpc
    '''
    
    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = dist.set_omega_k_0(cosmo)    
    
    return dist.angular_diameter_distance(z, z0=z0, **cosmo)

def critical_kappa( z_lens=0.3, z_source=1.0 ):
    
    '''
    Determine the critical kappa for a given lens
    setup

    Return the critical kappa in units of M_sun / pc^2
    '''

    #Globals in M_sun and pc
    G = grav_const()
    c = 3e5
    
    
    

    #Various angular diameter distance    
    Dl = ang_distance(z_lens )*1e6
    Dls = ang_distance(z_source, z0=z_lens )*1e6
    Ds = ang_distance(z_source)*1e6
    
    return c**2/(4*np.pi*G)  * Ds / (Dls*Dl)  


def critical_dens( z_lens=0.3):
    '''
    PURPOSE : Function that calculates the critical
            density at a given lens                  
            Assuming a cosmology of
            omega_matter = 0.3, omega_lambda = 0.7, 
            hubble=0.7
            and negligble radiation and curvature                                          

    Arguments :                                                                    
            z_lens : redshift of the lens 

    RETURNS : Critical density in units of M_sun / pc^3

    '''
    G = grav_const() # in (km/s)^2 pc M^-1
    hubble = 70.0/1e6 # In km/s/pc
  

    hubble_at_z = hubble*np.sqrt(0.3*(1.0+z_lens)**3+0.7)
    
    critical_density = 3.0/(8.0*np.pi*G) * hubble_at_z**2
    
    return critical_density


def intrinsic_ellipticity( shear, ellipticity_disp):
    '''

    Add on the intrinsic ellipticities to the shear of galaies
    with some shear.

    shear should be a complex number
    galaxies = e1 + ie2

    
    
    if ellipticity_disp != 0.0:
        intrinsic_e = 
    else:
        intrinsic_e = np.zeros( len(shear))
        
    intrinsic_t = np.random.random( len(shear) )*np.pi
    '''
    e1_intrinsic = ell_cut( 0.0, ellipticity_disp/np.sqrt(2.), \
                                len(shear), [-1,1] )
    e2_intrinsic = ell_cut( 0.0, ellipticity_disp/np.sqrt(2.), \
                                len(shear), [-1,1] )

    intrinsic =  1.*e1_intrinsic + 1j*e2_intrinsic

    chi = source_to_image( intrinsic, shear )
      
    return chi

def image_to_source( chi, shear):
    '''

    With the complex input of a shear put the input observed
    chi and put back to the source plane

    inpiut :
       chi : Complex observed ell of the galaxy in the image plane
       shear : the reduced shear of the observed ell of the galaxy

    '''

    return ( chi-2.*shear + shear**2*np.conj(chi) ) / \
      (1.0+np.abs(shear)**2-2.0*np.real(shear*np.conj(chi)))


def source_to_image( chi, shear ):
    
    '''

    With the complex input of a shear put the input intrinsic
    chi and put forward to the image plane

    inpiut :
       chi : Complex intrinisc ell of the galaxy in the source plane
       shear : the reduced shear of  of the galaxy

    '''

    return (chi+2.*shear + shear**2*np.conj(chi)) / \
      (1.0+np.abs(shear)**2+2.0*np.real(shear*np.conj(chi)))


    
def grav_const():
    '''
    Return grav constant in Pc / M_sun (km/s)^2
    '''
    
    return 0.0042992854

def ell_cut( mean, width, ngal, range_cut ):

    if width != 0.:
        dist = np.random.normal( mean, width, ngal)
    else :
        dist = np.zeros(ngal,float)
    
    return cut( mean, width, dist, range_cut ) 

def cut( mean, width, samples, range_cut ):

    redo = samples[ (samples < range_cut[0]) |
                (samples > range_cut[1]) ]
    
    if len(redo) > 0:
        redo_dist = np.random.normal( mean, width, len(redo))
        samples[ (samples < range_cut[0]) |
                    (samples > range_cut[1]) ] = redo_dist
        cut( mean, width, samples, range_cut )

    return samples

        

        
        
def einstein_radius( M, zl, zs ):
    '''compute the einstein radius of a point mass in pc
    '''

    dl = ang_distance( zl, z0=0.)*1e3
    ds = ang_distance( zs, z0=0.)*1e3
    dls = ang_distance( zs, z0=zl)*1e3

    k = dls / ds

    G = grav_const()
    c = 3e5

    return np.sqrt( 4.*G*M/c**2*k)*dl

def einstein_radius_SIS( sigma, zl, zs ):
    '''compute the einstein radius of an SIS in kpc                         
    '''
    #critical kapp in units of M_sun / pc^2  

    k = ang_distance( zs, z0=zl) / ang_distance(zs)
    c = 3e5
    einsteinRad = 4.*np.pi*(sigma/c)**2*k

    return einsteinRad*ang_distance(zl)*1e3

def virial_radius( mass, redshift, overdensity=200 ):
    '''
    Given the virial mass, work out the radius at the overdensity
    given

    return the virial radius in kpc
    '''
    return (mass/(4./3.*np.pi* \
                critical_dens( z_lens=redshift)\
                *overdensity))**(1./3.)/1e3
