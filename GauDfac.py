'''
'''
def GauDfac(bmaj1,bmin1,bpa1,bmaj2,bmin2,bpa2):
    
    r'''
    Determine the parameters of a gaussian deconvolved with another gaussian.
    
    Parameters
    ----------
    bmaj1 : float
        Major axis (FWHM) of the source
    bmin1 : float
        Minor axis (FWHM) of the source
    bpa1 : float
        Major axis position angle of the source, in degrees
    bmaj2 : float
        Major axis (FWHM) of the beam
    bmin2 : float
        Minor axis (FWHM) of the beam
    bpa2 : float
        Major axis position angle of the beam, in degrees
        
    Returns
    -------
    bmaj : float
        Major axis of resultant gaussian    
    bmin : float
        Minor axis of resultant gaussian 
    bpa : float
        Major axis position angle of resultant gaussian, in degrees
    ifail : int  
        Success status: 0   All OK.
                        1   Result is pretty close to a point source.
                        2   Illegal result.  
    
    Notes
    -----
    Required packages: numpy
    Author: R. C. Levy
    Based on GauDfac.m by A. D. Bolatto, which is based on Miriad
    Last updated: 2018-07-25
    
    Examples
    --------
    >>> import GauDfac
    >>> bmaj1 = 10.0
    >>> bmin1 = 8.0
    >>> bpa1 = 45.0
    >>> bmaj2 = 1.5
    >>> bmin2 = 1.2
    >>> bpa2 = 30.3
    >>> bmaj,bmin,bpa,flag = GauDfac.GauDfac(bmaj1, bmin1, bpa1, bmaj2, bmin2 , bpa2)
    >>> print(bmaj,bmin,bpa,flag,sep=', ')
    9.889554000017112, 7.906119255408783, 45.3227385587466, 0
    
    
    '''
    
    import numpy as np
    D2R=np.pi/180.
    R2D=1./D2R
    theta1 = bpa1 * D2R
    theta2 = bpa2 * D2R
    
    alpha  = (bmaj1*np.cos(theta1))**2 + (bmin1*np.sin(theta1))**2 - (bmaj2*np.cos(theta2))**2 - (bmin2*np.sin(theta2))**2
    beta   = (bmaj1*np.sin(theta1))**2 + (bmin1*np.cos(theta1))**2 - (bmaj2*np.sin(theta2))**2 - (bmin2*np.cos(theta2))**2
    gamma  = 2 * ( (bmin1**2-bmaj1**2)*np.sin(theta1)*np.cos(theta1) - (bmin2**2-bmaj2**2)*np.sin(theta2)*np.cos(theta2) )
    
    s = alpha + beta
    t = np.sqrt((alpha-beta)**2 + gamma**2)
    limit = np.min([bmaj1,bmin1,bmaj2,bmin2])
    limit = 0.1*limit*limit
    if alpha<0 or beta<0 or s<t:
        bmaj = 0.
        bmin = 0.
        bpa = 0.
        if 0.5*(s-t)<limit and alpha>-limit and beta>-limit:
            ifail = 1
        else:
            ifail = 2
    else:
        bmaj = np.sqrt(0.5*(s+t))
        bmin = np.sqrt(0.5*(s-t))
        if np.abs(gamma)+np.abs(alpha-beta)==0:
            bpa = 0.0
        else:
            bpa = 0.5 * np.arctan2(-gamma,alpha-beta) * R2D
        
        ifail = 0
        
    return bmaj,bmin,bpa,ifail        