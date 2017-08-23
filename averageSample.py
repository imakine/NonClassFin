# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:36:40 2017

@author: TH Lab
"""

def averageSample(m1,m2,fluxpic,variancepic,std,ERRORpic,n,Q,step,sigmoy,flux_local,s,variance,numscattot,numdeath,numesc,region,thickness):
    '''
    Knowing the contribution of the n particles, a average can be done on that. This function realized that.
    
    flux_local : Actual flux to update - List
    
    fluxpic : same than flux_local but just to compute variance - List
    
    std : Actual standart deviation of the flux - List
    
    ERRORpic : Actual statistical error of the flux - List
    
    n : number of particles that we have simulated - Integer
    
    Q : Intensity of the source - Float
    
    step : the precision on the z-axis, first bin is defined between 0 and step, 
    the seconde between step and 2*step,...etc  - float
    
    sigmatot : total cross section of the media - float
    
    s : 6 first moments of the sampling - List
    
    variance : Final variance
    
    iflux : Results collected for the current particle - List
    
    variancepic : Actual variance to update - List
    
    numesc : Actual number of escaped particles to update - Integer
    
    iesc : 0 or 1, to know if the current particle has escaped - Integer (0 or 1)
    
    numdeath : Actual number of absorbed particles to update - Integer
    
    numscattot ; Actual number of scattering - Integer
    
    iscat : number of scattering for the current particle - Integer
    
    ideath : 0 or 1, to know if the current particle has absorbed - Integer (0 or 1)
    
    Output : same but updated
    '''
    
    for i in range(len(fluxpic)):
        variancepic[i] = variancepic[i]/n;                          #Average number of collision square on the sample
        fluxpic[i] = fluxpic[i]/n;                                      #Average number of collision  on the sample
        variance[i] = variancepic[i]-(fluxpic[i])**2;                       #Variance
        std[i] = variance[i]**(0.5);                                      #Standart deviation
        if (fluxpic[i]!=0.):
            ERRORpic[i] = 2*std[i]/(fluxpic[i]*n**(0.5))              #STATISTICAL RELATIVE ERROR ON THE AVERAGE FLUX WITH 95% CONFIDENCE
    if region:
        flux_local[:] = [(x/n)*(Q)/step/sigmoy for x in flux_local];      #Average flux with good normalisation The integral of the source appears Q*1 (because Q is constant in a distance = 1)
    else:
        flux_local[:] = [(x/n)*(Q)/step/sigmoy*thickness for x in flux_local];
    s[:] = [x / (numscattot+numdeath+numesc) for x in s];       #Average of the moments
    return fluxpic,variancepic,std,ERRORpic,flux_local,s