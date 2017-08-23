# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:41:58 2017

@author: TH Lab
"""

def addParticleInSample(flux_local,fluxpic,iflux,variancepic,numesc,iesc,numdeath,numscattot,iscat,ideath):
    '''
    This function add the different results for the particle in the sample
    
    flux_local : Actual flux to update - List
    
    fluxpic : same than flux_local but just to compute variance - List
    
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
        fluxpic[i] = fluxpic[i] + iflux[i];                         #We add the collision of each particle for flux
        variancepic[i] = variancepic[i] + iflux[i]**2;       #We add the square power of each particle for variance
    numesc=numesc+iesc;                                             #Update of the number of lost outside
    numdeath=numdeath+ideath;                                       #Update of the number of qbsorption
    numscattot = numscattot + iscat;                                #Update of the number of scattering
    return fluxpic,variancepic,numesc,numdeath,numscattot 