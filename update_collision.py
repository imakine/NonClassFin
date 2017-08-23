# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:42:57 2017

@author: TH Lab
"""
import random

def update_collision(z,w,c,z1,step,flux_local,iflux,iscat,ideath,sigmascat,sigmatot):
    '''
    This function makes the event 'scattering' or 'absorption' in the life of the particle
    
    z : Actual position of the particle - float 
    
    w : Actual direction of the particle - float between -1 and 1
    
    c : scattering ratio - float between 0 and 1
    
    z1 : is the last known position of the particle - float 
    
    step : the precision on the z-axis, first bin is defined between 0 and step, 
    the seconde between step and 2*step,...etc  - float
    
    flux_local : Actual flux to update - List
    
    iflux : Results collected for the current particle - List
    
    iscat : number of scattering for the current particle - Integer
    
    ideath : 0 or 1, to know if the current particle has absorbed - Integer (0 or 1)
    '''
    index = int(z1/step)+1;                             #Count in the bins
    flux_local[index-1] = flux_local[index-1] + 1;      #Increment event for flux
    iflux[index-1] = iflux[index-1] +1;                 #Increment event for flux
    xi = random.random();
    if (xi < c):                                        #Check if it's a scattering with a probability of Sigmascat/sigmatot
        iscat = iscat + 1;                              #Increment of the number of scattering for the particle
        z=z1;                                           #Save the new position
        randnum=random.random()                         #Isotropic collision, direction random
        w=2*randnum-1;
    else:                                                #Check if it's an absorption with a probability Sigmaabs/sigmatot
        ideath = 1;                                      #Kill the particle
    return flux_local,iflux,iscat,z,w,ideath