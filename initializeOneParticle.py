# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:29:58 2017

@author: TH Lab
"""
from PosAndW import posInit,muInit

def initializeOneParticle(thickness,sourceregion,step,region):
    '''
    This function initialize each particle of history
    
    thickness : distance on the z-axis of the 1D media - Integer
    
    sourceregion: the distance where we have a constant source - Float
    
    region: True if regional source False if uniform source - Boolean
    
    step : the precision on the z-axis, first bin is defined between 0 and step, 
    the seconde between step and 2*step,...etc  - float
    '''
    z0 = posInit(thickness,sourceregion,region);
    w0 = muInit();
    z=z0;                                                           # position z follows the particle
    w=w0;                                                           # direction w follows the particle                                                       
    iesc=0;                                                         #Event escape value 0 --> particle does not escape, value 1 --> particle has escaped 
    ideath=0;                                                       #Event absorption by the media value 0 --> particle is not absorpbed, value 1 --> particle has been absorbed
    iscat=0;                                                        #Number of collision by particle
    iflux = [0]*(int(thickness/step));                              #Number of collision by particle
    return z,w,iesc,ideath,iscat,iflux 