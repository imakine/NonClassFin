# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:28:08 2017

@author: TH Lab
"""
import time 

def initialisation(thickness,step):
    
    '''
    This function initialize in the memory of the computer main parameters used in the code
    
    thickness : distance on the z-axis of the 1D media - integer
    
    step : the precision on the z-axis, first bin is defined between 0 and step, 
    the seconde between step and 2*step,...etc  - float
    '''

    initialtime = time.time();
    variance = [0]*(int(thickness/step));                                   #List of the variance of each element in flux_local
    std = [0]*(int(thickness/step));                                  #List of the standart deviation of each element in flux_local
    ERRORpic = [0]*(int(thickness/step));                               #List of the statistical error on the flux
    fluxpic = [0]*(int(thickness/step));
    variancepic = [0]*(int(thickness/step));                            #List of the sum of the square of each scattering
    s =[0]*6;                                                           #List of the 6 first moments of s
    flux_local = [0]*(int(thickness/step));                             #flux of neutrons = List of an estimator of the flux in the media ; collision between 0 and 1 is in first place, between 1 and 2 is second ...
    numesc=0;                                                           #Number of particle lost outside
    numdeath =0;                                                        #Number of particles absorbed in the media
    numscattot = 0;                                                     #Number of collisions for ALL particles  
    freepath =0;
    return initialtime,variance,std,ERRORpic,fluxpic,variancepic,s,flux_local,numdeath,numscattot,numesc,freepath
