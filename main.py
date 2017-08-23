# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:44:29 2017

@author: TH Lab
"""
from Initialisation import initialisation
from initializeOneParticle import initializeOneParticle
from lifeOfParticle import lifeOfParticle
from addParticleInSample import addParticleInSample
from averageSample import averageSample
import time
from output import output
from theoretical import theoretical
from initCrossSec import initCrossSec


def main(c,thickness,n):
    '''
    This function gives the flux for SP1 or SP2 or SP3 with Monte-Carlo Simulation
    c is the scattering ratio (ratio between the scattering cross section and the total cross-section)
    '''
    
    #Initialisation of the parameters of the simulation
    
    #thickness = 100;                     #Thickness of the mediq
    #m1 = 2.0004668136;
    #m2 = 8.79150636775;
    #m3 = 56.1525616557;
    #m4 = 457.65052155;
    #m5 = 4521.87604098;
    #m6 = 52969.8811984;
    m1=1
    m2=2
    m3=6
    m4=24
    m5=120
    m6=720
    #sigmatot = 1;                       #Total cross section
    #sigmaabs = (1-c)*sigmatot;          #Absorption cross-section of the media
    #sigmascat = c*sigmatot;             #Scattering cross-section of the media
    region = 1                      # If True the source is regional source on a distance sourceregion about the half of the media, if False, the source is uniform
    sourceregion = 1.0                  #thickness of center source region    
    Q = 1                              #Source of neutrons
    #n = 10000000;                       #Number of path
    step = 0.1;                         #Step for bins Number of bins will be Length/step
    sp = 1                             #SP =1,2,3 
    sigt1 = 1;
    l1 = 0.5;
    l2 = 0.5;
    sigmoy = sigt1*l1/(l1+l2)
    sigmatot,sigmaabs,lambda1,beta1,lambda2,beta2 = initCrossSec(m1,m2,m3,m4,m5,m6,sp,c)
    sigmascat = c*sigmatot;
    

    # Initialisation in the Memory.
    
    initialtime,variance,std,ERRORpic,fluxpic,variancepic,s,flux_local,numdeath,numscattot,numesc,freepath = initialisation(thickness,step)

    # Loop on each particle of the sample
    
    for i in range(1,n+1):
        
        # Interface with user
        
        if (i%10000 == 0):
            print(str(i/n*100) + "% of the running code done")  
            
        # Initialisation of the position, direction, etc. of the particle i
        
        z,w,iesc,ideath,iscat,iflux = initializeOneParticle(thickness,sourceregion,step,region)
        
        # Following of the life of the particle i
        
        s,iesc,ideath,iscat,flux_local,iflux=lifeOfParticle(sp,iesc,ideath,sigmatot,freepath,s,z,w,c,step,flux_local,iflux,iscat,thickness,m2,lambda1,beta1,lambda2,beta2,sigmascat)
        
        #Computation of interest elements for the sample
        
        fluxpic,variancepic,numesc,numdeath,numscattot = addParticleInSample(flux_local,fluxpic,iflux,variancepic,numesc,iesc,numdeath,numscattot,iscat,ideath)
    fluxpic,variancepic,std,ERRORpic,flux_local,s = averageSample(m1,m2,fluxpic,variancepic,std,ERRORpic,n,Q,step,sigmoy,flux_local,s,variance,numscattot,numdeath,numesc,region,thickness)
    finaltime = time.time(); 
    
    # Computation of the different theoretical moments to obtain
    
    s1,s2,s3,s4,s5,s6 = theoretical(sp,sigmatot)
    
    duration = finaltime-initialtime
    
    # Generate all the output files
    
    output(sp,c,s1,s2,s3,s4,s5,s6,s,finaltime,initialtime,flux_local,thickness,sigmaabs,sigmascat,sigmatot,Q,n,step,numscattot,numesc,numdeath,ERRORpic,std,variancepic)
    return flux_local,ERRORpic,std,duration

c = list(reversed([0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.999]))

l =  list(reversed([50,50,50,80,80,80,80,80,100,200,200,500,1500]))

n =  list(reversed([20000000,20000000,20000000,20000000,20000000,20000000,20000000,20000000,20000000,20000000,2000000,2000000,200000]))


for i in range(len(c)):
    a,error,std,delay = main(c[i],l[i],n[i])
    filename = str("TEST")
    file = open(filename + str(c[i]) + '_CRsp32' +'.txt', 'w')
    file.write('max of the flux = ' +  str(max(a)) + '\n')
    file.write("average on bins max = "+  str((a[a.index(max(a))]+a[a.index(max(a))+1])/2)+'\n')
    file.write("Statistical error on max = " + str(error[int(l[i]*10/2)]) + " and max error = " + str(max(error)) + 'n =' + str(n[i]) + '\n')
    file.write("Standart deviation on max = " + str(std[int(l[i]*10/2)])+'\n')
    file.write("Duration of computation =" + str(delay))
    file.close()