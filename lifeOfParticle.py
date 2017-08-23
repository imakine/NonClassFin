# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:38:42 2017

@author: TH Lab
"""
import random
from moment import moment
from check_out import check_out
from update_collision import update_collision
import scipy.optimize
import numpy as np

def lifeOfParticle(sp,iesc,ideath,sigmatot,freepath,s,z,w,c,step,flux_local,iflux,iscat,thickness,m2,lambda1,beta1,lambda2,beta2,sigmascat):
    '''
    This function follows the life of the particle, from birth to escape or absorption.
    
    sp : Choice of the user for SP1 or SP2 or SP3 - Integer 1 or 2 or 3
    
    flux_local : Actual flux to update - List
    
    iflux : Results collected for the current particle - List
    
    iscat : number of scattering for the current particle - Integer
    
    iesc : 0 or 1, to know if the current particle has escaped - Integer (0 or 1)
    
    ideath : 0 or 1, to know if the current particle has absorbed - Integer (0 or 1)
    
    sigmatot : total cross section of the media - float
    
    freepath : Actual distance to the next event for the particle - float
    
    s : 6 first moments of the sampling - List
    
    z : Actual position of the particle - float 
    
    w : Actual direction of the particle - float between -1 and 1
    
    c : scattering ratio - float between 0 and 1
    
    step : the precision on the z-axis, first bin is defined between 0 and step, 
    the seconde between step and 2*step,...etc  - float
    
    thickness : distance on the z-axis of the 1D media - integer
    
    '''
    
    
    
    while((iesc+ideath)==0):
        randnum=random.random();
        if (sigmatot !=0.):
            if (sp==1):
                #freepath = scipy.optimize.brentq(fsp,0.0,100.0,args=(randnum))/((3**0.5)*sigmatot) #Sampling of the distance travelled
                freepath = scipy.optimize.brentq(fsp,0.0,100.0,args=(randnum))/((6/m2)**0.5)
            elif (sp==2):
                if randnum < lambda1/(1+lambda1) : #4/9
                    freepath = 0;
                else:
                    
                    freepath = scipy.optimize.brentq(fsp,0.0,1000.0,args=((lambda1+1)*(1-randnum)))*((1+lambda1)*m2/(6*(1-beta1*(1-c))))**0.5
            elif (sp==3):
                try:
                    freepath = scipy.optimize.brentq(f,0,100,args=((randnum,m2,beta1,beta2,lambda1,lambda2,c)))/sigmatot
                    #freepath = scipy.optimize.brentq(g,0.0,100,args=((randnum)))/sigmatot
                except:
                    print('module scipy non efficient level 1')
                    #randnum=random.random();
                    freepath = secante(0,1,1e-6,randnum,sigmatot,m2,beta1,beta2,lambda1,lambda2,c)/sigmatot
                    #freepath = secante(0,100,1e-6,randnum,sigmatot)
            s = moment(s,freepath);
        z1=z+w*freepath;                                            #Update of the position
        iesc =check_out(z1,thickness,sigmatot);
        if (not iesc):
            flux_local,iflux,iscat,z,w,ideath = update_collision(z,w,c,z1,step,flux_local,iflux,iscat,ideath,sigmascat,sigmatot)
    return s,iesc,ideath,iscat,flux_local,iflux 

def fsp(x,xi):                          # function xi - (1+z)exp(-z)
    return  xi - (1+x)*np.exp(-x)

def f(x,xi,m2,beta1,beta2,lambda1,lambda2,c):
    z1 = 1+beta1*(1-c)
    z2 = 1-beta2*(1-c)
    alphap=((1/3*(lambda2+z1*z2)+(1/9*(lambda2+z1*z2)**2-4/9*(lambda2*z1-lambda1)*z2)**(0.5))*9/(lambda2*z1-lambda1)/m2)**(0.5)
    alpham=((1/3*(lambda2+z1*z2)-(1/9*(lambda2+z1*z2)**2-4/9*(lambda2*z1-lambda1)*z2)**(0.5))*9/(lambda2*z1-lambda1)/m2)**(0.5)
    a_plus = (0.5)*lambda1/z1/(z2-3/2/z1*(1/3*(lambda2+z1*z2)+(1/9*(lambda2+z1*z2)**2-4/9*(lambda2*z1-lambda1)*z2)**(0.5)))
    a_min = (0.5)*lambda1/z1/(z2-3/2/z1*(1/3*(lambda2+z1*z2)-(1/9*(lambda2+z1*z2)**2-4/9*(lambda2*z1-lambda1)*z2)**(0.5)))
    b = lambda1/(lambda2*z1-lambda1)
    eta = (2*a_min*(1+b)+b)/(a_min-a_plus)*3    
    mu = (2*a_plus*(1+b)+b)/(a_plus-a_min)*3
    y = mu/(alpham)**2/2
    z = eta/(alphap)**2/2
    return (xi - y*(1-fsp2(alpham*x)) - z*(1-fsp2(alphap*x)))

def secante(a,b,prec,xi,sigmatot,m2,beta1,beta2,lambda1,lambda2,c):
    while f(a*sigmatot,xi,m2,beta1,beta2,lambda1,lambda2,c) > prec:
        x = f(a*sigmatot,xi,m2,beta1,beta2,lambda1,lambda2,c)
        y = f(b*sigmatot,xi,m2,beta1,beta2,lambda1,lambda2,c)
        print(b-a)
        a = a-x*(b-a)/(y-x)
        print(a)
    return a 

def secante2(a,b,prec,xi,sigmatot):
	while g(a*sigmatot,xi) > prec:
		a = a-g(a*sigmatot,xi)*(b-a)/(g(b*sigmatot,xi)-g(a*sigmatot,xi))
	return a

def g(x,xi):
    a = 5.642025;
    b = 0.469086;
    l = (5+2*(10/3)**(0.5))**(0.5);
    m = (5-2*(10/3)**(0.5))**(0.5);
    y = a/(l**2);
    z = b/(m**2);
    return (xi - y*(1-fsp2(l*x)) - z*(1-fsp2(m*x)))

def fsp2(x): # function xi - (1+z)exp(-z)
    return (1+x)*np.exp(-x)