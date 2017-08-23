# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 10:30:39 2017

@author: TH Lab
"""

def initCrossSec(m1,m2,m3,m4,m5,m6,sp,c,lambda1=None,beta1=None,lambda2=None,beta2=None):
    if sp ==1:
        sigmatot = 2*m1/m2
        sigmaabs = (1-c)/m1
        #q = Q
    if sp==2:
        lambda1 = 3*m4/(10*m2**2)-m3/(3*m1*m2)
        beta1=m3/(3*m1*m2)-1
        sigmatot = 2*m1/m2
        sigmaabs = ((1-c)/m1)*(1-beta1*(1-c))/(1+lambda1*(1-c))
        #q = Q*(1-c)/m1*(1-beta1*(1-c))/(1+lambda1(1-c))
    if sp==3:
        lambda1 = 3*m4/(10*m2**2)-m3/(3*m1*m2)
        beta1=m3/(3*m1*m2)-1 
        lambda2 =(9*m5/5-27*m1*m6/(21*m2)+3*m3*m4/m2-10*m3**2/(3*m1))/(10*m2*m3-9*m1*m4)
        beta2 = (10*m3**2/(3*m1)-9*m5/5)/(10*m2*m3-9*m1*m4)-1
        sigmatot = 2*m1/m2
        sigmaabs = ((1-c)/m1)*(1-beta1*(1-c))/(1+lambda1*(1-c))
        #q =Q /(1+beta1*(1-c))
    return sigmatot,sigmaabs,lambda1,beta1,lambda2,beta2