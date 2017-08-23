# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:43:55 2017

@author: TH Lab
"""

def moment(s,freepath):
    '''
    This function gives the different moment for s (until sixth)
    
    s : 6 first moments of the sampling - List
    
    freepath : Actual distance to the next event for the particle - float
    
    '''
    
    s[0] = s[0] + freepath;                                 #First moment estimator
    s[1] = s[1] + freepath**2;                              #Second moment estimator
    s[2] = s[2] + freepath**3;                              #Third moment estimator
    s[3] = s[3] + freepath**4;                              #..
    s[4] = s[4] + freepath**5;                              #..
    s[5] = s[5] + freepath**6;                              #.
    return s