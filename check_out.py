# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:43:27 2017

@author: TH Lab
"""

def check_out(z1,thickness,sigmatot):
    '''
    This function check if the particle is escaped (1) outside or not (0).
    
    z1 : position to check 
    
    thickness : distance on the z-axis of the 1D media - integer
    
    sigmatot : total cross section of the media, if it is equal to 0, it is like a case without media - float
    '''
    if(z1 >= thickness or sigmatot==0 or z1 < 0):                         
        return 1;
    return 0