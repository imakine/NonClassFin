# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:31:03 2017

@author: TH Lab
"""
import random

def posInit(thickness,sourceregion,region):
    '''
    This function initialize the first position of the particle according to the source and its kind
    
    thickness : distance on the z-axis of the 1D media - Integer
    
    sourceregion: the distance where we have a constant source - Float
    
    region: True if regional source False if uniform source - Boolean
    '''
    randnum = random.random();
    if region:
        z0=thickness/2+(sourceregion*(2*randnum-1))/2;  #Initial rqndom position of particles according to the regional source : Source = Q between -0.5 + thickness/2 and +0.5 + thickness/2
    else:
        z0=thickness*randnum;                           #Uniform source, the particle can be born everywhere in the media
    return z0

def muInit():
    '''
    This function gives the first direction of travel for the particle, this element is given by mu. The value of mu is between -1 and 1.
    '''
    randnum = random.random()
    w0=2*randnum-1   #Initial direction of the particle randnum --> between 0 & 1, modify the expression to have -1 to 1
    return w0