# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:20:15 2017

@author: TH Lab
"""

def theoretical(sp,sigmatot):
    '''
    This function gives the theoretical moments for the SPn equations - These results come from computation
    
    sp : Choice of the user for SP1 or SP2 or SP3 - Integer 1 or 2 or 3
    
    sigmatot : total cross section of the media, if it is equal to 0, it is like a case without media - float
    '''
    if sp==1:
        s1 = 1.1547/sigmatot;
        s2 = 2/sigmatot**2;
        s3 = 4.6188/sigmatot**3;
        s4 = 13.3333/sigmatot**4;
        s5 = 48.188/sigmatot**5;
        s6 = 186.66667/sigmatot**6;
    elif sp==2:
        s1 = 0.8606/sigmatot;
        s2 = 2/sigmatot**2;
        s3 = 6.1968/sigmatot**3;
        s4 = 24/sigmatot**4;
        s5 = 111.5419/sigmatot**5;
        s6 = 604.8/sigmatot**6;
    elif sp==3:
        s1 = 1.042533/sigmatot;
        s2 = 2/sigmatot**2;
        s3 = 5.94625/sigmatot**3;
        s4 = 24/sigmatot**4;
        s5 = 120.734028/sigmatot**5;
        s6 = 720/sigmatot**6;
    return s1,s2,s3,s4,s5,s6    