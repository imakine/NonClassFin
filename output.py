# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:34:19 2017

@author: TH Lab
"""
import matplotlib.pyplot as plt

def output(sp,c,s1,s2,s3,s4,s5,s6,s,finaltime,initialtime,flux_local,thickness,sigmaabs,sigmascat,sigmatot,Q,n,step,numscattot,numesc,numdeath,ERRORpic,std,variancepic):
    '''
    This file generate the output files
    '''
    
    filename = 'realend';
    file = open(filename + ".txt", "w")
    fluxPSP1 = open("flux_spCRsec" + str(sp) + str(c) + ".txt","w")
    coordPSP1 = open("coord_spCRsec" + str(sp) +str(c) + ".txt","w")
    variancePSP1 = open("variance_spCRsec" + str(sp) +str(c) + ".txt","w");
    deviationPSP1 = open("deviation_spCRsec" + str(sp) +str(c) + ".txt","w");
    ERRORPSP1 = open("ERRORPSP1_spCRsec" + str(sp) +str(c) + ".txt","w");
    print("Time elapsed during the running of the code : ",finaltime - initialtime, "seconds");
    print(" ");
    print(" ");
    print('max of the flux = ', max(flux_local))
    print("average on bins max = ", (flux_local[flux_local.index(max(flux_local))]+flux_local[flux_local.index(max(flux_local))+1])/2)
    print(" ");
    print("Maximum for the standart deviation = "+ str(max(std)) + " at the position " + str(std.index(max(std))) )
    print("s = ",s)
    plt.plot(flux_local)
    plt.ylabel('Flux(z)')
    plt.show()

    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ------------------------------------------------------------\n")
    file.write("   |                                                          |\n")
    file.write("   |                 MONTE-CARLO CODE FOR                     |\n")
    file.write("   |           SP1 EQUATION (1D-Ponctual source)              |\n")
    file.write("   |             Done at Berkeley by Ilker Makine             |\n")
    file.write("   |                      OUTPUT FILE                         |\n")
    file.write("   |                                                          |\n")
    file.write("   |                                                          |\n")
    file.write("   ------------------------------------------------------------\n")
    file.write("   ---------------------------INPUT PARAMETERS---------------------\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   " + " Thickness of the media = " + str(thickness) + "\n")
    file.write("   " + " Absorption cross-section = " + str(sigmaabs)+ "\n")
    file.write("   " + " Scattering cross-section = " + str(sigmascat)+ "\n")
    file.write("   " + " Total cross-section = " + str(sigmatot)+ "\n")
    file.write("   " + " Source of neutrons = " + str(Q)+ "\n")
    file.write("   " + " Number of particles = " + str(n)+ "\n")
    file.write("   " + " Precision in the media (step) = " + str(step)+ "\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------TIME ELAPSED---------------------\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   " + str(finaltime - initialtime) + " seconds are elapsed during the running of the code")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------MOMENTS--------------------------\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^m         | Theoretical values | Computational values |        Errors        |\n")
    file.write("   |                  |    (Classical)     |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^1         | "+str(s1)+" "*(22-len(str(s1))) +         " |"      +    str(s[0])+" "*(22-len(str(s[0])))+"|"          +str((s[0]- s1)/s1) + " "*(22-len(str((s[0]- s1)/s1)))+     "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^2         | "+str(s2)+" "*(22-len(str(s2))) +         " |"      +    str(s[1])+" "*(22-len(str(s[1])))+"|"        +str((s[1]- s2)/s2) +" "*(22-len(str((s[1]- s2)/s2))) +    "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^3         | "+str(s3)+" "*(22-len(str(s3))) +         " |"      +    str(s[2])+" "*(22-len(str(s[2])))+"|"           +str((s[2]- s3)/s3) +  " "*(22-len(str((s[2]- s3)/s3))) +      "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^4         | "+str(s4)+" "*(22-len(str(s4))) +         " |"      +    str(s[3])+" "*(22-len(str(s[3])))+"|"         +str((s[3]- s4)/s4)+    " "*(22-len(str((s[3]- s4)/s4)))  +     "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^5         | "+str(s5)+" "*(22-len(str(s5))) +         " |"      +    str(s[4])+" "*(22-len(str(s[4])))+"|"         +str((s[4]- s5)/s5)+  " "*(22-len(str((s[4]- s5)/s5)))   +      "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^6         | "+str(s6)+" "*(22-len(str(s6))) +         " |"      +    str(s[5])+" "*(22-len(str(s[5])))+"|"         +str((s[5]- s6)/s6)+ " "*(22-len(str((s[5]- s6)/s6))) +        "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------EVENTS COUNTER---------------------\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   " + str(numscattot) + " scatterings occurs during the simulation\n")
    file.write("   " + str(numesc) + " escapes occurs during the simulation\n")
    file.write("   " + str(numdeath) + " absorptions occurs during the simulation\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------FLUX(Z)----------------------------\n")
    file.write(" --- Z ---         ---FLUX(Z)--- \n")
    z=0;
    file.close();
    for i in range(len(flux_local)):
        fluxPSP1.write("   "+ str(flux_local[i])+"\n");
        coordPSP1.write("   "+ str(z)+"\n");
        z = z + step
        variancePSP1.write("   "+ str(variancepic[i])+"\n"); #THIS IS NOT VARIANCE; IT IS THE SECOND MOMENT OF FLUX
        deviationPSP1.write("   "+ str(std[i])+"\n");
        ERRORPSP1.write("   "+ str(ERRORpic[i])+"\n")
    deviationPSP1.close()
    variancePSP1.close()
    fluxPSP1.close();
    coordPSP1.close();
    ERRORPSP1.close();