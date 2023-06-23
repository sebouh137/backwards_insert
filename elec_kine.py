import numpy as np

#Return theta in rad
def etatotheta(eta):
    theta = 2*np.arctan(np.exp(-1.*eta))
    return theta

#Return theta in deg
def etatothetadeg(eta):
    theta = np.rad2deg( etatotheta(eta) )
    return theta 

def Q2toEp(Q2,E,theta):
    cos2theta2 = np.cos(theta/2.)*np.cos(theta/2.)
    Ep = Q2/(4.*E*cos2theta2)
    return Ep

def ytoEp(y,E,theta):
    costheta = np.cos(theta)
    Ep = (1.-y)*(2.*E)/(1.-costheta)
    return Ep

def EptoQ2(E,Ep,theta):
    cos2theta2 = np.cos(theta/2.)*np.cos(theta/2.)
    Q2 = 4.*E*Ep*cos2theta2
    return Q2

def ytoQ2(y,E,theta):
    cos2theta2 = np.cos(theta/2.)*np.cos(theta/2.)
    Ep = ytoEp(y,E,theta)
    Q2 = 4.*E*Ep*cos2theta2
    return Q2

def Eptoy(E,Ep,theta):
    costheta = np.cos(theta)
    y = 1. - Ep*(1.-costheta)/(2.*E)
    return y

def Q2toy(Q2,E,theta):
    costheta = np.cos(theta)
    Ep = Q2toEp(Q2,E,theta)
    y = 1. - (Ep)/(2.*E)*(1.-costheta)
    return y


    
