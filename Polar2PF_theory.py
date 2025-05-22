# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:29:18 2020

@author: Maciej Lipok
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
pi=np.pi
cos=np.cos
sin=np.sin


#objective collecting parameters, calculated basing on the theoretical functions available in repository in file "Fitting model description" 
K1=2.945
K2=0.069
K3=1.016

#dichroic mirror parameters
delta=0.98845 #for excitation wavelength 810nm, for wavelength 750nm, delta=0.99662. Delta is an ellipticity of used dichroic mirror.
gamma=0.01 #dichroism of dichroic mirror

#Exemplary parameters depending on the orientation of molecule and fluorophore
psi=1*pi/180 # Angle between a fluorophore and matrix molecule (from 0 to 180 degrees)
dpsi=1*pi/180 #Deviation of PSI angle of fluorophore
phi=10*pi/180 #Orientation(angle) of a matrix (from 0 to 180)
#all the angles are calculated from the begining of the plus x axis



def f(ith,iphi):
    return np.exp(-((ith-psi)/dpsi)**2)/(2*pi*sin(psi)*dpsi*pi**0.5)
def uax(ith, iphi,phi):
    return -sin(phi)*sin(ith)*sin(iphi)+cos(phi)*cos(ith)
def uay(ith,iphi,phi):
    return cos(phi)*sin(ith)*sin(iphi)+sin(phi)*cos(ith)
def uaz(ith,iphi): 
    return -sin(ith)*cos(iphi)
def JX(ith,iphi,phi):
    return K1*uax(ith,iphi,phi)**2+K2*uay(ith,iphi,phi)**2+K3*uaz(ith,iphi)**2
def JY(ith,iphi,phi):
    return K2*uax(ith,iphi,phi)**2+K1*uay(ith,iphi,phi)**2+K3*uaz(ith,iphi)**2
def fxxxxx(phi):
    return integrate.dblquad(lambda ith, iphi: JX(ith,iphi,phi)*uax(ith,iphi,phi)**4*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fxyyyy(phi):
    return integrate.dblquad(lambda ith, iphi: JX(ith,iphi,phi)*uay(ith,iphi,phi)**4*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fxxyyy(phi):
    return integrate.dblquad(lambda ith, iphi: JX(ith,iphi,phi)*uax(ith,iphi,phi)*uay(ith,iphi,phi)**3*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fxxxyy(phi):
    return integrate.dblquad(lambda ith, iphi: JX(ith,iphi,phi)*uax(ith,iphi,phi)**2*uay(ith,iphi,phi)**2*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fxxxxy(phi):
    return integrate.dblquad(lambda ith, iphi: JX(ith,iphi,phi)*uay(ith,iphi,phi)*uax(ith,iphi,phi)**3*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fyyyyy(phi):
    return integrate.dblquad(lambda ith, iphi: JY(ith,iphi,phi)*uay(ith,iphi,phi)**4*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fyxxxx(phi):
    return integrate.dblquad(lambda ith, iphi: JY(ith,iphi,phi)*uax(ith,iphi,phi)**4*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fyxxyy(phi):
    return integrate.dblquad(lambda ith, iphi: JY(ith,iphi,phi)*uax(ith,iphi,phi)**2*uay(ith,iphi,phi)**2*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fyxxxy(phi):
    return integrate.dblquad(lambda ith, iphi: JY(ith,iphi,phi)*uay(ith,iphi,phi)*uax(ith,iphi,phi)**3*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]
def fyxyyy(phi):
    return integrate.dblquad(lambda ith, iphi: JY(ith,iphi,phi)*uax(ith,iphi,phi)*uay(ith,iphi,phi )**3*f(ith,iphi)*sin(ith),0,2*pi,0,pi)[0]

data_size=180 #size of data column
m=np.arange(1, data_size)
l=np.arange(1, data_size)
o=np.arange(1, data_size)
    
alph=2*o*pi/data_size
iphi=l*pi/data_size # angle of orintation of a dipole from z direction
ith=2*m*pi/data_size # angle of orientation of a dipole in xy plane

Ex4,Ey4,Ex2Ey2,Ex3Ey,ExEy3=np.arange(1.0, data_size),np.arange(1.0, data_size),np.arange(1.0,data_size),np.arange(1.0, data_size),np.arange(1.0, data_size)

for x in range(0, data_size-1):
        Ex4[x]=integrate.quad(lambda fi: ((1-gamma)*cos(alph[x])*cos(fi))**4,0,2*pi)[0]
        Ey4[x]=integrate.quad(lambda fi:(sin(alph[x])*cos(fi+delta))**4,0,2*pi)[0]
        Ex2Ey2[x] =integrate.quad(lambda fi:(((1-gamma)*cos(alph[x])*cos(fi))**2)*((sin(alph[x])*cos(fi+delta))**2),0,2*pi)[0]
        Ex3Ey[x] =integrate.quad(lambda fi:((((1-gamma)*cos(alph[x])*cos(fi))**3)*sin(alph[x])*cos(fi+delta)),0,2*pi)[0]
        ExEy3[x] =integrate.quad(lambda fi:((1-gamma)*cos(alph[x])*cos(fi)*(sin(alph[x])*cos(fi+delta))**3),0,2*pi)[0]

Px=fxxxxx(phi)*Ex4+fxyyyy(phi)*Ey4 +6*fxxxyy(phi)*Ex2Ey2+4*fxxxxy(phi)*Ex3Ey+4*fxxyyy(phi)*ExEy3
Py=fyyyyy(phi)*Ey4+fyxxxx(phi)*Ex4+6*fyxxyy(phi)*Ex2Ey2+4*fyxyyy(phi)*ExEy3+4*fyxxxy(phi)*Ex3Ey


plt.figure()
plt.polar(alph,Px,'-r',alph,Py,'-b')
plt.title('Px - red line Py - blue line')



