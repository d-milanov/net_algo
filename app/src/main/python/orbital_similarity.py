# Python implementations of orbital similarity criteria:
# - two variants of D-criterion by Southworth & Hawkins
# - D-criterion by Drummond
# - rho_2 distance by Kholshevnikov
# - DH-criterion by Jopek
#
# Please cite our article 'Relaxed triangle inequality for the orbital similarity criterion by 
# Southworth and Hawkins and its variants' (https://link.springer.com/article/10.1007/s10569-019-9884-6)
# if you use this code for your public work.
#
# Author: Yulia Milanova (yulia.milanova@gmail.com)
# Year:   2018

import math

PRECISION = 1e-5

# calculates D-criterion by Southworth & Hawkins
# Southworth R., Hawkins G.: Statistics of meteor streams. Smithsonian Contributions
# to Astrophysics 7, 261 (1963)
def DSH(q1, e1, i1, Omega1, omeg1, q2, e2, i2, Omega2, omeg2):
    a = AngleBetweenPlanes(i1,Omega1,i2,Omega2)
    b = AngleBetweenPP(i1,omeg1,Omega1,i2,omeg2,Omega2,a)
    return math.sqrt((e2-e1)**2 + (q2-q1)**2 + (2*math.sin(a/2))**2 + (((e2+e1)/2)**2)*((2*math.sin(b/2))**2))

# calculates the aproximation of D-criterion by Southworth & Hawkins in the case of small i1, i2
# Southworth R., Hawkins G.: Statistics of meteor streams. Smithsonian Contributions
# to Astrophysics 7, 261 (1963)
def DSH2(q1, e1, i1, Omega1, omeg1, q2, e2, i2, Omega2, omeg2):
    a = AngleBetweenPlanes(i1,Omega1,i2,Omega2)
    b = ((Omega1+omeg1)-(Omega2+omeg2))
    return math.sqrt((e2-e1)**2 + (q2-q1)**2 + (2*math.sin(a/2))**2 + (((e2+e1)/2)**2)*((2*math.sin(b/2))**2))

# calculates D-criterion of Drummond
# Drummond J.: A test of comet and meteor shower associations. Icarus 45(3), 545-553 (1981)
def DD(q1, e1, i1, Omega1, omeg1, q2, e2, i2, Omega2, omeg2):
    a = AngleBetweenPlanes(i1,Omega1,i2,Omega2)
    b = AngleBetweenPPO(i1,omeg1,Omega1,i2,omeg2,Omega2)
    return math.sqrt(((e2-e1)/(e2+e1))**2 + ((q2-q1)/(q2+q1))**2 + (a/(math.pi))**2 + (((e2+e1)/2)**2)*((b/math.pi)**2))

# calculates the rho_2 distance by Kholshevnikov
# Kholshevnikov K., Kokhirova G., Babadzhanov P., Khamroev U.: Metrics in the space of orbits 
# and their application to searching for celestial objects of
# common origin. Monthly Notices of the Royal Astronomical Society 462(2), 2275-2283 (2016), formula (15)
def RHO2(q1, e1, i1, Omega1, omeg1, q2, e2, i2, Omega2, omeg2):
    p1 = q1*(1+e1)
    p2 = q2*(1+e2)
    rI = cosI(i1,Omega1,i2,Omega2)
    rP = cosP(i1,Omega1,omeg1,i2,Omega2,omeg2)
    R = (1 + e1**2)*p1 + (1 + e2**2)*p2 - 2*(math.sqrt(p1*p2))*(rI + e1*e2*rP)
    if (R < 0)and(R>-PRECISION):
        R = 0
    return math.sqrt(R)

# calculates DH-criterion by Jopek
# Jopek T.: Remarks on the meteor orbital similarity D-criterion. Icarus 106(2), 603-607 (1993)
def DH(q1, e1, i1, Omega1, omeg1, q2, e2, i2, Omega2, omeg2):
    a = AngleBetweenPlanes(i1,Omega1,i2,Omega2)
    b = AngleBetweenPP(i1,omeg1,Omega1,i2,omeg2,Omega2,a)
    return math.sqrt((e2-e1)**2 + ((q2-q1)/(q2+q1))**2 + (2*math.sin(a/2))**2 + (((e2+e1)/2)**2)*((2*math.sin(b/2))**2))

# calculates the angle between the planes on which each of the two orbits lie (I_21)
def AngleBetweenPlanes(i1,Omega1,i2,Omega2):
    Z = math.cos(i1)*math.cos(i2)+math.sin(i1)*math.sin(i2)*math.cos(Omega2-Omega1)
    if (Z>1) and (Z < 1 + PRECISION):
        Z = 1
    return math.acos(Z)

# calculates the angle between the two orbits respective perihelion points (Pi_21)
def AngleBetweenPP(i1,omeg1,Omega1,i2,omeg2,Omega2,func1):
    if math.fabs(Omega2-Omega1) <= math.pi:
        P = omeg2 - omeg1 + 2*math.asin(math.cos((i2+i1)/2)*math.sin((Omega2-Omega1)/2)*(1/math.cos(func1/2)))
        return P
    else:
        P = omeg2 - omeg1 - 2*math.asin(math.cos((i2+i1)/2)*math.sin((Omega2-Omega1)/2)*(1/math.cos(func1/2)))
        return P
    

def cosI(i1,Omega1,i2,Omega2):
    return math.cos(i1)*math.cos(i2) + math.sin(i1)*math.sin(i2)*math.cos(Omega1-Omega2)

def cosP(i1,Omega1,omeg1,i2,Omega2,omeg2):
    return math.sin(i1)*math.sin(i2)*math.sin(omeg1)*math.sin(omeg2)+(math.cos(omeg1)*
                    math.cos(omeg2) + math.cos(i1)*math.cos(i2)*math.sin(omeg1)*
                    math.sin(omeg2))*math.cos(Omega1-Omega2) + (math.cos(i2)*math.cos(omeg1)*
                            math.sin(omeg2) - math.cos(i1)*math.cos(omeg2)*math.sin(omeg1))*math.sin(Omega1-Omega2)


# the angle between the perihelion points on each orbit    (Theta_21)    
def AngleBetweenPPO(i1,omeg1,Omega1,i2,omeg2,Omega2):
    l1 = Omega1 + math.atan(math.cos(i1)*math.tan(omeg1))
    l2 = Omega2 + math.atan(math.cos(i2)*math.tan(omeg2))
    if math.cos(omeg1) < 0 and math.cos(omeg2) < 0:
        l2 = l2 + math.pi
        l1 = l1 + math.pi
    elif math.cos(omeg1) < 0 and math.cos(omeg2) >= 0:
        l1 = l1 + math.pi
    elif math.cos(omeg1) >= 0 and math.cos(omeg2) < 0:
        l2 = l2 + math.pi
    T = math.acos(math.sin(math.asin(math.sin(i1)*math.sin(omeg1)))*
                  math.sin(math.asin(math.sin(i2)*math.sin(omeg2)))+
                  math.cos(math.asin(math.sin(i1)*math.sin(omeg1)))*
                  math.cos(math.asin(math.sin(i2)*math.sin(omeg2)))*math.cos(l2 - l1))
    return T
 
