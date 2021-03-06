import numpy as np
import math

phiA = 50.25
phiB = 50.00
phiC = 50.25
phiD = 50.00
lamA = 20.75
lamB = 20.75
lamC = 21.25
lamD = 21.25

#obliczenie punktu sredniej szerokosci
def MidLatitudePoint(pA, pB, pC, pD, lA, lB, lC, lD):
    midPhi = (pA + pB + pC + pD) / 4
    midLam = (lA + lB + lC + lD) / 4
    return midPhi, midLam

MidLatitudePoint(phiA, phiB, phiC, phiD, lamA, lamB, lamC, lamD)

def Vincent(lam2, phi2, lam1, phi1):
    a = 6378137
    e2 = 0.00669437999013
    #1
    b = a * np.sqrt(1 - e2)
    f = 1 - (b / a)
    #2
    L = np.deg2rad(lam2 - lam1) #lambD - lambA
    #3
    Ua = np.arctan((1 - f) * np.tan(np.deg2rad(phi1)))
    Ub = np.arctan((1 - f) * np.tan(np.deg2rad(phi2)))
    #4
    L1 = L

    while True:
        #5
        sinSigma = np.sqrt((((np.cos(Ub)) * np.sin(L1)) ** 2) + (((np.cos(Ua) * np.sin(Ub))
                                                                 - (np.sin(Ua) * np.cos(Ub) * np.cos(L1))) ** 2))
        #6
        cosSigma = (np.sin(Ua) * np.sin(Ub)) + (np.cos(Ua) * np.cos(Ub) * np.cos(L1))
        #7
        Sigma = np.arctan(sinSigma / cosSigma)
        #8
        sinAlfa = (np.cos(Ua) * np.cos(Ub) * np.sin(L1)) / sinSigma
        #9
        cos2Alfa = 1 - (sinAlfa ** 2)
        #10
        cos2sig = cosSigma - ((2 * np.sin(Ua) * np.sin(Ub)) / cos2Alfa)
        #11
        C = (f / 16) * cos2Alfa * (4 + f * (4 - 3 * cos2Alfa))
        #12
        L2 = L + (1 - C) * f * sinAlfa * (Sigma + C * np.sin(Sigma) * (cos2sig + C * cosSigma * (-1 + 2 * (cos2sig ** 2))))
        if abs(L2 - L1) < np.deg2rad(0.000001 / 3600):
            break
        else:
            L1 = L2

    #13
    u2 = (((a ** 2) - (b ** 2)) / (b ** 2)) * cos2Alfa
    #14
    A1 = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    #15
    B1 = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    #16
    deltaSigma = B1 * sinSigma * (cos2sig + (B1 * (cosSigma * (-1 + 2 * (cos2sig ** 2)))) / 4 -
                                  (B1 * cos2sig * (-3 + 4 * (sinSigma ** 2)) * (-3 + 4 * (cos2sig ** 2))) / 6)

    #17
    Sab = (b * A1 * (Sigma - deltaSigma))

    num1 = np.cos(Ub) * np.sin(L1)
    denum1 = (np.cos(Ua) * np.sin(Ub)) - (np.sin(Ua) * np.cos(Ub) * np.cos(L1))
    if (num1 > 0 and denum1 > 0):
        Aab = np.rad2deg(np.arctan(num1 / denum1))
    elif (num1 > 0 and denum1 < 0):
        Aab = np.rad2deg(np.arctan(num1 / denum1) + np.pi)
    elif (num1 < 0 and denum1 < 0):
        Aab = np.rad2deg(np.arctan(num1 / denum1) + np.pi)
    elif (num1 < 0 and denum1 > 0):
        Aab = np.rad2deg(np.arctan(num1 / denum1) + 2 * np.pi)

    num2 = np.cos(Ua) * np.sin(L1)
    denum2 = (-1 * np.sin(Ua) * np.cos(Ub)) + (np.cos(Ua) * np.sin(Ub) * np.cos(L1))
    if (num2 > 0 and denum2 > 0):
        Aba = np.rad2deg(np.arctan(num2 / denum2) + np.pi)
    elif (num2 > 0 and denum2 < 0):
        Aba = np.rad2deg(np.arctan(num2 / denum2) + 2 * np.pi)
    elif (num2 < 0 and denum2 < 0):
        Aba = np.rad2deg(np.arctan(num2 / denum2) + 2 * np.pi)
    elif (num2 < 0 and denum2 > 0):
        Aba = np.rad2deg(np.arctan(num2 / denum2) + 3 * np.pi)
    return Sab, Aab, Aba

def Kivioj():
    a = 6378137
    e2 = 0.00669437999013
    Sab, Aab, Aba = Vincent(lamD, phiD, lamA, phiA)
    Sab = Sab / 2
    ds = 1111
    n = Sab / ds #ds = 1111
    B = np.deg2rad(phiA)
    L = np.deg2rad(lamA)
    Aab = np.deg2rad(Aab)

    for i in range (int(n)):
        M = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(B) ** 2)) ** 3))
        N = a / (np.sqrt(1 - e2 * (np.sin(B) ** 2)))
        dB = (ds * np.cos(Aab)) / M
        dA = (np.sin(Aab) * np.tan(B) * ds) / N
        phiM = B + dB / 2
        Am = Aab + dA / 2
        Mm = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(phiM) ** 2)) ** 3))
        Nm = a / (np.sqrt(1 - e2 * (np.sin(phiM) ** 2)))
        dphiM = (ds * np.cos(Am)) / Mm
        dlamM = (ds * np.sin(Am)) / (Nm * np.cos(phiM))
        dAm = (np.sin(Am) * np.tan(phiM) * ds) / Nm
        B = B + dphiM
        L = L + dlamM
        Aab = Aab + dAm


    ds1 = Sab % 1111
    M = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(B) ** 2)) ** 3))
    N = a / (np.sqrt(1 - e2 * (np.sin(B) ** 2)))
    dB = (ds1 * np.cos(Aab)) / M
    dA = (np.sin(Aab) * np.tan(B) * ds1) / N
    phiM = B + dB / 2
    Am = Aab + dA / 2
    Mm = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(phiM) ** 2)) ** 3))
    Nm = a / (np.sqrt(1 - e2 * (np.sin(phiM) ** 2)))
    dphiM = (ds1 * np.cos(Am)) / Mm
    dlamM = (ds1 * np.sin(Am)) / (Nm * np.cos(phiM))
    dAm = (np.sin(Am) * np.tan(phiM) * ds1) / Nm
    B1 = B + dphiM
    L1 = L + dlamM
    Aab1 = Aab + dAm
    print()

    return np.rad2deg(B1), np.rad2deg(L1), np.rad2deg(Aab1)

midPhi, midLam = MidLatitudePoint(phiA, phiB, phiC, phiD, lamA, lamB, lamC, lamD)
B, L, Aab = Kivioj()

def area(p1,p2,l1,l2):
    a = 6378137
    e2 = 0.00669437999013
    e = np.sqrt(e2)
    b = a * np.sqrt(1 - e2)
    p1 = np.deg2rad(p1)
    p2 = np.deg2rad(p2)
    l1 = np.deg2rad(l1)
    l2 = np.deg2rad(l2)
    phi1 = (np.sin(p1) / (1 - e2 * (np.sin(p1) ** 2))) + (1/(2 * e) * np.log((1 + e * np.sin(p1)) / (1 - e * np.sin(p1))))
    phi2 = (np.sin(p2) / (1 - e2 * (np.sin(p2) ** 2))) + (1/(2 * e) * np.log((1 + e * np.sin(p2)) / (1 - e * np.sin(p2))))
    P = (((b ** 2) * (l1 - l2)) / 2) * (phi2 - phi1)
    print("Pole czworok??ta:", P)

def dectodeg(x):
    degrees = int(x)
    minutes = int((x - degrees) * 60)
    seconds = (x - degrees - minutes / 60) * 3600
    seconds = round(seconds, 5)
    return degrees,minutes,seconds

print("Punkt ??redniej szeroko??ci:", MidLatitudePoint(phiA, phiB, phiC, phiD, lamA, lamB, lamC, lamD))
print("Azymuty AD:", Vincent(lamD, phiD,lamA,phiA))
print("Wsp????rzedne punktu ??rodkowego:", B, L)
print("Odleg??o???? mi??dzy punktem ??r. i ??r. szeroko??ci i azymuty:", Vincent(L, B, midLam, midPhi))
area(phiD,phiA,lamD,lamA)


print("Zamiana na stopnie:", dectodeg(127.68147015704359))
print(dectodeg(308.0651945576038))
print(dectodeg(50.12527054195198))
print(dectodeg(21.00065090208314))
print(dectodeg(57.11621259889227))
print(dectodeg(237.11671213128352))



