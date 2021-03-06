import numpy as np
from shapely.geometry import Polygon

phiA = 50.25
phiB = 50.00
phiC = 50.25
phiD = 50.00
phiS = 50.125
phiM = 50.12527054195198
lamA = 20.75
lamB = 20.75
lamC = 21.25
lamD = 21.25
lamS = 21.0
lamM = 21.00065090208314

def GK(phi, lamb, lamb0):
    a = 6378137
    e2 = 0.00669437999013
    a2 = a ** 2
    b2 = a2 * (1 - e2)
    ep2 = (a2 - b2) / b2
    lamb0 = np.deg2rad(lamb0)
    lamb = np.deg2rad(lamb)
    phi = np.deg2rad(phi)
    deltalamb = lamb - lamb0
    t = np.tan(phi)
    eta2 = ep2 * (np.cos(phi) ** 2)
    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    a0 = 1 - (e2/4) - ((3 * (e2 ** 2)) / 64) - ((5 * (e2 ** 3)) / 256)
    a02 = 3/8 * (e2 + ((e2 ** 2) / 4) + ((15 * (e2 ** 3)) / 128))
    a4 = 15/256 * ((e2 ** 2) + ((3 * (e2 ** 3)) / 4))
    a6 = (35 * (e2 ** 3)) / 3072
    sigma = a * (a0 * phi - a02 * np.sin(2 * phi) + a4 * np.sin(4 * phi) - a6 * np.sin(6 * phi))
    xkg = sigma + ((deltalamb ** 2) / 2) * N * np.sin(phi) * np.cos(phi) * (1 + ((deltalamb ** 2) / 12) * \
          (np.cos(phi) ** 2) * (5 - t ** 2 + 9 * eta2 + (4 * (eta2 ** 2))) + ((deltalamb ** 4) / 360) * \
          (np.cos(phi) ** 4) * (61 - (58 * (t ** 2)) + (t ** 4) + (270 * eta2) - (330 * eta2 * (t ** 2))))
    ykg = deltalamb * N * np.cos(phi) * (1 + ((deltalamb ** 2) / 6) * (np.cos(phi) ** 2) * (1 - t ** 2 + eta2) + \
                                         ((deltalamb ** 4) / 120) * (np.cos(phi) ** 4) * (5 - 18 * (t **2) + (t ** 4) + \
                                                                                          14 * eta2 - 58 * eta2 * (t ** 2)))
    return xkg, ykg

def fromGK(xgk, ygk, lamb0):
    a = 6378137
    e2 = 0.00669437999013
    a2 = a ** 2
    b2 = a2 * (1 - e2)
    ep2 = (a2 - b2) / b2
    lamb0 = np.deg2rad(lamb0)
    a0 = 1 - (e2 / 4) - ((3 * (e2 ** 2)) / 64) - ((5 * (e2 ** 3)) / 256)
    a02 = 3 / 8 * (e2 + ((e2 ** 2) / 4) + ((15 * (e2 ** 3)) / 128))
    a4 = 15 / 256 * ((e2 ** 2) + ((3 * (e2 ** 3)) / 4))
    a6 = (35 * (e2 ** 3)) / 3072
    phi1 = xgk / (a * a0)
    sigma = a * (a0 * phi1 - a02 * np.sin(2 * phi1) + a4 * np.sin(4 * phi1) - a6 * np.sin(6 * phi1))

    while True:
        phi2 = phi1 + (xgk - sigma) / (a * a0)
        sigma = a * (a0 * phi2 - a02 * np.sin(2 * phi2) + a4 * np.sin(4 * phi2) - a6 * np.sin(6 * phi2))
        M = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(phi2) ** 2)) ** 3))
        N = a / (np.sqrt(1 - e2 * (np.sin(phi2) ** 2)))
        t = np.tan(phi2)
        eta2 = ep2 * (np.cos(phi2) ** 2)
        if abs(phi2 - phi1) < np.deg2rad(0.000001 / 3600):
            break
        else:
            phi1 = phi2

    phi1 = phi2 - ((ygk ** 2) * t / (2 * M * N)) * (1 - ygk ** 2 / (12 * (N ** 2)) * (5 + 3 * (t ** 2) + eta2 - 9 * eta2 *
                                                                                  (t ** 2) - 4 * (eta2 ** 2)) +
                                                  (ygk ** 4) / (360 * (N ** 4)) * (61 + 90 * (t ** 2) + 45 * (t ** 4)))
    lamb1 = lamb0 + (ygk / (N * np.cos(phi2))) * (1 - ygk ** 2 / (6 * (N ** 2)) * (1 + 2 * (t ** 2) + eta2) +
                                                (ygk ** 4) / (120 * (N ** 4)) * (5 + 28 * (t **2) + 24 * (t ** 4) + \
                                                                                          6 * eta2 + 8 * eta2 * (t ** 2)))
    phi1 = np.rad2deg(phi1)
    lamb1 = np.rad2deg(lamb1)
    return phi1, lamb1

xgka, ygka = GK(phiA, lamA, 19)
xgkb, ygkb = GK(phiB, lamB, 19)
xgkc, ygkc = GK(phiC, lamC, 19)
xgkd, ygkd = GK(phiD, lamD, 19)
xgks, ygks = GK(phiS, lamS, 19)
xgkm, ygkm = GK(phiM, lamM, 19)

def zone_number(lamb):
    zn = 0
    if 13.5 <= lamb < 16.5:
        zn = 5
    if 16.5 <= lamb < 19.5:
        zn = 6
    if 19.5 <= lamb < 22.5:
        zn = 7
    if 22.5 <= lamb < 25.5:
        zn = 8
    return zn

def u1992(xgk,ygk):
    m0 = 0.9993
    x1992 = m0 * xgk - 5300000
    y1992 = m0 * ygk + 500000
    return x1992, y1992

def u2000(phi,lamb):
    m0 = 0.999923
    xgk, ygk = GK(phi, lamb, 21)
    zn = zone_number(21)
    x2000 = m0 * xgk
    y2000 = m0 * ygk + zn * 1000000 + 500000
    return x2000, y2000

def u1992toGK(xgk, ygk):
    x,y = u1992(xgk,ygk)
    m0 = 0.9993
    xgk = (x + 5300000) / m0
    ygk = (y - 500000) / m0
    return xgk, ygk

def u2000toGK(phi, lamb):
    x,y = u2000(phi,lamb)
    m0 = 0.999923
    zn = zone_number(21)
    xgk = x / m0
    ygk = (y - (zn * 1000000) - 500000) / m0
    return xgk, ygk

def scale1992(xgk, ygk):
    m0 = 0.9993
    x, y = u1992toGK(xgk, ygk)
    phi, lamb = fromGK(x, y, 19)
    a = 6378137
    e2 = 0.00669437999013
    M = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(phi) ** 2)) ** 3))
    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    R = np.sqrt(M * N)
    m = 1 + (y ** 2) / (2 * R ** 2) + (y ** 2)/(24 * R ** 4)
    m1992 = m0 * m
    kappa = (m1992 - 1) * 1000
    return m1992, kappa

def scale2000(phi, lamb):
    m0 = 0.999923
    a = 6378137
    e2 = 0.00669437999013
    x, y = u2000toGK(phi,lamb)
    phi, lamb = fromGK(x, y, 21)
    M = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(phi) ** 2)) ** 3))
    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    R = np.sqrt(M * N)
    m = 1 + (y ** 2) / (2 * R ** 2) + (y ** 2) / (24 * R ** 4)
    m2000 = m0 * m
    kappa = (m2000 - 1) * 1000
    return m2000, kappa

def GKscale(xgk, ygk):
    phi, lamb = fromGK(xgk, ygk, 19)
    a = 6378137
    e2 = 0.00669437999013
    M = (a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(phi) ** 2)) ** 3))
    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    R = np.sqrt(M * N)
    m = 1 + (ygk ** 2) / (2 * R ** 2) + (ygk ** 2) / (24 * R ** 4)
    kappa = (m - 1) * 1000
    return m, kappa

def doublescaleha(m):
    m = m ** 2
    kappa = (m - 1) * 10000
    return m, kappa

xgk = [5570120.59683, 5542315.02551, 5543273.89200, 5571077.96021, 5570120.59683] 
ygk = [124812.22774, 125464.20085, 161308.28341, 160469.90666, 124812.22774]
x1992 = [266221.51256, 238435.40499, 239393.60028, 267178.20564, 266221.51256]
y1992 = [624724.85918, 625376.37591, 661195.36761, 660357.57773, 624724.85918]
x2000 = [5568256.03001, 5540450.35001, 5540450.35001, 5568256.03001, 5568256.03001]
y2000 = [7482170.56246, 7482077.45155,  7517922.54845,7517829.43754, 7482170.56246]

pgk = Polygon(zip(xgk, ygk))
p1992 = Polygon(zip(x1992, y1992))
p2000 = Polygon(zip(x2000, y2000))

print("A:", "GK:", GK(phiA, lamA, 19), "1992:", u1992(xgka, ygka), "2000:", u2000(phiA, lamA))
print("B:", "GK:", GK(phiB, lamB, 19), "1992:", u1992(xgkb, ygkb), "2000:", u2000(phiB, lamB))
print("C:", "GK:", GK(phiC, lamC, 19), "1992:", u1992(xgkc, ygkc), "2000:", u2000(phiC, lamC))
print("D:", "GK:", GK(phiD, lamD, 19), "1992:", u1992(xgkd, ygkd), "2000:", u2000(phiD, lamD))
print("S:", "GK:", GK(phiS, lamS, 19), "1992:", u1992(xgks, ygks), "2000:", u2000(phiS, lamS))
print("M:", "GK:", GK(phiM, lamM, 19), "1992:", u1992(xgkm, ygkm), "2000:", u2000(phiM, lamM))
print("A:", "mGK & kGK:", GKscale(xgka,ygka), "m1992 & k1992:", scale1992(xgka,ygka), "m2000 & k2000:", scale2000(phiA, lamA))
print("B:", "mGK & kGK:", GKscale(xgkb,ygkb), "m1992 & k1992:", scale1992(xgkb,ygkb), "m2000 & k2000:", scale2000(phiB, lamB))
print("C:", "mGK & kGK:", GKscale(xgkc,ygkc), "m1992 & k1992:", scale1992(xgkc,ygkc), "m2000 & k2000:", scale2000(phiC, lamC))
print("D:", "mGK & kGK:", GKscale(xgkd,ygkd), "m1992 & k1992:", scale1992(xgkd,ygkd), "m2000 & k2000:", scale2000(phiD, lamD))
print("S:", "mGK & kGK:", GKscale(xgks,ygks), "m1992 & k1992:", scale1992(xgks,ygks), "m2000 & k2000:", scale2000(phiS, lamS))
print("M:", "mGK & kGK:", GKscale(xgkm,ygkm), "m1992 & k1992:", scale1992(xgkm,ygkm), "m2000 & k2000:", scale2000(phiM, lamM))
print("A:", "mGK & kGK:", doublescaleha(1.0001927579507983), "m1992 & k1992:", doublescaleha(0.9994926230202327), "m2000 & k2000:", doublescaleha(0.9999269337500253))
print("B:", "mGK & kGK:", doublescaleha(1.0001945981454208), "m1992 & k1992:", doublescaleha(0.999494461926719), "m2000 & k2000:", doublescaleha(0.999926971293782))
print("C:", "mGK & kGK:", doublescaleha(1.000318628948981), "m1992 & k1992:", doublescaleha(0.9996184059087166), "m2000 & k2000:", doublescaleha(0.9999269337500253))
print("D:", "mGK & kGK:", doublescaleha(1.0003216713537757), "m1992 & k1992:", doublescaleha(0.999621446183828), "m2000 & k2000:", doublescaleha(0.999926971293782))
print("S:", "mGK & kGK:", doublescaleha(1.0002530136645005), "m1992 & k1992:", doublescaleha(0.9995528365549353), "m2000 & k2000:", doublescaleha(0.999923))
print("M:", "mGK & kGK:", doublescaleha(1.000253175765519), "m1992 & k1992:", doublescaleha(0.9995529985424831), "m2000 & k2000:", doublescaleha(0.9999230000267982))
print("GK:", pgk.area, "1992:", p1992.area, "2000:", p2000.area)
