import numpy as np

PhiA, LamA = (50.25), (20.75)
PhiB, LamB = (50.00), (20.75)
PhiC, LamC = (50.25), (21.25)
PhiD, LamD = (50.00), (21.25)
PhiS, LamS = (50.125), (21.0)
PhiM, LamM = (50.12527054195198), (21.00065090208314)

#dane dla grs80
ag = 6378137
e2g = 0.00669437999013
#dane dla elipsoidy krasowskiego
ak = 6378245
e2k = 0.0066934215520398155
kx0 = -33.4297
ky0 = 146.5746
kz0 = 76.2865
kk = 0.8407728 * 0.000001
kex = (-0.35867 / 3600)
key = (-0.05283 / 3600)
kez = (0.84354 / 3600)


def geo_to_xyz(phi, lamb, H, a, e2):
    phi = np.deg2rad(phi)
    lamb = np.deg2rad(lamb)
    N = a / np.sqrt(1 - e2 * (np.sin(phi) ** 2))
    x = (N + H) * np.cos(phi) * np.cos(lamb)
    y = (N + H) * np.cos(phi) * np.sin(lamb)
    z = (N * (1 - e2) + H) * np.sin(phi)
    return x, y, z


def hirvonen(x, y, z, a, e2):
    r = np.sqrt((x ** 2) + (y ** 2))
    fi1 = np.arctan((z/r) * ((1-e2) ** -1))
    n = a / np.sqrt(1 - e2 * np.sin(fi1) ** 2)
    h = (r / np.cos(fi1)) - n
    fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    epsilon = np.deg2rad(0.00005/3600)

    while abs(fi2-fi1) > epsilon:
        fi1 = fi2
        n = a / np.sqrt(1 - e2 * np.sin(fi1) ** 2)
        h = (r / np.cos(fi1)) - n
        fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)

    n = a / np.sqrt(1 - e2 * np.sin(fi2) ** 2)
    h = (r / np.cos(fi2)) - n
    lam = np.arctan(y/x)

    fi_d = np.rad2deg(fi2)
    lam_d = np.rad2deg(lam)
    return fi_d, lam_d, h


def transformacja(xp, yp, zp, x0, y0, z0, ex, ey, ez, k):
    ex = np.deg2rad(ex)
    ey = np.deg2rad(ey)
    ez = np.deg2rad(ez)
    m_wektorow = np.array([[x0], [y0], [z0]])
    m_wspolrzednych = np.array([[xp], [yp], [zp]])
    m_obrotu = np.array([[k, ez, -ey], [-ez, k, ex], [ey, -ex, k]])
    m_koncowa = m_wspolrzednych + m_obrotu.dot(m_wspolrzednych) + m_wektorow
    return np.transpose(m_koncowa)

def dectodeg(x):
    degrees = int(x)
    minutes = int((x - degrees) * 60)
    seconds = (x - degrees - minutes / 60) * 3600
    seconds = round(seconds, 5)
    return degrees,minutes,seconds

ax, ay, az = geo_to_xyz(PhiA, LamA, 0, ag, e2g)
bx, by, bz = geo_to_xyz(PhiB, LamB, 0, ag, e2g)
cx, cy, cz = geo_to_xyz(PhiC, LamC, 0, ag, e2g)
dx, dy, dz = geo_to_xyz(PhiD, LamD, 0, ag, e2g)
sx, sy, sz = geo_to_xyz(PhiS, LamS, 0, ag, e2g)
mx, my, mz = geo_to_xyz(PhiM, LamM, 0, ag, e2g)
print('GRS80')
print('xyz A na GRS80:', ax, ay, az)
print('xyz B na GRS80:', bx, by, bz)
print('xyz C na GRS80:', cx, cy, cz)
print('xyz D na GRS80:', dx, dy, dz)
print('xyz S na GRS80:', sx, sy, sz)
print('xyz M na GRS80:', mx, my, mz)
print('Krasowski')
print('xyz A na e Krasowskiego', transformacja(ax, ay, az, kx0, ky0, kz0, kex, key, kez, kk))
print('xyz B na e Krasowskiego', transformacja(bx, by, bz, kx0, ky0, kz0, kex, key, kez, kk))
print('xyz C na e Krasowskiego', transformacja(cx, cy, cz, kx0, ky0, kz0, kex, key, kez, kk))
print('xyz D na e Krasowskiego', transformacja(dx, dy, dz, kx0, ky0, kz0, kex, key, kez, kk))
print('xyz S na e Krasowskiego', transformacja(sx, sy, sz, kx0, ky0, kz0, kex, key, kez, kk))
print('xyz M na e Krasowskiego', transformacja(mx, my, mz, kx0, ky0, kz0, kex, key, kez, kk))
pa, la, ha = hirvonen(3821428.58996192, 1447942.18763557, 4880698.98862972, ak, e2k)
pb, lb, hb = hirvonen(3841385.34575167, 1455503.06544474, 4862870.95955055, ak, e2k)
pc, lc, hc = hirvonen(3808648.7665734, 1481235.17275336, 4880699.04979542, ak, e2k)
pd, ld, hd = hirvonen(3828538.78247279, 1488969.91604633, 4862871.02103567, ak, e2k)
ps, ls, hs = hirvonen(3825045.96875475, 1468430.08850005, 4871796.54802818, ak, e2k)
pm, lm, hm = hirvonen(3825007.72752856, 1468465.26620259, 4871815.84098864, ak, e2k)
print('.........')
print('A', dectodeg(pa), dectodeg(la), ha)
print('B', dectodeg(pb), dectodeg(lb), hb)
print('C', dectodeg(pc), dectodeg(lc), hc)
print('D', dectodeg(pd), dectodeg(ld), hd)
print('S', dectodeg(ps), dectodeg(ls), hs)
print('M', dectodeg(pm), dectodeg(lm), hm)




