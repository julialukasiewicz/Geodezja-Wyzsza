import numpy as np
import matplotlib.pyplot as plt

a = 6378137
e2 = 0.00669437999013

flightphi = []
flightlamb = []
flighth = []
with open("danelotu.txt", "r") as file:
    for line in file:
            phia, lamba, ha = line.split()
            flightphi.append(float(phia))
            flightlamb.append(float(lamba))
            flighth.append(float(ha))
    file.close()
    pass

#wspolrzedne lotniska
ap_phi = 52.4265
ap_lamb = 16.8012
ap_h = 328

def gtoxyz(phi, lamb, h):
    phi = np.deg2rad(phi)
    lamb = np.deg2rad(lamb)
    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    x = (N + h) * np.cos(phi) * np.cos(lamb)
    y = (N + h) * np.cos(phi) * np.sin(lamb)
    z = (N * (1 - e2) + h) * np.sin(phi)
    return x, y, z

def gtoneu(phi1, phi2, lamb1, lamb2, h1, h2):
    i = np.array(gtoxyz(phi1, lamb1, h1))
    j = np.array(gtoxyz(phi2, lamb2, h2))
    phi = np.deg2rad(phi1)
    lamb = np.deg2rad(lamb1)
    R = np.array([[-np.sin(phi) * np.cos(lamb), -np.sin(lamb), np.cos(phi) * np.cos(lamb)],
                  [-np.sin(phi) * np.sin(lamb), np.cos(lamb), np.cos(phi) * np.sin(lamb)],
                  [np.cos(phi), 0, np.sin(phi)]])

    R = np.transpose(R)
    delta = np.array([j[0] - i[0],
                      j[1] - i[1],
                      j[2] - i[2]])
    neu = R @ delta
    return neu

plt.plot(flightlamb, flightphi)
plt.show()

def azimuth(n, e):
    for i in range (len(n)):
        Az = np.arctan(e/n)
        Az = np.rad2deg(Az)
        if e[i] > 0 and n[i] > 0:
            Az = Az
        elif e[i] > 0 and n[i] < 0:
            Az = Az + 180
        elif e[i] < 0 and n[i] < 0:
            Az = Az + 180
        elif e[i] < 0 and n[i] > 0:
            Az = Az + 360
    return Az

def slant_distance(n, e, u):
    s = np.sqrt(n ** 2 + e ** 2 + u ** 2)
    return s

def zenith_distance(n, e, u):
    zen = np.arccos(u / np.sqrt(n ** 2 + e ** 2 + u ** 2))
    zen = np.rad2deg(zen)
    return zen
n, e, u = gtoneu(ap_phi, flightphi, ap_lamb, flightlamb, ap_h, flighth)

zen = zenith_distance(n, e, u)
zen2 = zen[zen > 90]
overhorizon = zen2[0]
index = np.where(zen == overhorizon)
nz = n[index]
ez = e[index]
uz = u[index]
print(nz,ez,uz)

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(flightlamb, flightphi, flighth, color="b", s=40)
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(218101.25473561,-310604.53269594, -328.94546456, color='r', s=90)
ax1.text(218101.25473561,-310604.53269594, -328.94546456, "Samolot znika za horyzontem")
ax1.scatter(e, n, u, color='b', s=40)
plt.show()

def dectodeg(x):
    degrees = int(x)
    minutes = int((x - degrees) * 60)
    seconds = (x - degrees - minutes / 60) * 3600
    seconds = round(seconds, 5)
    return degrees,minutes,seconds

print(dectodeg(zen[-1]))
az = azimuth(n,e)
sd = slant_distance(n,e,u)
print(dectodeg(az[-1]))
print(sd[-1])
