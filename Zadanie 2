import numpy as np
import math
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def juliandate(y,m,d,h):
    if (m <= 2):
       y = y -1
       m = m + 12
    jd = math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + d + h / 24 - 1537.5
    return jd

#Greenwich Mean Sidereal Time
def GMST(y,m,d,h):
    T = (juliandate(y,m,d,h) - 2451545) / 36525
    g = 280.46061837 + 360.98564736629 * (juliandate(y,m,d,0) -2451545) + 0.000387933 * (T ** 2) - (T ** 3) / 38710000
    g = np.remainder(g, 360)
    return g

#zamiana jednostki rektastensji na radiany
a = 21.09911884
right_ascension = np.deg2rad(a)

#zamiana deklinacji na radiany
b = -17.232862
declination = np.deg2rad(b)

#obliczenie kata godzinnego
def hour_angle(y,d,m,h,city_lambd,right_as):
    g = GMST(y,m,d,0)
    UT1 = h * 1.002737909350795
    S = UT1 * 15 + city_lambd + g
    t = S - np.rad2deg(right_as) * 15
    return t

#rozwiazanie trojkata paralaktycznego
def zenith_distance(city_phi,dec,t):
    #zamiana szerokosci geograficznej miasta na radiany i kata godz
    city_phi = np.deg2rad(city_phi)
    t = np.deg2rad(t)
    cos_z = (math.sin(city_phi) * math.sin(dec)) + (math.cos(city_phi) * math.cos(dec) * math.cos(t))
    cos_z = np.rad2deg(np.arccos(cos_z))
    return cos_z

def star_azimuth(dec,t,city_phi):
    city_phi = np.deg2rad(city_phi)
    t = np.deg2rad(t)
    num = (-(math.cos(dec)) * math.sin(t))
    denum = (math.cos(city_phi) * math.sin(dec)) - (math.sin(city_phi) * math.cos(dec) * math.cos(t))
    tg_A = np.arctan(num/denum)
    tg_A = np.rad2deg(tg_A)
    if num > 0 and denum > 0:
        tg_A = tg_A
    elif num > 0 and denum < 0:
        tg_A = tg_A + 180
    elif num < 0 and denum < 0:
        tg_A = tg_A + 180
    elif num < 0 and denum > 0:
        tg_A = tg_A + 360
    return tg_A

#transformacja wspolrzednych
def transform_coordinates(cos_z,tg_A):
    cos_z = np.deg2rad(cos_z)
    tg_A = np.deg2rad(tg_A)
    x = 1 * math.sin(cos_z) * math.cos(tg_A)
    y = 1 * math.sin(cos_z) * math.sin(tg_A)
    z = 1 * math.cos(cos_z)
    x = np.rad2deg(x)
    y = np.rad2deg(y)
    z = np.rad2deg(z)
    return x,y,z

def height(phi, dec, t):
    h = 90 - zenith_distance(phi, dec, t)
    return h

#glasgow
tgA_array = []
cosz_array = []
hglas_array = []
city_lambd = -4.2576300
city_phi = 55.8651500
for h in range(24):
    hA = hour_angle(2021,21,11,h,city_lambd,right_ascension)
    cos_z = zenith_distance(city_phi,declination,hA)
    tg_A = star_azimuth(declination,hA,city_phi)
    h_glasgow = height(city_phi, declination, hA)
    tgA_array.append(tg_A)
    cosz_array.append(cos_z)
    hglas_array.append(h_glasgow)

hours = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
hours = np.asarray(hours)
plt.plot(hours, cosz_array)
plt.xlabel('Czas')
plt.ylabel('Odległość zenitalna')
plt.title('Odległość zenitalna gwiazdy w zależności od czasu (Glasgow)')
plt.show()

x_array = []
y_array = []
z_array = []
for i in range(24):
    x,y,z = transform_coordinates(cosz_array[i],tgA_array[i])
    x_array.append(x)
    y_array.append(y)
    z_array.append(z)



plt.plot(hours, hglas_array)
plt.xlabel('Czas')
plt.ylabel('Wysokość gwiazdy nad horyzontem')
plt.title('Wysokość gwiazdy nad horyzontem w zależności od czasu (Glasgow)')
plt.show()

r = 56
phi, theta = np.mgrid[0:np.pi:30j, 0:2*np.pi:30j]
xs = r*np.sin(phi)*np.cos(theta)
ys = r*np.sin(phi)*np.sin(theta)
zs = r*np.cos(phi)


fig = pyplot.figure()
ax = Axes3D(fig)
ax.scatter(x_array,y_array,z_array)
ax.plot_surface(xs, ys, zs, alpha=0.2)
plt.show()

#melbourne
tgA_array2 = []
cosz_array2 = []
hmelb_array = []
city_lambd2 = 144.9633200
city_phi2 = -37.8140000
for h in range(24):
    hA2 = hour_angle(2021,21,11,h,city_lambd2,right_ascension)
    cos_z2 = zenith_distance(city_phi2,declination,hA2)
    tg_A2 = star_azimuth(declination,hA2,city_phi2)
    h_melbourne = height(city_phi2, declination, hA2)
    tgA_array2.append(tg_A2)
    cosz_array2.append(cos_z2)
    hmelb_array.append(h_melbourne)

plt.plot(hours, cosz_array2)
plt.xlabel('Czas')
plt.ylabel('Odległość zenitalna')
plt.title('Odległość zenitalna gwiazdy w zależności od czasu (Melbourne)')
plt.show()

x_array2 = []
y_array2 = []
z_array2 = []
for i in range(24):
    x,y,z = transform_coordinates(cosz_array2[i],tgA_array2[i])
    x_array2.append(x)
    y_array2.append(y)
    z_array2.append(z)


h_melbourne = height(city_phi2, declination, hA2)
plt.plot(hours, hmelb_array)
plt.xlabel('Czas')
plt.ylabel('Wysokość gwiazdy nad horyzontem')
plt.title('Wysokość gwiazdy nad horyzontem w zależności od czasu (Melbourne)')
plt.show()

r = 56
phi, theta = np.mgrid[0:np.pi:30j, 0:2*np.pi:30j]
xs = r*np.sin(phi)*np.cos(theta)
ys = r*np.sin(phi)*np.sin(theta)
zs = r*np.cos(phi)

fig = pyplot.figure()
ax = Axes3D(fig)
ax.scatter(x_array2,y_array2,z_array2)
ax.plot_surface(xs, ys, zs, alpha=0.2)
plt.show()

#singapore
tgA_array3 = []
cosz_array3 = []
hsing_array = []
city_lambd3 = 103.8500700
city_phi3 = 1.2896700
for h in range(24):
    hA3 = hour_angle(2021,21,11,h,city_lambd3,right_ascension)
    cos_z3 = zenith_distance(city_phi3,declination,hA3)
    tg_A3 = star_azimuth(declination,hA3,city_phi3)
    h_singapore = height(city_phi3, declination, hA3)
    tgA_array3.append(tg_A3)
    cosz_array3.append(cos_z3)
    hsing_array.append(h_singapore)

plt.plot(hours, cosz_array3)
plt.xlabel('Czas')
plt.ylabel('Odległość zenitalna')
plt.title('Odległość zenitalna gwiazdy w zależności od czasu (Singapur)')
plt.show()

x_array3 = []
y_array3 = []
z_array3 = []
for i in range(24):
    x,y,z = transform_coordinates(cosz_array3[i],tgA_array3[i])
    x_array3.append(x)
    y_array3.append(y)
    z_array3.append(z)


h_singapore = height(city_phi3, declination, hA3)
plt.plot(hours, hsing_array)
plt.xlabel('Czas')
plt.ylabel('Wysokość gwiazdy nad horyzontem')
plt.title('Wysokość gwiazdy nad horyzontem w zależności od czasu (Singapur)')
plt.show()

r = 56
phi, theta = np.mgrid[0:np.pi:30j, 0:2*np.pi:30j]
xs = r*np.sin(phi)*np.cos(theta)
ys = r*np.sin(phi)*np.sin(theta)
zs = r*np.cos(phi)

fig = pyplot.figure()
ax = Axes3D(fig)
ax.scatter(x_array3,y_array3,z_array3)
ax.plot_surface(xs, ys, zs, alpha=0.2)
plt.show()


