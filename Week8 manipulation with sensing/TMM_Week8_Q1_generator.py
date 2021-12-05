# -*- coding: utf-8 -*-
"""
Create a sequence of wrenches at the center of a hemisphere associated
with a moving contact. In this case we assume that it is a point
contact sliding with friction.

We use the ISO standard spherical coordinates per Wikipedia 
 https://en.wikipedia.org/wiki/Spherical_coordinate_system 
theta = polar angle (down from Z), phi = azimuth (anticlockwise from X)
Started 29Feb 2020 - Cutkosky
3March - Wilson Ruotolo additions
25Nov2021 - Modifications for TMM 2021 cutkosky
"""

from pprint import pprint
import numpy as np
from WrenchUtils import Rotx, Roty, Rotz, Rcross  # from week2
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)

"""
First, a few utility functions that may be useful for converting between
spherical and Cartesian coordinates and transforming contact forces
to body wrenches.
"""
# Convert from Cartesian to spherical coordinates
# [x,y,z] to [r,phi,theta] where phi=azimuth measured from X,
# theta=inclination measured down from Z (radians)


def xyz2rphitheta(point):
    x, y, z = point[0], point[1], point[2]
    r = (x**2+y**2+z**2)**0.5
    theta = np.arctan2(np.sqrt(x**2+y**2), z)
    phi = np.arctan2(y, x)
    spherecoords = np.zeros(3)
    spherecoords[0], spherecoords[1], spherecoords[2] = r, phi, theta
    return spherecoords

# Convert from spherical to Cartesian coordinates
# [r,phi,theta] to [x,y,z] where phi=azimuth measured from X,
# theta=inclination from Z (radians)


def rphitheta2xyz(spherecoords):
    r, phi, theta = spherecoords[0], spherecoords[1], spherecoords[2]
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    cpoint = np.zeros(3)
    cpoint[0], cpoint[1], cpoint[2] = x, y, z
    return cpoint

# Transformation between two frames, the second being
# translated by [rx,ry,rz] and then rotated by [thetaz,thetay,thetax].
# In contrast to the RPY definition in WrenchUtils.py,
# for the ISO spherical coordinates it is convenient to use
# a mapping where we rotate first about Z, then Y, then X


def CartesmapZYX(rx, ry, rz, thetaz, thetay, thetax):
    Rx = Rotx(thetax)
    Ry = Roty(thetay)
    Rz = Rotz(thetaz)
    Amat = Rz @ Ry @ Rx
    Rskew = Rcross(rx, ry, rz)
    UL = np.transpose(Amat)
    UR = np.dot(UL, np.transpose(Rskew))
    LL = np.zeros((3, 3))
    LR = UL
    Jb = np.block([[UL, UR], [LL, LR]])
    return Jb


"""
Define and plot the unit radius sphere
"""
# Setup figure
fig = plt.figure(figsize=plt.figaspect(.75))  # Square figure
ax = fig.add_subplot(111, projection='3d')

# Draw a unit sphere using numpy.mgrid() and surface plot
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi/2:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_surface(x, y, z, color="r", alpha=.5)

# Make axis labels
for i in ["x", "y", "z"]:
    eval("ax.set_{:s}label('{:s}')".format(i, i))
ax.set_ylim(1, -1)

"""
Some tests to make sure utilities are working corectly
"""
# Try an example at the north pole
#spherept = np.array([1,0,0])
#xyzpt = rphitheta2xyz(spherept)
# newspherept = xyz2rphitheta(xyzpt)  #Should be same as spherept

# another example, half way down and on -x side
#spherept = np.array([1,np.pi,np.pi/4])
#xyzpt = rphitheta2xyz(spherept)
#newspherept = xyz2rphitheta(xyzpt)

#Jb = CartesmapZYX(xyzpt[0],xyzpt[1],xyzpt[2],0,0,0)
#Jbtran = Jb.transpose()
#wrench = Jbtran @ fcontact


"""
Create a sequence of wrenches measured at the center of a sphere,
corresponding to a trajectory of contact forces for a sliding 
contact with friction. 
(With no friction we do not get any moments at the center
of the sphere.)
Also create a second version of the wrenches where we add some
noise, as would be measured by a real sensor.
Write the results to two text files.
"""
thetas = np.arange(0, (np.pi/2), np.pi/20)
phis = np.arange(0, ((np.pi/2)+np.pi/20), np.pi/20)
wrenches = np.zeros((thetas.size, 6))
r = 1.0
scoords = np.zeros(3)
xyz = np.zeros(3)
plotvec = np.zeros((thetas.size, 6))

"""
For a force trajectory we take a contact force (assumed constant
in local [x,y,z] coordinates in this example) and map it to
the corresponding series of wrenches at the center of the sphere
while traversing a sequence of angles (phi,theta).
plotvec puts each angled fcontact in world coordinates for potting
"""
fcontact = np.array([-0.5, 0.5, -1, 0, 0, 0])
for i in range(thetas.size):  # need phis[] to have same size as thetas[]
    scoords[0], scoords[1], scoords[2] = r, thetas[i], phis[i]
    xyz = rphitheta2xyz(scoords)
    Jb1 = CartesmapZYX(0, 0, 0, phis[i], thetas[i], 0)
    Jb2 = CartesmapZYX(0, 0, r, 0, 0, 0)
    # ideally i'd rewrite CartesmapZYX to do rotation, then translation -mrc
    Jbtotal = Jb2 @ Jb1
    Jbtran = Jbtotal.transpose()
    wrenches[i] = Jbtran @ fcontact
    plotvec[i] = np.array(
        [xyz[0], xyz[1], xyz[2], wrenches[i, 0], wrenches[i, 1], wrenches[i, 2]])

pprint(plotvec[:, :3])
# Make a new version with 5% normally distributed noise
np.random.seed(123)
noiselevel = 0.05*np.linalg.norm(fcontact)
noise = np.random.normal(0, noiselevel, wrenches.shape)
noisywrenches = wrenches + noise


# Plot the sequence of contact forces on the sphere as arrows
# unzip: [X,Y,Z] has point on sphere; [U,V,W] has force components
X, Y, Z, U, V, W = zip(*plotvec)
ax.quiver(X, Y, Z, U, V, W, pivot='tip', length=.2, color='y')

plt.show()
plt.savefig("Week8 manipulation with sensing/intrinsic-input.png", dpi=600)


"""
Save the wrenches data to text files, tab delimited.
"""
headerstring = "Wrenches [fx,fy.fz,mx,my,mz] for trajectory of contact forces on sphere"
np.savetxt('Week8 manipulation with sensing/wrench-sequence.txt', wrenches, header=headerstring,
           delimiter='\t', newline='\n', fmt='%5.3f')

headerstring = "Wrenches [fy,fz,mx,my,mz] + noise for trajectory of contact forces on sphere"
np.savetxt('Week8 manipulation with sensing/noisywrench-sequence.txt', noisywrenches,
           header=headerstring, delimiter='\t', newline='\n', fmt='%5.3f')

"""
The saved files are the result of running the 'forward' calculation. But we need to solve the
inverse problem: given a series of wrenches (e.g. as in the files above) deduce what
the sequence of contact forces must have been, using Intrinsic Tactile Sensing.
"""
