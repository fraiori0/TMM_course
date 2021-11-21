# -*- coding: utf-8 -*-
"""
Limit surface method on a planar sliding object with discrete
frictional contacts. Adapted from Howe & Cutkosky, "Practical
Force-Motion Models for Sliding Manipulation," IJRR 1996. The script 
is written so it's easy to change contacts or add new ones.

Initial numbers are from H. Sakurai "Automatic Setup Planning and Fixture 
Design for Machining", PhD thesis, MIT, 1990, p93. He gives an example with 
3 contact points and an external force that is like pulling on a string 
with unknown magnitude but fixed ratios of forces.
See 'SakuraiExampleSlides.pdf' slides in the Week6 folder for setup.

We create an array for the sliding contacts with each row as [x,y,mu*fn] 
giving the location and friction force of each point.
The external force has a given pulling direction and contact location, but
unkown magnitude. So it creates a unit wrench [fx,fy,mz] about the origin
which will be multiplied by some scaling factor, rho, to be found.

1. First, if the number of contacts is small (like 4 or fewer) there 
is a good chance it will rotate about one of the contacts. This corresponds
to a facet on the friction limit surface. In this case, with an assumed COR,
we can solve for the force and moment with linear equations (see slides
in Week6 presentation).

2. Otherwise we fit an approximating ellipsoid using selected possible
COR locations. Any combination of rho*[fx,fy,mz] can then be tested against
this ellipsoid to check for slipping and associated unit twist [vx,vy,omega]

26Feb2020 MRC: Fixed typo in LSwrench() line 59 which makes the moment
too high (unless mu*fn=1)
18Nov2021 -mrc Additions for sliding twist at bottom and 
to consider if rotation is clockwise.
"""

import matplotlib.pyplot as plt
import numpy as np
from WrenchUtils import PTrans

'''
Define a couple of useful utilities:
'''
# Given an array of sliding points and a COR location (rcx,rcy),
# compute total friction wrench w.r.t. origin, assuming
# anticlockwise rotation. This gives one point on a LS.
# Each row of contacts[,,] should have [x,y,mu*fn] for a contact point.
# rotation < 0 means clockwise; else anticlockwise (default)


def LSwrench(rcx, rcy, contacts, rotation):
    npoints = np.shape(contacts)[0]
    fwrench = np.zeros(3)
    tiny = 1.0e-8  # in case COR is essentially at this point
    if (rotation < 0):
        rotation = -1  # default
    else:
        rotation = 1

    for i in range(npoints):
        mufn = contacts[i, 2]
        rx = contacts[i, 0]-rcx
        ry = contacts[i, 1]-rcy
        rmag = np.sqrt(rx**2+ry**2)
        if rmag > tiny:
            fx = mufn*(ry/rmag)
            fy = mufn*(-rx/rmag)
        else:
            fx, fy = 0, 0

        mz = (-fx*contacts[i][1]+fy*contacts[i][0])
        fwrench[0] += fx
        fwrench[1] += fy
        fwrench[2] += mz

    return rotation * fwrench

# friction-weighted center of pressure of a set of contacts


def Centroid(contacts):
    center = np.zeros(2)
    npoints = np.shape(contacts)[0]
    cx, cy, fmag = 0, 0, 0
    for i in range(npoints):
        fmag += contacts[i, 2]
        cx += contacts[i, 0]*contacts[i, 2]
        cy += contacts[i, 1]*contacts[i, 2]
    center[0] = cx/fmag
    center[1] = cy/fmag

    return center


"""
PART ONE 
Could it be rotating about one of the contact points? 
If so (see slides) this corresponds to a facet and it's
most accurate to solve the sliding equilibrium equations directly.
    
Equilibrium equations in the plane (see slides for diagram) are:
[Qmat]*[fcx,fcy,rho]' + [frx,fry,frm]' = [0,0,0]'
Where [frx,fry,frm] are the net frictional forces and moment from
the sliding points. 
Invert Qmat to solve for the unkowns, [fcx,fcy,rho] 
"""
# Values here taken from Sakurai p. 93. (modify to suit)
# External force and where it's attached.
# We don't know magnitude yet, so use unit magnitude.
fext = np.array([0, -1, 0])
pex, pey = 2, -2

# Get the equivalent wrench [fx,fy,m] at origin and make a unit wrench.
fetrans = PTrans(pex, pey, 0)  # from WrenchUtils.py
Je = fetrans.transpose()
ewrench = Je.dot(fext)
uwrench = ewrench/np.linalg.norm(ewrench)

# Create an array of data for sliding contacts. Each row: [px, py, mu*fn]
# In this example mu*fn = 1 for all points (values from Sakurai).
mu = 0.5
fn = 1.0
p1x, p1y, ft1 = 2, 1, 2*mu*fn  # we'll see if it rotates about this one
p2x, p2y, ft2 = 2, -2, mu*fn
p3x, p3y, ft3 = -2, -2, mu*fn
p4x, p4y, ft4 = -2, 1, 2*mu*fn
# add additional points as needed...
# Put sliding points in an array (all except the CoR point)
# scontacts = np.array([[p1x, p1y, ft1], [p2x, p2y, ft2], [p3x, p3y, ft3]])
scontacts = np.array([[p1x, p1y, ft1], [p2x, p2y, ft2], [p4x, p4y, ft4]])

# The non-sliding point is our candidate COR
# rcx, rcy, ftmax = p4x, p4y, ft4
rcx, rcy, ftmax = p3x, p3y, ft3

# [frx, fry, frm] are the net friction forces and moment for all sliding
# contacts (i.e. other than the point about which we're rotating)
frwrench = LSwrench(rcx, rcy, scontacts, -1)


# Qmat * [fcx, fcy, rho] = -[frwrench] where fcx,fcy are the unknown
# force components at the center of rotation (rcx, rcy), and rho is the
# unknown scaling factor for uwrench. frwrench = [frx, fry, frm] is the
# friction force and moment w.r.t. origin from the sliding contacts for
# anticlockwise rotation about COR.
# So Qmat^-1 * [frwrench] = -[fcx,fcy,rho]
Qmat = np.array(
    [[1, 0, uwrench[0]], [0, 1, uwrench[1]], [-rcy, rcx, uwrench[2]]])
Qinv = np.linalg.inv(Qmat)  # could be singular if problem ill-posed
result = -Qinv.dot(frwrench)
fcx, fcy, rho = result[0], result[1], result[2]

# Finally compute what the pulling force would be
if np.sqrt(fcx**2+fcy**2) < ftmax:
    Jinv = np.linalg.inv(Je)
    fpull = Jinv.dot(rho*uwrench)
    print('Part 1 External pulling wrench: ', fpull)
else:
    print('Friction limit exceeded at assumed COR')
# If it fails, then our assumption of COR at this point is wrong.
# We could try one of the other points as COR...

"""
PART TWO
Now compute more generally what the limit surface looks like using
the Howe/Lee/Cutkosky ellipsoid.
    
First, compute the centroid of the pressure distribution and 
move the origin there. This will often give a more symmetric LS.
(If you do this, you may need to shift your external wrench also,
depending on your problem setup.)
"""

# Augment scontacts[] to now include previous CoR point and move origin
extrapoint = np.array([rcx, rcy, ftmax])
scontacts = np.vstack([extrapoint, scontacts])  # add in the 1st point
npts = np.shape(scontacts)[0]
origin = Centroid(scontacts)

# If you want to recenter the points about new origin:
symcontacts = np.zeros((npts, 3))
for i in range(npts):
    symcontacts[i, 0] = scontacts[i, 0] - origin[0]
    symcontacts[i, 1] = scontacts[i, 1] - origin[1]
    symcontacts[i, 2] = scontacts[i, 2]
# symcontacts = scontacts   #If don't want to recenter


# We should recenter our pulling force too. Depending on setup
# this may or may not make a difference to uwrench
newfetrans = PTrans(pex-origin[0], pey-origin[1], 0)
newJe = newfetrans.transpose()
newewrench = newJe.dot(fext)
newuwrench = newewrench/np.linalg.norm(ewrench)

"""
Next compute the apex of the ellipsoid. The maximum possible
moment should now be for a COR at origin. fx, fy will ideally
be small for this case, meaning the LS is nearly symmetric.

For convenience of plotting, we take moment to be positive.
If your problem has a negative moment and clockwise sliding
you may want to flip sign on the moment and then remember 
to flip back when computing the sliding twist at the end.
"""
# Applied moment is opposite of sliding friction moment
maxmoment = -LSwrench(0, 0, symcontacts, 1)

# If the ellipsoid is not tilted much we have its axes:
# (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
c = maxmoment[2]
a = np.sum(symcontacts[:, 2])  # sum of the tangential forces
b = a  # assuming isotropic friction

"""
Let us plot the ellipsoid in [fx, fy, mz] space to see how it looks
"""
#from mpl_toolkits.mplot3d import Axes3D

theta = np.linspace(0, np.pi, 20)
phi = np.linspace(0, 2*np.pi, 20)

x = a * np.outer(np.cos(phi), np.sin(theta))
y = b * np.outer(np.sin(phi), np.sin(theta))
z = c * np.outer(np.ones_like(phi), np.absolute(np.cos(theta)))

fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
ax = fig.add_subplot(111, projection='3d')

# Plot:
ax.plot_wireframe(x, y, z, color='y')

# Make axis labels
for i in ["x", "y", "z"]:
    eval("ax.set_{:s}label('{:s}')".format(i, i))

"""
Now we can test any applied wrench rho*[ufx,ufy, umz] to 
see if it would cause slipping (i.e., if it intersects the ellipsoid).
From the point of intersection we get the sliding twist,
[vx,vy,omega], given by the unit normal to the LS at the intersection. 
(See Table 2 in Howe & Cutkosky)
"""


""" ----------------------------------------------------------------
Example: Having constructed the ellipsoid we can use it.
"""
# Suppose we apply a new external wrench that we suspect does not cause
# rotation about one of the contact points (in which case the PART ONE
# method above would be better). For example we might have:
newwrench = np.array([2, 2, 3])


"""
Note that we have been assuming anti-clockwise rotation. If our wrench has
a negative moment we can flip the sign and remember that +z is actually
clockwise for our case
"""

# We can find where this wrench intersects
# the ellipsoid
phi = np.arctan2(newwrench[1], newwrench[0])
r = np.sqrt(newwrench[0]**2+newwrench[1]**2)
theta = np.arctan2(c*r, a*newwrench[2])

# Using ellipsoid definition:
x = a*np.cos(phi)*np.sin(theta)
y = b*np.sin(phi)*np.sin(theta)
z = c*np.absolute(np.cos(theta))
slidewrench = np.array([x, y, z])
# slidewrench should be a scaled version of newwrench that just barely
# intersects the ellipsoidal shell.

# We can plot this point
ax.scatter(x, y, z, marker='o')
print('\nPart 2 New sliding wrench (see blue dot)', slidewrench)

"""
Sliding velocity
"""
# vx,vy components of the sliding twist will be parallel to fslipx, fslipy
Lam = c/a  # ratio of ellipse axes
vx = slidewrench[0]
vy = slidewrench[1]
vtan = np.sqrt(vx**2+vy**2)
omegaz = slidewrench[2]/(Lam**2 * vtan)  # see Howe & Cutkosky Table 2.

slidetwist = np.array([vx, vy, omegaz])
mag = np.linalg.norm(slidetwist)
print("unit sliding twist:", slidetwist/mag)


plt.show()
