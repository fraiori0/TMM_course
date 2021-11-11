#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Planar grasp quality example that roughly follows the process in
Ferrari&Canny and Miller&Allen GraspIt! to compute the quality
of a grasp:
(1) Compute the wrenches for a set of frictional contacts at
    various locations on an object.
(2) Compute corresponding 3D convex hull in [fx,fy,mz] wrench space using
    using Minkowski (convex) sum of contact wrenches.
(3) Find minimum distance from the origin to enclosing convex hull.
    This is the worst-case wrench for this grasp (bigger is better).

Notes:
* We scale moments and forces equally 1:1 but you might want a factor
  corresponding to object diameter.
* First check that the contacts satisfy Force Closure.
* We take all contact forces as unit magnitude. Alternative would be
  they sum to some value while satisfying equilibrium constraints.
* This script is meant to be stepped through, a few lines at a time,
  looking at results after each step.
  
@Started: cutkosky 26Dec2019; Added 2D projections to check that
  origin inside hull -- Rachel Thomasson 28Jan2020;
  Made separate file for Minkowski Sum metric --Cutkosky1Feb2020
"""

import numpy as np
from scipy.spatial import ConvexHull
from sympy import Plane, Point3D

from WrenchUtils import PTrans     #in local directory


#Convex sum of two arrays of column vectors (e.g. wrenches, twists)
def ConvexSum(wrenches1,wrenches2):     
    dim1,n1 = wrenches1.shape
    dim2,n2 = wrenches2.shape
    if(dim2 != dim1):
        raise Exception('wrenches should have same dimension')
    
    csum = np.zeros((dim1,n1*n2))
    for i in range(n1):
        for j in range(n2):
            csum[:,i*n2+j] = wrenches1[:,i]+wrenches2[:,j]
            
    return csum



'''
Example: suppose the object is a Trapezoid with vertices at
[1,-1,0], [3,1,0], [-3,1,0], [-1,-1,0]
and suppose we start with 3 fingers, one at midpoints of the bottom and sides.
A pictorial representation of this can be found in the assignment document.
'''

#TODO: Change this value to 4 in Question 3.2
n=3  # 3 contacts

#TODO: Modify this value in Question 3.1
mu = 0.5

phi = np.arctan(mu)  # friction cone half-angle


#1.1 Create coordinate frame for each contact, oriented so that
# local X axis is along outward normal. Start at left, proceed anticlockwise.
# Again, you can reference the figure in the assignment doc

#TODO: Change these frames to place the 4 fingers in Question 3.2
frames = np.zeros((n,3))
frames[0,:] = [-2,0,5*np.pi/4]
frames[1,:] = [0,-1,-np.pi/2]
frames[2,:] = [2,0,-np.pi/4]
#frames[3,:] = [...]


#1.2 Get local contact wrenches and corresponding global wrenches
#Inward forces along left, right edges of a friction cone,
#assuming coordinate frame with X axis pointing outward
#and unit normal force along -X. So the 2 vectors are: 
fl = np.array([-np.cos(phi), -np.sin(phi), 0])
fr = np.array([-np.cos(phi),np.sin(phi), 0])

#Get the set W1 corresponding to the first contact
W1 = np.zeros((3,3))
Jb = PTrans(frames[0,0],frames[0,1],frames[0,2]) 
Jbtrans = Jb.transpose()
W1[0] = Jbtrans.dot(fl)
W1[1] = Jbtrans.dot(fr)
W1[2] = [0,0,0]

#Get the set W2 corresponding to the second contact
W2 = np.zeros((3,3))
Jb = PTrans(frames[1,0],frames[1,1],frames[1,2]) 
Jbtrans = Jb.transpose()
W2[0] = Jbtrans.dot(fl)
W2[1] = Jbtrans.dot(fr)
W2[2] = [0,0,0]

#Get the set W3 corresponding to the third contact
W3 = np.zeros((3,3))
Jb = PTrans(frames[2,0],frames[2,1],frames[2,2]) 
Jbtrans = Jb.transpose()
W3[0] = Jbtrans.dot(fl)
W3[1] = Jbtrans.dot(fr)
W3[2] = [0,0,0]

#Repeat as needed for additional fingers

wrenches = np.concatenate((W1,W2,W3))

#check: Planar wrenches array had better have rank 3
#(Remember, force closure is necessary, but not sufficient.)
np.linalg.matrix_rank(wrenches)
#For Force Closure we could also check if there is a
#solution with normal forces all being positive inward
#per Mason & Salisbury. But the metric below will more or 
#less catch this as well.

###############Now compute the Minkowski sums#################


csum12 = ConvexSum(W1.transpose(),W2.transpose())
csum123 = ConvexSum(csum12,W3.transpose())
#csum1234 = ConveSum(csum123,W4.transpose())   #etc.

msumwrenches = csum123.transpose()


#2. Assuming  wrench matrix is fine, get Convex Hull of
#the points corresponding to wrenches in (fx,fy,mz)
hull = ConvexHull(msumwrenches)

#Check if the convex hull vertices look right
np.set_printoptions(precision=2)
for s in hull.vertices:
    print(s,msumwrenches[s,:],'mag:', "%.2f" % np.linalg.norm(msumwrenches[s,:]))


#3. Find distance from origin to a plane that contains
#each triangular facet of the convex hull
origin = Point3D(0,0,0)
nfacets = np.shape(hull.simplices)[0] #how many facets
hullwrenchmags = np.zeros(nfacets)
i=0
for s in hull.simplices:
    triangle = msumwrenches[s]   #convert from simplices to 3D points
    point1 = Point3D(triangle[0])
    point2 = Point3D(triangle[1])
    point3 = Point3D(triangle[2])
    theplane = Plane(point1,point2,point3)
    planedistance =  theplane.distance(origin)
    hullwrenchmags[i]=planedistance
    i=i+1

leastwrench = np.amin(hullwrenchmags)
print('least wrench (if enclosing), Minkowski hull:',"%.2f" % leastwrench)

#The above distance calculation assumes the convex hull encloses
#the origin. We should check to be sure that is true! An easy
#way is to plot orthogonal projections.

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig3D = plt.figure()
ax3D = fig3D.add_subplot(111, projection='3d')

figxz = plt.figure()
axxz = figxz.add_subplot(111, projection='3d')

figyz = plt.figure()
axyz = figyz.add_subplot(111, projection='3d')

figxy = plt.figure()
axxy = figxy.add_subplot(111, projection='3d')

# For those points in hull.simplices...
for s in hull.simplices:
    s = np.append(s, s[0])  # Cycle back to the first coordinate
    ax3D.plot(msumwrenches[s, 0], msumwrenches[s, 1], msumwrenches[s, 2], "r-")
    ax3D.scatter(msumwrenches[s, 0], msumwrenches[s, 1], msumwrenches[s, 2], marker='o')

    axxz.plot(msumwrenches[s, 0], msumwrenches[s, 2], msumwrenches[s, 1], "r-")
    axxz.scatter(msumwrenches[s, 0], msumwrenches[s, 2], msumwrenches[s, 1], marker='o')

    axyz.plot(msumwrenches[s, 1], msumwrenches[s, 2], msumwrenches[s, 0], "r-")
    axyz.scatter(msumwrenches[s, 1], msumwrenches[s, 2], msumwrenches[s, 0], marker='o')

    axxy.plot(msumwrenches[s, 0], msumwrenches[s, 1], msumwrenches[s, 2], "r-")
    axxy.scatter(msumwrenches[s, 0], msumwrenches[s, 1], msumwrenches[s, 2], marker='o')

    # also plot the origin
    ax3D.scatter([0], [0], [0], marker='x')
    axxz.scatter([0], [0], [0], marker='x')
    axyz.scatter([0], [0], [0], marker='x')
    axxy.scatter([0], [0], [0], marker='x')


axxz.view_init(azim=0, elev=90)
axyz.view_init(azim=0, elev=90)
axxy.view_init(azim=0, elev=90)

# Make axis label
for i in ["x", "y", "z"]:
    eval("ax3D.set_{:s}label('{:s}')".format(i, i))

plt.show()