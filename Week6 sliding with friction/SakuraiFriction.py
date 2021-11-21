#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Using the method of Hiroshi Sakurai (1990 PhD MIT) to solve for the motion
and limiting force of a rigid plate held down by N strap clamps with Coulomb
friction. As in Peshkin, etc. he treats it as a Maximum Work or Power problem. (See Goyal
& Ruina for a general discussion.). Numbers and definitions from Sakurai thesis
(MIT 1990) pp.87-100. This script is needs some minor additions at the end.
@author: cutkosky 20Jan2020; Minor updates 12Nov2021
"""
from scipy.optimize import linprog
from WrenchUtils import PTrans, Rcross
import numpy as np
from sympy import symbols, Matrix, latex
from pprint import pprint
from math import pi, sin, cos

fcx, fcy, mcz = symbols('fcx,fcy,mcz', real=True)

f1x, f1y = symbols('f1x,f1y', real=True)
f2x, f2y = symbols('f2x,f2y', real=True)
f3x, f3y = symbols('f3x,f3y', real=True)

p1x, p1y = symbols('p1x,p1y', real=True)
p2x, p2y = symbols('p2x,p2y', real=True)
p3x, p3y = symbols('p3x,p3y', real=True)


# With the following definitions you can match algebra in Sakurai eq (4.2.9), p. 87
Wrenchmat = Matrix([[1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1],
                   [-p1y, p1x, -p2y, p2x, -p3y, p3x]])

Fpoints = Matrix([f1x, f1y, f2x, f2y, f3x, f3y])

Fext = Matrix([fcx, fcy, mcz])

wext = Fext.normalized()

wextcross = Matrix([[0, -wext[2], wext[1]],
                    [wext[2], 0, -wext[0]],
                    [-wext[1], wext[0], 0]])
# Could take wextcross * Wrenchmat * Fpoints  = 0  but it's a mess...

"""
Now enter numbers per example in Sakurai p. 92.
We don't know magnitude of Fext, but we know its direction
and point of application. See slides in Week6 presentation for details.
Let angle of pulling be measured theta anticlockwise from horizontal.
Let point of force application be pe = (px, py).
"""


# Attachment is at (px,py) and rotated theta to pull:
px, py, theta = 0.0, 0, np.pi/2
fext = np.array([1, 0, 0])  # initially pulling in the vertical direction
Jb = PTrans(px, py, theta)
Jbt = np.transpose(Jb)
extwrench = Jbt.dot(fext)
unitwrench = extwrench/np.linalg.norm(extwrench)

np.set_printoptions(precision=3, suppress=1)
pprint(unitwrench)  # check that it looks OK

# Enter the coordinates for the contact points p1, p2, p3
W = Wrenchmat.subs({p1x: -2, p1y: -1, p2x: 2, p2y: -1, p3x: 0, p3y: 1})

# Convert to np.array
Wmat = np.array(W).astype(np.float64)
pprint(Wmat)

# Per Sakurai eq (4.2.20), for equilibrium we require
#  unitwrench x (Wrenchmat*Fpoints) = 0
# For linprog() we would like to convert this into form Aeq*f = beq = 0
# So we make a 3x3 skew-symmetric cross product matrix
unitcross = Rcross(unitwrench[0], unitwrench[1], unitwrench[2])

Aeq = unitcross.dot(Wmat)
beq = np.zeros(3)

f = unitwrench.dot(Wmat)   # unitwrench'*W = work in sliding

# We approximate the friction limits sqrt(fx^2+fy^2) <= 1 with
# a crude friction square. We fudge the limits a bit so the
# square fits the bounding circle a bit better on average.
# Sakurai uses 36-sided polygons.
mu = 1.0

bounds = [[-mu*0.9, mu*0.9] for i in range(6)]

sides = 36
thetas = np.linspace(0, 2*pi, sides)


def one_side(th):
    return np.array((
        (cos(th), sin(th), 0, 0, 0, 0),
        (0, 0, cos(th), sin(th), 0, 0),
        (0, 0, 0, 0, cos(th), sin(th)),
    ))


Aub = np.concatenate([one_side(th) for th in thetas], axis=0)
# print(Aub)
bub = mu*np.ones(Aub.shape[0])
# print(bub.shape)

sol = linprog(c=f, A_eq=Aeq, b_eq=beq, A_ub=Aub, b_ub=bub, method='simplex')
sol = linprog(c=-f, A_eq=Aeq, b_eq=beq, A_ub=Aub,
              b_ub=bub, method='interior-point')
sol = linprog(c=f, A_eq=Aeq, b_eq=beq, bounds=bounds, method='simplex')
print(sol)

# you can do the others...

"""
Now plug stuff into linprog() and see what happens. Your returned
array can be compared with the friction forces in Sakurai. See
slides for problem setup.

Review Week2 assignment and example (linprogexamples.py) to get
the right arguments for linprog():
scipy.optimize.linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, 
bounds=None, method='simplex', callback=None, options=None)

Note that linprog() might warn Aeq is not full rank (which happens to
be true in this example).                                                     
"""
