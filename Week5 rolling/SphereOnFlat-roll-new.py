#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Sphere rolling on a surface
@started: yuanshenli Thu 30Jan2020 to check Montana eq (40-43)
3Nov2021 cutkosky - Modified heavily for sphere (obj1) rolling on a plane (obj2).
"""

from sympy import sin, cos, tan, Symbol, Matrix, diag, eye, zeros
from sympy import simplify
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint

"""
Montana: 'We describe the motion of obj1 relative to obj2
using local coordinate frames Cl1(t) and CL2(t).
Let vx,vy,vz be the components of translational velocity
of Cl1(t) relative to CL2(t) at time t...''
Montana: obj1 rolls on obj2
"""

R1 = Symbol('R1', positive=True, real=True)
u1 = Symbol('u1', real=True)
v1 = Symbol('v1', real=True)

R2 = Symbol('R2', positive=True, real=True)
u2 = Symbol('u2', real=True)
v2 = Symbol('v2', real=True)

"""
What happens if obj1 is a sphere, obj 2 is flat?
"""
# KTM for unit sphere: R=1, [K] = [I], Torsion, Metric as in
# Montana eq. (15). See MontanaKmat.py for derivation
Kmat1 = eye(2)
Tmat1 = Matrix([0, tan(u2)]).transpose()
Mmat1 = diag(1, cos(u2))

# Make obj2 be a flat surface
# KTM for flat surface: no curvature, no torsion, Metric = [I]
# One can also derive  these using fvec = Matrix([u,v,w]) and
# following process of MontanaKmat.py
Kmat2 = zeros(2)
Tmat2 = zeros(1, 2)
Mmat2 = eye(2)

# Now define symbols associated with relative velocity of obj1 on obj2
psi = Symbol('psi', real=True)  # alignment of frames on obj1 and obj2
omegax = Symbol('omegax', real=True)
omegay = Symbol('omegay', real=True)
omegaz = Symbol('omegaz', real=True)  # could be 0 if 'soft finger'
vx = Symbol('vx', real=True)  # 0 if no sliding
vy = Symbol('vy', real=True)  # 0 if no sliding
#vz = 0

"""
Set up Montana equations (16-20) for iteration as we roll
"""
# Eq 16
Rpsi = Matrix([[cos(psi), -sin(psi)], [-sin(psi), -cos(psi)]])
K2_tilde = Rpsi*Kmat2*Rpsi

# Eq 17
Krel = Kmat1 + K2_tilde
v1gen = Matrix([-omegay, omegax]) - K2_tilde * Matrix([vx, vy])
du1 = simplify(Mmat1.inv() * Krel.inv() * v1gen)

# Eq 18
v2gen = Matrix([-omegay, omegax]) + Kmat1 * Matrix([vx, vy])
du2 = simplify(Mmat2.inv() * Rpsi * Krel.inv() * v2gen)


# Eq 19
dpsi = omegaz + (Tmat1 * Mmat1 * du1 + Tmat2 * Mmat2 * du2)[0]

# If there is rolling with no relative sliding then
# vx=0, vy =0 (but note that dpsi can still be nonzero)
# Eq 42
du1 = du1.subs([(vx, 0), (vy, 0)])
du2 = du2.subs([(vx, 0), (vy, 0)])
dpsi = dpsi.subs([(vx, 0), (vy, 0)])

pprint(du1)
pprint(du2)
pprint(dpsi)

"""
Now consider what happens if we roll with some velocities about
the (moving) x and y axes. Create variables with suffix  '_t'
that will represent the current value of these quantities as a 
function of time.
"""


# Specify some rolling velocities about x and y
numsteps = 100
stepsize = 0.1
omegax_t = 1.0
omegay_t = 0.0
omegaz_t = 2.0

# These variables will hold the current values of du1, du2, dpsi
du1_t = du1.subs([(omegax, omegax_t), (omegay, omegay_t), (omegaz, omegaz_t)])
du2_t = du2.subs([(omegax, omegax_t), (omegay, omegay_t), (omegaz, omegaz_t)])
dpsi_t = dpsi.subs(
    [(omegax, omegax_t), (omegay, omegay_t), (omegaz, omegaz_t)])


# Suppose we start at origin on a plane and no angular misalignment:
u1_t = 0
v1_t = 0
u2_t = 0
v2_t = 0
psi_t = 0

# Iterate, updating du1, du2, dpsi with each step and incrementing
# the vectors u1 and u2 on the two bodies and their relative angle, psi.
# The points we want are [u2,v2] on obj2 (the flat surface). You can
# look at [u1,v1] as well, but they are harder to interpret.
plotpts = np.zeros((numsteps, 3))
for count in range(numsteps):
    du1_new = du1_t.subs(
        [(u1, u1_t), (v1, v1_t), (u2, u2_t), (v2, v2_t), (psi, psi_t)])
    du2_new = du2_t.subs(
        [(u1, u1_t), (v1, v1_t), (u2, u2_t), (v2, v2_t), (psi, psi_t)])
    dpsi_new = dpsi_t.subs(
        [(u1, u1_t), (v1, v1_t), (u2, u2_t), (v2, v2_t), (psi, psi_t)])
    u1_t += stepsize * du1_new[0]
    v1_t += stepsize * du1_new[1]
    u2_t += stepsize * du2_new[0]
    v2_t += stepsize * du2_new[1]
    psi_t += stepsize * dpsi_new
    plotpts[count, 0] = u2_t
    plotpts[count, 1] = v2_t
    plotpts[count, 2] = psi_t

plt.figure(1)
fig1 = plt.gcf()
plt.axes().set_aspect('equal', 'datalim')  # square and limited by data
plt.title('trajectory u2, v2')
plt.plot(plotpts[:, 0], plotpts[:, 1], color='b')
plt.scatter(plotpts[0, 0], plotpts[0, 1], marker='o')
plt.scatter(plotpts[numsteps-1, 0], plotpts[numsteps-1, 1], marker='X')
plt.show()

"""
u1 and v1 are the final coordinates on the sphere; u2 and v2 are the coordinates
on the plane
"""
print("omegax", omegax_t, "  omegay", omegay_t, "  omegaz", omegaz_t)
print('u1 %5.3f \t v1 %5.3f' % (u1_t, v1_t))
print('u2 %5.3f \t v2 %5.3f' % (u2_t, v2_t))
