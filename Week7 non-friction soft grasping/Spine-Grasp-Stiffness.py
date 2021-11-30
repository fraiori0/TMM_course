#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adapting from Week 3 Grasp Stiffness materials for the case of a gripper
with spines grasping a block from above and from the right.
22Nov2021 -Cutkosky
"""

from sympy import sin, cos, pi, Matrix, Symbol, symbols, simplify, pprint, lambdify
from sys import exit

# poor man's debugging


def yoh():
    exit('Yeeet')


#################################
# Define 6 element wrench and 6 element twist
fx, fy, fz, mx, my, mz = symbols('fx,fy,fz,mx,my,mz', real=True)
dbx, dby, dbz = symbols('dbx,dby,dbz', real=True)
qbx, qby, qbz = symbols('qbx,qby,qbz', real=True)

wrench = Matrix([fx, fy, fz, mx, my, mz])
twist = Matrix([dbx, dby, dbz, qbx, qby, qbz])

###############################################################
# Jbtran is a 6x6 matrix that, when multipied by a wrench
# in a body contact frame, gives the corresponding body wrench.
# Contact is assumed to be translated by rx,ry,rz from the body
# frame and then rotated by thetax,thetay,thetaz (in order).
# Transpose of Jbtran maps a twist from the body frame to the
# corresponding contact frame.
rx, ry, rz = symbols('rx, ry, rz', real=True, positive=True)
thetax, thetay, thetaz = symbols('thetax,thetay,thetaz', real=True)

Rotx = Matrix([[1, 0, 0], [0, cos(thetax), -sin(thetax)],
              [0, sin(thetax), cos(thetax)]])
Roty = Matrix([[cos(thetay), 0, sin(thetay)], [
              0, 1, 0], [-sin(thetay), 0, cos(thetay)]])
Rotz = Matrix([[cos(thetaz), -sin(thetaz), 0],
              [sin(thetaz), cos(thetaz), 0], [0, 0, 1]])
Amat = Rotx*Roty*Rotz  # 3x3 orthonormal orientation

# cross product matrix
Rskew = Matrix([[0, -rz, ry], [rz, 0, -rx], [-ry, rx, 0]])

top = Amat.row_join(Matrix.zeros(3))
bottom = (Rskew*Amat).row_join(Amat)
Jbtran = top.col_join(bottom)

#########################################################
# Create the H matrix that filters out rotational terms at contacts
# H*twist = twist_transmitted
Hmat = (Matrix.eye(3)).row_join(Matrix.zeros(3))  # 3x6

#########################################################
# Specialize Jbtran for our object and contact coordinate frames
# contact1: translate along z by w/2 and rotate about x by pi
# contact2: translate along x by w/2 and rotate about y by -pi/2
w = Symbol('w', real=True, positive=True)
Jb1t = Jbtran.subs([(rx, 0), (ry, 0), (rz, w/2),
                   (thetax, pi), (thetay, 0), (thetaz, 0)])
Jb2t = Jbtran.subs([(rx, w/2), (ry, 0), (rz, 0),
                   (thetax, 0), (thetay, -pi/2), (thetaz, 0)])


# Build our fingertip spine stiffness matrices in local coordinates
ksl = Symbol('ksl', real=True, positive=True)
ksn = Symbol('ksn', real=True, positive=True)
Kfp1 = Matrix.diag(ksl, ksl, ksn)
Kfp2 = Matrix.diag(ksl, ksl, ksn)

"""
Now apply a small (x,z) displacement of the body. Get the corresponding
small displacements at the fingers and the corresponding forces
from the (inverted) fingertip stiffness matrices.
"""

dbody = twist.subs([(dby, 0), (qbx, 0), (qby, 0), (qbz, 0)])

dxf1 = Hmat*Jb1t.T*dbody  # Motion of the fingertip
dxf2 = Hmat*Jb2t.T*dbody

df1 = Kfp1*dxf1  # Corresponding change in force
df2 = Kfp2*dxf2

print('df1:')
pprint(df1)


print('\ndf2:')
pprint(df2)


""" 
These force vectors now can be added to a bias force, if any, 
and the sum can be compared against the limit surfaces
for the corresponding spine contacts to check for failure.

Note that a bias force must be the same for both fingers and must only
have equal X and Z components to keep it in the null space of
the grasp:  [fbias,0,fbias]
"""

# Already in local contact coordinats
fbias = Symbol('fbias', real=True)
fbs1 = Matrix((fbias, 0, fbias))
fbs2 = Matrix((fbias, 0, fbias))


f1 = df1 + fbs1
f2 = df2 + fbs2

pprint(f1)

ksl_v, ksn_v = 500, 100
dbx_v, dbz_v = -2e-3, -2e-3
fbias_v = 1.1

print("\n----- Q1.1 -----")
print("\n--- f1: ")
pprint(f1.subs([[ksn, ksn_v], [ksl, ksl_v], [
       dbx, dbx_v], [dbz, dbz_v], [fbias, fbias_v]]))
pprint(f1.subs([[ksn, ksn_v], [ksl, ksl_v], [
       dbx, dbx_v], [dbz, dbz_v], [fbias, fbias_v]]))
print("\n--- f2: ")
pprint(f2.subs([[ksn, ksn_v], [ksl, ksl_v], [
       dbx, dbx_v], [dbz, dbz_v], [fbias, fbias_v]]))
