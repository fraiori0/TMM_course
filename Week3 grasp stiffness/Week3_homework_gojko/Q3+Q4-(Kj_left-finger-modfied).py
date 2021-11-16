#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 21:00:29 2020; Minor edits 24Oct2021 - cutkosky
Compute the "geometric" effect of a compliant grasp that can cause 
instability when multiplied by a substantial grasp force.
The approach loosely follows Kao & Cutkosky (1989) and the Cutkosky
book/thesis.
Basic idea:
    * Make a small motion of the body.
    * Get resulting small motions of the body_contact frames. 
      This is just rigid body motion:  bpJ*dxbody = dxcontact
    * Get small resulting translations on the finger side (no rotations).
    * Get total finger movements using the finger jacobians, Jq.
    * Difference = how much the fingertips moved w.r.t their initial position 
      and orientation (which initially matched body_contact frames). 
    * Build a new 6x6 differential matrix, deltaJtran, that corresponds 
      to this small additional rotation and (for rolling contacts) translation. 
      So, whereas the mapping from grasp force to body force was previously 
       bpJtran*fgrasp = fbody, now we have: pbJtran*deltaJtran*fgrasp

The code below is for the left finger only. You will add the right finger
and then sum for the total effect at the end. Interestingly, from
symmetry, some stuff cancels so the final result is simpler after 
adding the second finger. Final result should match 
Kj in eq (30) of Cutkosky&Kao, for Example 2.
Select and evaluate a few lines at a time, checking intermediate results as
you go to make sure they make sense.
"""

from sympy import sin, cos, pi, Matrix, Symbol, symbols, simplify, pprint

#################################
# Define a 6 element wrench and 6 element body twist
fx, fy, fz, mx, my, mz = symbols('fx,fy,fz,mx,my,mz', real=True)
dbx, dby, dbz = symbols('dbx,dby,dbz', real=True)
qbx, qby, qbz = symbols('qbx,qby,qbz', real=True)
wrench = Matrix([fx, fy, fz, mx, my, mz])
bodytwist = Matrix([dbx, dby, dbz, qbx, qby, qbz])


# We start from the body and work out toward the joints
########################################################
# Jbtran is a 6x6 matrix that, when multipied by a wrench
# in a body contact frame, gives the corresponding body wrench.
# The contact is assumed to be translated by rx,ry,rz from the body
# frame and rotated by thetax,thetay,thetaz (in order).
# Transpose of Jbtran maps a twist from the body frame to the
# corresponding contact frame.
# The general symbolic form is messy (although it looks simple
# in block form as in Cutkosky thesis Appendix A). However,
# often several of the terms [rx,ry,rz,thetax,thetay,thetaz] will be zero.
rx, ry, rz = symbols('rx, ry, rz', real=True, positive=True)
thetax, thetay, thetaz = symbols('thetax,thetay,thetaz', real=True)

Rotx = Matrix([[1, 0, 0], [0, cos(thetax), -sin(thetax)],
              [0, sin(thetax), cos(thetax)]])
Roty = Matrix([[cos(thetay), 0, sin(thetay)], [
              0, 1, 0], [-sin(thetay), 0, cos(thetay)]])
Rotz = Matrix([[cos(thetaz), -sin(thetaz), 0],
              [sin(thetaz), cos(thetaz), 0], [0, 0, 1]])
Amat = Rotx*Roty*Rotz

Rskew = Matrix([[0, -rz, ry], [rz, 0, -rx], [-ry, rx, 0]])

top = Amat.row_join(Matrix.zeros(3))
bottom = (Rskew*Amat).row_join(Amat)
Jbtran = top.col_join(bottom)
#########################################################
# Do some quick tests to check that it looks OK
# Jtest1 = Jbtran.subs([(thetax,0),(thetay,0),(thetaz,0)]) #only translated
# Jtest2 = Jbtran.subs([(thetax,0),(thetay,0),(rz,0)])     #planar offset


#########################################################
# Create the H matrix that filters out rotational terms at contacts
# H*twist = twist_transmitted  (linear motion only unless using
# "soft" fingers). See Appendix I of Cutkosky&Kao
# Htran*wrench_transmitted = contact wrench  (nonzero forces only)
Hmat = (Matrix.eye(3)).row_join(Matrix.zeros(3))  # 3x6


#########################################################
# Jacobian for left, right fingers per Appendix II B, Example 2
# of Cutkosky&Kao 1989
link = Symbol('link', real=True, positive=True)  # link length
Jq1 = Matrix([[0, -link, -link], [link, 0, 0], [0, link, 0],
             [0, 0, 0], [0, 1, 1], [-1, 0, 0]])
Jq2 = Matrix([[0, -link, -link], [link, 0, 0],
             [0, -link, 0], [0, 0, 0], [0, 1, 1], [1, 0, 0]])

# Helps to have a little coord. frame toy when building these

# Check if looks OK, using the left finger:
dq1, dq2, dq3 = symbols('dq1,dq2,dq3', real=True)  # small joint motions
djoints = Matrix([dq1, dq2, dq3])
dfcontact = Jq1*djoints    # total twist motion of the fingertip
dtrans = Hmat*Jq1*djoints  # 3 element vector of velocities in (l,m,n) frame


#########################################################
# Now specialize bpJtran for our object and contact frames.
# Helps to have a little coord. frame toy when building these:
# From (x,y,z) to (l1,m1,n1) we translate -w and rotate -pi/2
# about y, and then -pi/2 about z.
w = Symbol('w', real=True, positive=True)
Jb1t = Jbtran.subs([(rx, -w), (ry, 0), (rz, 0), (thetax, 0),
                   (thetay, -pi/2), (thetaz, -pi/2)])
Jb2t = Jbtran.subs([(rx, w), (ry, 0), (rz, 0), (thetax, 0),
                   (thetay, pi/2), (thetaz, pi/2)])

"""
Up to this point everything has been the same as for getting [Kb] the "direct"
stiffness term. Now we start to look at the effects of differential motions.
"""

"""
Using method in Cutkosky thesis/book/report Appendix A we can
define (Jtran_new - Jtran) = deltaJtran where Jtran_new is after 
a small additional rotation and translation. Here we take the 
version from Cutkosky&Kao Appendix I, eq (42) to work in 
local instead of world frame.  Note that there's a typo in the paper
just before eq (42): what are called deltaA' and deltaR' are actually
not the transposes but are the blocks that go directly into (42).

We build deltaJtran by blocks (upper left, upper right, etc):
"""

deltax, deltay, deltaz = symbols('deltax,deltay,deltaz', real=True)
dqx, dqy, dqz = symbols('dqx,dqy,dqz', real=True)  # small rotations
# Note that deltax, dqx etc are in the local contact coordinate frame
deltaA = Matrix([[0, -dqz, dqy], [dqz, 0, -dqx], [-dqy, dqx, 0]])
deltaR = Matrix(
    [[0, -deltaz, deltay], [deltaz, 0, -deltax], [-deltay, deltax, 0]])
dJtop = deltaA.row_join(Matrix.zeros(3))
dJbot = deltaR.row_join(deltaA)
deltaJtran = dJtop.col_join(dJbot)
deltaJtran.shape  # should be 6x6

"""
Make a small motion of the body and find the corresponding
motions of the fingers
"""
bctwist1 = Jb1t.T*bodytwist  # contact frames on the body
bctwist2 = Jb2t.T*bodytwist

ctwist1 = Hmat*bctwist1  # what gets transmitted through the contact
ctwist2 = Hmat*bctwist2

Jq1inv = (Hmat*Jq1)**-1
Jq2inv = (Hmat*Jq2)**-1

djoints1 = Jq1inv * ctwist1  # corresponding joint motions
djoints2 = Jq2inv * ctwist2

ftwist1 = Jq1*djoints1  # Work back outward to get the total motion
ftwist2 = Jq2*djoints2

# Simplify a bit by letting 'link' = 1 as in example
ftwist1 = ftwist1.subs(link, 1)
ftwist2 = ftwist2.subs(link, 1)

# Relative motion of fingertip to object in (l,m,n) frame
diff1 = ftwist1-bctwist1
diff2 = ftwist2-bctwist2

# Put the elements of diff1 and diff2 into our differential transforms:
dJtran1p = deltaJtran.subs([(dqx, diff1[3]), (dqy, diff1[4]),
                           (dqz, diff1[5]), (deltax, 0), (deltay, 0), (deltaz, 0)])
dJtran2p = deltaJtran.subs([(dqx, diff2[3]), (dqy, diff2[4]),
                           (dqz, diff2[5]), (deltax, 0), (deltay, 0), (deltaz, 0)])

"""Modifications for Q4"""
r = Symbol('r', real=True, positive=True)
dJtran1pQ4 = deltaJtran.subs([(dqx, diff1[3]), (dqy, diff1[4]), (
    dqz, diff1[5]), (deltax, r*diff1[4]), (deltay, -r*diff1[3]), (deltaz, 0)])
dJtran2pQ4 = deltaJtran.subs([(dqx, diff2[3]), (dqy, diff2[4]), (
    dqz, diff2[5]), (deltax, r*diff2[4]), (deltay, -r*diff2[3]), (deltaz, 0)])
# Note: This is where we would change things for rolling (e.g., 'deltax' not zero)
# Note also that there could be rolling in two orthogonal directions.

# Now impose our grasp forces (in contact frame)
fn = Symbol('fn', real=True, positive=True)
fgrasp = Matrix([0, 0, -fn, 0, 0, 0])

dfbody1 = -Jb1t*dJtran1p*fgrasp
dfbody2 = -Jb2t*dJtran2p*fgrasp

dfbody1Q4 = -Jb1t*dJtran1pQ4*fgrasp
dfbody2Q4 = -Jb2t*dJtran2pQ4*fgrasp

#print("Change in force on body due to small geometry change [dfbody]")
# pprint(dfbody1)

"""
You will add the corresponding terms for the right finger above and then
get the total change in the body force. Then you can use sympy to take
the jacobian of dfbody with respect to bodytwist to get [Kj]
"""
dfbody = simplify(dfbody1+dfbody2)
Kj = dfbody.jacobian(bodytwist)  # Matches Kj in eq (30) Cutkosky&Kao
print("Geometric stifness:")
pprint(Kj)

dfbodyQ4 = simplify(dfbody1Q4+dfbody2Q4)
KjQ4 = dfbodyQ4.jacobian(bodytwist)
print("Geometric stifness for hemispherical fingertip:")
pprint(KjQ4)

print("Critical value r=w, Kj >= 0")
