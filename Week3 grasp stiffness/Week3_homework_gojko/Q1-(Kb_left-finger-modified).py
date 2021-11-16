
"""
Created on Fri Jan 20 21:00:29 2020 @author: cutkosky
Create the "direct" grasp stiffness that arises from joint stiffness
and passive structural stiffness in the fingers.
The approach roughly follows Cutkosky & Kao (1989) and the Cutkosky
book/thesis.
Basic idea: f = Kx where forces and small motions are mapped using
transformation matrices and jacobians among body, contacts and joints.

This file has the instructions for the left finger. The total stiffness
can be computed by adding in expressions for the right finger and summing.
Recommended use is to evaluate a few lines at a time, checking intermediate
results as you go, to confirm they look correct and don't have sign errors, etc.
"""

from sympy import sin, cos, pi, Matrix, Symbol, symbols, simplify

#################################
# Define 6 element wrench and 6 element body twist
fx,fy,fz,mx,my,mz = symbols('fx,fy,fz,mx,my,mz',real=True)
dbx,dby,dbz = symbols('dbx,dby,dbz',real=True)
qbx,qby,qbz = symbols('qbx,qby,qbz',real=True)
wrench = Matrix([fx,fy,fz,mx,my,mz])
bodytwist = Matrix([dbx,dby,dbz,qbx,qby,qbz])

"""
We start from the body and work out toward the joints
First we create [Jbtran].
Jbtran is a 6x6 matrix that, when multipied by a wrench 
in a body_contact frame, gives the corresponding wrench
at the center of mass. 
Contact is assumed to be translated by rx,ry,rz from the body 
frame and rotated by thetax,thetay,thetaz (in order).
Transpose of Jbtran maps a twist from the body frame to the
corresponding body_contact frame.
The general symbolic form is messy (although it looks simpler
in block form as in Cutkosky thesis Appendix A). 

We will specialize Jbtran for the left and right finger later.
Often several of the terms [rx,ry,rz,thetax,thetay,thetaz] will be zero
in applications and examples.
"""
rx, ry, rz = symbols('rx, ry, rz', real = True, positive = True)
thetax,thetay,thetaz = symbols('thetax,thetay,thetaz', real = True)

Rotx = Matrix([[1,0,0],[0,cos(thetax),-sin(thetax)],[0,sin(thetax),cos(thetax)]])
Roty = Matrix([[cos(thetay),0,sin(thetay)],[0,1,0],[-sin(thetay),0,cos(thetay)]])
Rotz = Matrix([[cos(thetaz),-sin(thetaz),0],[sin(thetaz),cos(thetaz),0],[0,0,1]])
Amat = Rotx*Roty*Rotz                                #3x3 orthonormal orientation    

Rskew = Matrix([[0,-rz,ry],[rz,0,-rx],[-ry,rx,0]])   #cross product matrix

top = Amat.row_join(Matrix.zeros(3))
bottom = (Rskew*Amat).row_join(Amat)
Jbtran = top.col_join(bottom)
#########################################################
#Do some quick tests to check that it looks OK
#Jtest1 = Jbtran.subs([(thetax,0),(thetay,0),(thetaz,0)]) #if only translated
#Jtest2 = Jbtran.subs([(thetax,0),(thetay,0),(rz,0)])     #if planar offset


"""
Create the H matrix that filters out rotational terms at contacts
H*twist = twist_transmitted  (linear motion only unless using
"soft" fingers). See Appendix I p162 of Cutkosky&Kao
Htran*wrench_transmitted = contact wrench  (nonzero forces only)
"""
Hmat = (Matrix.eye(3)).row_join(Matrix.zeros(3))   #3x6


"""
Jacobian for left finger per Appendix III B, Example 2, p163
of Cutkosky&Kao 1989 assuming links of length "link"
"""
link = Symbol('link',real=True,positive=True)    #link length
Jq1 = Matrix([[0,-link,-link],[link,0,0],[0,link,0],[0,0,0],[0,1,1],[-1,0,0]])
Jq2 = Matrix([[0,-link,-link],[link,0,0],[0,-link,0],[0,0,0],[0,1,1],[1,0,0]])

#It helps to have a little coordinate frame toy handy when building these.
#Check if looks OK:
dq1,dq2,dq3 = symbols('dq1,dq2,dq3',real=True)    #small joint motions
djoints = Matrix([dq1,dq2,dq3])
dfcontact = Jq1*djoints    # total twist motion of the fingertip
dtrans = Hmat*Jq1*djoints  # 3 element vector of velocities in (l,m,n) frame
#dtrans is the vector of twist elements that are transmitted through a
#point contact for the left finger. 


"""
Now specialize bpJtran for our object and coordinate frames.
It helps to have a little coord. frame toy when building these:
From (x,y,z) to (l1,m1,n1) we translate -w and rotate -pi/2
about y, and then -pi/2 about z.
"""
w = Symbol('w',real=True,positive=True)
Jb1t = Jbtran.subs([(rx,-w),(ry,0),(rz,0),(thetax,0),(thetay,-pi/2),(thetaz,-pi/2)]) 
Jb2t = Jbtran.subs([(rx,w),(ry,0),(rz,0),(thetax,0),(thetay,pi/2),(thetaz,pi/2)])

"""
Now do the joint space compliance. Assume a diagonal joint stiffness matrix
This is upper left block of Ktheta in Appendix III B, Example 2.
"""
ka,kb,kc,kq = symbols('ka,kb,kc,kq',real=True,positive=True)

Ktheta = Matrix.eye(3)
Ktheta[0,0] = ka; Ktheta[1,1] = kb; Ktheta[2,2] = kc
#We could make it even simpler: ka=kb=kc=kq
#Ktheta = kq*Matrix.eye(3)


"""
Get the Kp (body_contact stiffness matrix) for each finger
"""
Ctheta = Ktheta**-1
Cf1 = Jq1*Ctheta*Jq1.T  # 6x6 compliance matrix (usually singular unless structural compliance)
Cf2 = Jq2*Ctheta*Jq2.T

#Things are getting a bit messy so lets set 'link' = 1 as in example
Cf1 = Cf1.subs(link,1)
Cf2 = Cf2.subs(link,1)
Kfp1 = simplify((Hmat*Cf1*Hmat.T)**-1)  # 3x3 linear stiffness matrix at contact
Kfp2 = simplify((Hmat*Cf2*Hmat.T)**-1)

Kp1 = Hmat.T*Kfp1*Hmat
Kp2 = Hmat.T*Kfp2*Hmat

########################################################
#Map the body_contact stiffness to the body center frame
Kb1 = simplify(Jb1t*Kp1*Jb1t.T)
Kb2 = simplify(Jb2t*Kp2*Jb2t.T)

Kbtotal = Kb1+Kb2
#You can compare this with eq (29) in Cutkosky & Kao
#where 'w' is called 'r' in the paper.

print(Kbtotal)