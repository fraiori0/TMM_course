# -*- coding: utf-8 -*-
"""
Examples of using linprog(), as well as checking the rank, to see if
a wrench matrix of contacts satisfies Force Closure while keeping 
normal forces positive inward.
13Jan2020; minor edits 15Oct2021 -mrc
"""

import numpy as np
from scipy.linalg import null_space
from numpy.linalg import matrix_rank
from scipy.optimize import linprog


#First try numbers to match Example 12.7 in Lynch & Park
F = np.array([[0,0,-1,2],[-1,0,1,0],[0,-1,0,1]])  #planar wrench matrix
Frank = matrix_rank(F)  #Needs to be 3 in the plane
nullF = null_space(F)

#scipy.optimize.linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, 
#bounds=None, method='simplex', callback=None, options=None)
f= np.array([1,1,1,1])   #weights
b = -f  #This is for the inequality constraints (unisense forces)
A = -np.identity(4)
Aeq = F  #This is for the homogeneous solution
beq = np.array([0,0,0])
#Minimize f'*k 
#Want F*k = 0 with all k positive (unisense) for internal forces.
#f = weights, A*k <= b gives the unisense constraint 
# b=-1 is convenient choice and makes k >= 1
#F is the homogeneous grasp solution: Aeq=F so Aeq*k = beq = 0
result=linprog(f,A,b,Aeq,beq)

"""
Do a planar example for 3 fingers holding a planar block.
Fingers are at locations [0,-1], [-1,1] and [1,1] in the XY plane.
Assume y1 is upward, y2 & y3 point downward. 
Assume force fx1 to the right, fx2,fx3 to left
"""

#Unlike Lynch&Park I prefer my planar wrenches to be 
# [fx,fy,mz]' 
#So for frictionless fingers we have 3 columns in [Gtn]
#These are normal to the block and sign must be preserved (unisense)
Gtn = np.array([[0,0,0],[1,-1,-1],[0,1,-1]])
matrix_rank(Gtn)  #rank = 2 (inadequate)
null_space(Gtn)   #indeed, we have a problem...

#Now add the tangential wrenches due to friction
#These can be positive or negative sense
Gtf = np.array([[1,-1,-1],[0,0,0],[1,1,1]])
Gt = np.column_stack((Gtn,Gtf))
matrix_rank(Gt)   #rank = 3 (good)

#Full rank is necessary but not sufficient for Force Closure.
#See if solution exists with first 3 wrenches positive
#Follow the approach of example from Lynch & Park
w= np.array([1,1,1,1,1,1])  #weights
Aeq = Gt
beq = np.zeros(3)
#inequality constraint on the 1st three elements
b=np.array([-1,-1,-1,0,0,0])

Atop = np.column_stack((-np.identity(3),np.zeros((3, 3))))
Abottom = np.zeros((3,6))
A = np.concatenate((Atop,Abottom))

result1=linprog(w,A,b,Aeq,beq)
#Check w to see if a solution exists. 
#If it does we have Force Closure
np.set_printoptions(precision=3,suppress=1)
print(result1.x)

########
#An alternative, which might be more intuitive, is to set bounds directly
f1ybound = (1,None)
f2ybound = (1,None)
f3ybound = (1,None)
f1xbound = (None,None)  #Could make it bounded by friction if we want...
f2xbound = (None,None)
f3xbound = (None,None)
bounds = (f1ybound, f2ybound, f3ybound, f1xbound, f2xbound, f3xbound)

#find x that minimizes w'*x  subject to Aeq*x=beq and bounds
result2 = linprog(w,None,None,Aeq,beq,bounds)
print(result2.x)
