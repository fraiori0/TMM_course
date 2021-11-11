#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 17:51:21 2019. @author: cutkosky
Sympy adaptation of earlier Mathematica script to follow procedure
in Montana 'Kinematics of Contact and Grasp', IJRR 1988. We follow
Equations (6-15) to derive the [K] matrix for the case of a sphere.
3Nov2021 minor edits
"""
from sympy import sin, cos, diff, Symbol, Matrix, diag
from sympy import simplify, refine, Q, pprint

R = Symbol('R', positive=True, real=True)
u = Symbol('u', real=True)
v = Symbol('v', real=True)

# Use f() from equation 11 in Montana for a sphere
# Let (u,v) be latitude and longitude, with u pointing north.
fvec = Matrix([R*cos(u)*cos(v), -R*cos(u)*sin(v), R*sin(u)])

# Create x(u), y(u), z(u) per equation 6.
partialu = diff(fvec, u)
partialv = diff(fvec, v)

# Get norms ||f_u(u)|| etc.
normpu = simplify(partialu.norm())
normpv = simplify(partialv.norm())
# Help sympy to recognize some algebraic simplifications
# Get rid of abs(cos(u)) stuff knowing that u is from -pi/2 to pi/2
normpvpos = refine(normpv, Q.positive(cos(u)))

# xu, yu, zu should match equation 14
xu = partialu/normpu
yu = partialv/normpvpos
zu = simplify(xu.cross(yu))

pprint(zu)

# Assemble K per equation 7
Ka0 = xu.transpose()
Ka1 = yu.transpose()
Ka = Ka0.col_join(Ka1)

Kb0 = simplify(diff(zu, u))/normpu
Kb1 = diff(zu, v)/normpvpos
Kb = Kb0.row_join(Kb1)

# Kmat should match the [K] in equation 15
Kmat = simplify(Ka*Kb)

pprint(Kmat)

# We leave the derivation of the [T] and [M] matrices as an exercise
# for the reader :-) You will similarly use eq (8) and (9)
xuu = simplify(diff(xu, u))
xuv = simplify(diff(xu, v))

Tb = Matrix([(xuu/normpu).transpose(),
            (xuv/normpvpos).transpose()]).transpose()

T = yu.transpose() * Tb
T = simplify(T)
print("\nT = ")
pprint(T)

M = diag(normpu, normpvpos)
print("\nM = ")
pprint(M)

"""
To actually do rolling you will need to have 2 sets of these matrices
for the two coordinate frames, and keep updating
in a loop using equations 16-20
"""
