#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 19:54:23 2019; minor edits 15Oct2021.
Utilties for wrenches and twists. Much of this surely exists already in some Scipy or other libraries
but it sometimes is less convenient to use. Anyway, we use the following extensively for constructing
grasp wrench matrices and moving twists and wrenches from one frame to another.
@author: cutkosky
"""

import numpy as np

#3x3 orthonormal rotation matrices
def Rotx(thetax):
    st = np.sin(thetax)
    ct = np.cos(thetax)
    Rot = np.array([[1,0,0],[0,ct,-st],[0,st,ct]])
    return Rot
def Roty(thetay):
    st = np.sin(thetay)
    ct = np.cos(thetay)
    Rot = np.array([[ct,0,st],[0,1,0],[-st,0,ct]])
    return Rot
def Rotz(thetaz):
    st = np.sin(thetaz)
    ct = np.cos(thetaz)
    Rot = np.array([[ct,-st,0],[st,ct,0],[0,0,1]])
    return Rot

#3x3 skew symmetric cross product matrix: Rcross*vector = R x vector
def Rcross(rx,ry,rz):
    Skew = np.array([[0,-rz,ry],[rz,0,-rx],[-ry,rx,0]])
    return Skew

"""
6x6 cartesian matrix Jb that maps a velocity or incremental
motion (vx,vy,vz,omegax,omegay,omegaz) from one coordinate
frame to another, which is translated by [rx,ry,rz]
and (then) rotated by RPY angles [thetax,thetay,thetaz]
with respect to the first.
Transpose of Jb will map a wrench from second frame to the first.
Notation matches Appendix A of Cutkosky thesis.
"""
def Cartesmap(rx,ry,rz,thetax,thetay,thetaz):
   Rx = Rotx(thetax);   Ry = Roty(thetay);   Rz = Rotz(thetaz)
   Amat = np.linalg.multi_dot((Rx,Ry,Rz))
   Rskew = Rcross(rx,ry,rz)
   UL = np.transpose(Amat)
   UR = np.dot(UL,np.transpose(Rskew))
   LL = np.zeros((3,3))
   LR = UL
   Jb = np.block([[UL,UR],[LL,LR]])
   return Jb

"""
Planar version for twist (vx, vy, omega) with two frames,
the second being located [rx,ry,thetaz] with respect to first
"""
def PTrans(x,y,theta):
    ct = np.cos(theta)
    st = np.sin(theta)
    Jbp = np.array([[ct, st, (x*st-y*ct)],[-st, ct, (x*ct+y*st)],[0,0,1]])
    return Jbp
