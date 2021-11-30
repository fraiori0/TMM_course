#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from sympy import sin, cos, pi, Matrix, Symbol, symbols, simplify, pprint, lambdify, latex
import numpy as np
import scipy.spatial as sp_spatial
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import os

save = True
fig_name = 'Q23_intersection'

# external wrench
fx, fy, fz, mx, my, mz = symbols('fx,fy,fz,mx,my,mz', real=True)
wrench = Matrix([fx, fy, fz, mx, my, mz])

# transformation between two coordinate frames
rx, ry, rz = symbols('rx, ry, rz', real=True)
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

# utility function, call *.subs with the given value


def Jbitran(rx_v, ry_v, rz_v, thetax_v, thetay_v, thetaz_v):
    return Jbtran.subs([
        [rx, rx_v],
        [ry, ry_v],
        [rz, rz_v],
        [thetax, thetax_v],
        [thetay, thetay_v],
        [thetaz, thetaz_v],
    ])


# Q2.1

Jb1t = Jbitran(0.1, 0.2, 0, 0, 0, 0)
Jb2t = Jbitran(0.1, -0.2, 0, 0, 0, 0)
Jb3t = Jbitran(-0.1, -0.2, 0, 0, 0, 0)
Jb4t = Jbitran(-0.1, 0.2, 0, 0, 0, 0)

w_ext = wrench.subs([
    [mx, 0],
    [my, 0],
    [mz, 0],
])

w1 = Jb1t * w_ext
w2 = Jb2t * w_ext
w3 = Jb3t * w_ext
w4 = Jb4t * w_ext

print("\n---- Q2.1")
print("\n w1:")
pprint(w1)
# pprint(latex(w1))
print("\n w2:")
pprint(w2)
# pprint(latex(w2))
print("\n w3:")
pprint(w3)
# pprint(latex(w3))
print("\n w4:")
pprint(w4)
# pprint(latex(w4))


# Q2.2
# vertices of the limit surface of a single two-tile unit
p_surface = np.array([
    [40, 0, 0],
    [0, 20, 0],
    [-40, 0, 0],
    [0, -20, 0],
    [0, 0, -20]
])

# given that the load is equally distributed on each unit, the limit surface is
# a scaled version of the one for a single unit
# (also because they are all oriented in the same way)
n_units = 4
points = p_surface * n_units
# here we can follow the same approach as in Vector-in-ConvexHull.py
hull = sp_spatial.ConvexHull(points)
indices = hull.simplices
faces = points[indices]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')


for f in faces:
    face = a3.art3d.Poly3DCollection([f])
    face.set_edgecolor('k')
    face.set_alpha(0.3)
    ax.add_collection3d(face)

# For those points in hull.simplices...
for s in hull.simplices:
    s = np.append(s, s[0])  # Cycle back to 1st coordinate
    ax.scatter(points[s, 0], points[s, 1],
               points[s, 2], marker='.', color='b')

# Q2.3
extf = np.array((7, 7, -5))
dhull = sp_spatial.Delaunay(points)
growvector = np.zeros(3)
steps = 100
i = 0
while (dhull.find_simplex(growvector) >= 0):
    print(i)
    growvector = (i/steps)*extf
    i += 1
# print('vector hull intersection:', growvector)

print("\n----- Q2.3")
print("limit for the external force")
pprint(growvector)

ax.scatter(growvector[0], growvector[1], growvector[2],
           marker='o', color='r', s=25)

# set the view and save
ax.view_init(elev=5, azim=5)
plt.show()

if save:
    save_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'figs',
        fig_name+'.png'
    )
    fig.savefig(save_path, dpi=600)
