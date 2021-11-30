'''
Compute the convex hull of a bunch of [x,y,z] points and test whether
a given point is inside the hull.
Then plot the hull and vector from the origin to the same point.
And plot the intersection point to  face, if there is one.
21Nov2021 Cutkosky
'''

import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy import spatial as sp_spatial
import numpy as np
import os

"""
Define the vertices (x,y,z) of a poyhedron.
"""
points = np.array([[0, 0, 0], [0, -1, 0], [0, 1, 0], [-1, 1, 2],
                   [-1, -1, 2], [-1, 0, 2], [3, -1, 2], [3, 1, 2],
                   [3, -1, -0.75], [3, 1, -0.75], [3, 0, -0.25], [3, 0, 2]])


# Get the convex hull of the points
hull = sp_spatial.ConvexHull(points)
indices = hull.simplices
faces = points[indices]

"""
Define our vector here (change this to suit)
"""
# extf = np.array([1, 0, 3])
extf = np.array([0.1, 0, 1.3])

save = True
fig_name = 'Q12_polyhedron'

# https://www.py4u.net/discuss/12821
dhull = Delaunay(points)
print('extf inside convex hull?', dhull.find_simplex(extf) >= 0)

"""
If the vector falls outside the convex hull, we can find the 
intersection point (brute force method).
"""
growvector = np.zeros(3)
if((dhull.find_simplex(extf) >= 0) == False):
    i = 0
    steps = 100
    while((dhull.find_simplex(growvector) >= 0) == True):
        growvector = (i/steps)*extf
        i += 1
    print('vector hull intersection:', growvector)

"""
Plot the convex hull and a vector from [0,0,0] to extf and the intersection
point (if there is one).
"""

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
    ax.scatter(points[s, 0], points[s, 1], points[s, 2], marker='.', color='b')

ax.scatter(extf[0], extf[1], extf[2], marker='o', color='r')
ax.plot([0, extf[0]], [0, extf[1]], [0, extf[2]], color='r')
ax.scatter(growvector[0], growvector[1], growvector[2], marker='o', color='k')

ax.view_init(elev=5, azim=-90)

plt.show()

if save:
    save_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'figs',
        fig_name+'.png'
    )
    fig.savefig(save_path, dpi=600)
