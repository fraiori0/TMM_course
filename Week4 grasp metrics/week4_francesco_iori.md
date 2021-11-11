# Week 4

Francesco Iori

### Q2

###### Q 2.1

To include the origin in the convex hull the fingers must be able to exert a (total) force in all directions, requiring 
$$
\mu > 1.0
$$
Physically, this means that the fingers on the sides of the trapezoid can also exert a force with a negative $y$-axis component.

###### Q 2.2

Let’s consider that the resulting force on the object, $f_b$:

- does not depend on the position of the contact (upper part of the grasp matrix);
- due to friction, it depends on how the contact surfaces are oriented;

and that the resulting torque on the object, $\tau_b$, is linearly affected by the position of the contacts AND the orientation of the contact surfaces, due to friction.

To maximize the convex hull, we can consider to be in the limit case where each finger is exerting the maximum amount of force (max $f_{ix}$ and $f_{iy}$, i.e.  limit of the friction cone).

To maximize the resulting $f_b$ on the object we would like our contact to be able to exert a force equal in all direction (e.g. if the shape was a cube, placing one finger on each side). 

Placing one finger on top and one on the bottom solve this issue for the $y$ component. For the $x$ component we would like to have one finger on both the lateral sides of the trapezoid, even if they are not oriented normally to the $x$-axis.

For $\tau_b$ we have an heuristic “opposite” to the one for resisting forces on the $y$-axis, i.e. placing the fingers as far as possible from the origin, and oriented such at the limit (max $f_{ix}$ and $f_{iy}$, limit of the friction cone). This would mean placing two fingers on the upper corner of the object, far away from the origin.

After some experimenting, the highest value I found was with

```python
frames[0, :] = [-3, 1, 3*np.pi/4]
frames[1, :] = [3, 1, np.pi/4]
frames[2, :] = [-3, 1, -3*np.pi/4]
frames[3, :] = [3, 1, -np.pi/4]
```

with a grasp quality given by the Convex Hull of the Union **= 0.6729**

This is found assuming that, on the corner of the trapezoid, the finger can be placed with any orientation that is in the range given by the normals to the two surfaces forming that corner.

![optimal_grasp](./optimal_grasp.png)

###### Q 2.3

For the same placement, the quality of the grasp given by the Convex Hull of the Minkowski sum is **= 1.41**

This value is higher because, to generate the convex hull, it considers all the possible combination of the contact forces. On the other hand, the union operation is very conservative, as (at most, considering that all points generate a vertex of the convex hull) it considers only a weighted mean of the contact forces (instead of their sum).

