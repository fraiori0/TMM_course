# Week 2

**Francesco Iori**

#### Q1

$$ w_1 = \left[0, -1, 0, 0,0,1\right]  $$



#### Q2


$$
W_n = \left[\begin{array}{ccc}

0 & 0 & 0 \\
-1 & -1 & 1 \\
0 & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0 \\
-1 & 1 & 0

\end{array}
\right]
$$

#### Q3

- rank$(W_n) = 2$
- The unconstrained degrees of freedom (DOFs) are 
  - translation on the x and z axis
  - rotation around the x and y axis. 
- It does not satisfy *force closure*, as the object still has unconstrained DOFs ($W_n$ is not full rank)

#### Q4

$$
W_f
\begin{bmatrix}
  1  &  0  &  0  &  1  &  0  &  0  &  1  &  0  &  0 \\
  0  &  0  & -1  &  0  &  0  & -1  &  0  &  0  &  1 \\
  0  & -1  &  -0  &  0  & -1  &  -0  &  0  &  1  &  -0 \\
  0  & -1  &  -0  &  0  & -1  &  -0  &  0  & -1  &  0 \\
  0  &  1  &  0  &  0  & -1  &  -0  &  0  &  0  &  0 \\
 -1  &  0  & -1  & -1  &  -0  &  1  &  1  &  0  &  0 
\end{bmatrix}
$$

#### Q5

- rank$(W_f) = 6$ 	(full rank, null-space is null)
- This matrix is a good candidate to satisfy *force closure*.

#### Q6


$$
G = \left[
\begin{array}{ccc|ccc|ccc}
  1  &  0  &  0  &  1  &  0  &  0  &  1  &  0  &  0 \\
  0  &  0  & -1  &  0  &  0  & -1  &  0  &  0  &  1 \\
  0  & -1  &  0  &  0  & -1  &  0  &  0  &  1  &  0 \\
  \hline
  0  & -1  &  0  &  0  & -1  &  0  &  0  & -1  &  0 \\
  0  &  1  &  0  &  0  & -1  &  0  &  0  &  0  &  0 \\
 -1  &  0  & -1  & -1  &  0  &  1  &  1  &  0  &  0 \\
 \hline
 -2  &  0  &  0  &  2  &  0  &  0  &  0  &  0  &  0 \\
 -1  & -2  &  0  &  0  &  0  &  0  &  1  &  2  &  0 \\
  0  &  0  &  0  & -1  &  2  &  0  &  1  & -2  &  0 
\end{array}
\right] \;\;\;
\text{for }
k = \begin{pmatrix}
	f_{x1}\\
	f_{y1}\\
	f_{z1}\\
	f_{x2}\\
	f_{y2}\\
	f_{z2}\\
	f_{x3}\\
	f_{y3}\\
	f_{z3}\\
\end{pmatrix}
$$

- Simple example of a force that does not satisfy frictional constraint: an external force acting only on the x-axis of the object.

#### Q7

- Solving the problem for $f_{zi} \geq 1.0$ and for null external forces brings as a result:
  $$
  k_0 = (0,0,1,0,0,1,0,0,2)^T
  $$

- Repeating, considering an external force of $mg = (0,0,-5)^T$, and considering that the grasp should generate an opposite force, we get:
  $$
  k_{mg} = (0, -1.3, 1,0,-1.3,1,0,2.5,2)^T
  $$

  - Frictional forces $(f_{xi}, f_{yi})$ appear, as they are needed to hold the object against an external z-axis force).
  - However, the frictional constraints are not respected (considering $\mu < 1$).

#### Q8

$$
G \cdot k_{mg} = (0,0,5,0,0,0 |0,7.5,-7.5)^T
$$



- Internal forces act on the object. It is reasonable, as the object is squeezed by the $f_{zi}$.