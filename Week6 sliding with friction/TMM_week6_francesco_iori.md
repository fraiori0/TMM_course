# Week 6

Francesco Iori

### Q1

###### Q1.1

There are 2 unknowns for each contact point (the value of the components of the friction force in that point, $f_{ix}$ and $f_{iy}$). In this case, the total is 6.

###### Q1.2

There are 2 equalities equations and 3 inequality equations.

In the general case, there are:

- 3 equality equations, for the 2D equilibrium of the body;
- $n$ inequality equations, due to frictional constraint ($\sqrt{f_{ix}^2+f_{iy}^2}\leq\mu f_{in}$).

In the special case being considered in Fig. 4.2.7 of Sakurai's Thesis the equality equations reduce to 2, owing to the particular choice of $f_c$ and the position of the contact points.

###### Q1.3

```python
# We approximate the friction limits sqrt(fx^2+fy^2) <= 1 with
# a crude friction square. We fudge the limits a bit so the
# square fits the bounding circle a bit better on average.
# Sakurai uses 36-sided polygons.
mu = 1.0

bounds = [[-mu*0.9, mu*0.9] for i in range(6)]

sol = linprog(c=f, A_eq=Aeq, b_eq=beq, bounds=bounds, method='simplex')
print(sol)
```

which outputs

```console
  sol = linprog(c=f, A_eq=Aeq, b_eq=beq, bounds=bounds, method='simplex')
     con: array([-0.,  0., -0.])
     fun: -3.3941125496954285
 message: 'Optimization terminated successfully.'
     nit: 13
   slack: array([], dtype=float64)
  status: 0
 success: True
       x: array([-0.9, -0.6, -0. , -0.9,  0.9, -0.9])
```

which is consistently different from the result in Sakurai's Thesis.

The difference may be due to the 4-sided polygon approximation of the frictional constraints.

It is also interesting to note that choosing a different optimization method (e.g. `interior-point` ) gives a different solution.

###### Q1.4

The external force does not produce any torque on the body; furthemore, the problem becomes symmetric on the $y$-axis. This means that we should not expect any rotation or motion on the $x$-axis.

Running the script again produces

```console
sol = linprog(c=f, A_eq=Aeq, b_eq=beq, bounds=bounds, method='simplex')
     con: array([ 0., -0.,  0.])
     fun: -2.7
 message: 'Optimization terminated successfully.'
     nit: 12
   slack: array([], dtype=float64)
  status: 0
 success: True
       x: array([-0.9, -0.9,  0.9, -0.9,  0. , -0.9])
```

The solution is intuitevely correct for the $f_{iy}$ component, but not necessarily for $f_{1x}$ and $f_{2x}$. In fact, the equilibrium equation on the $x$-axis just requires that these two forces are equal with respect to each other. Using a better approximation for the friction constraint should give $f_{1x}, f_{2x}\rightarrow$ 0 . 

The magnitude of the external force needed to initiate the sliding in this case is $f_c = 3\,\mu\,f_n$ .

In this case the answer should be more accurate because the direction of the friction force coincide with one of the side of the pyramid used to approximate the constraint.

### Q2

Here we are considering a rotation around the first contact point.

We get linear equations as we are consideing the COR to be known (and so we know the total wrench due to friction on the other contact points). As the friction constraint is not violated, the solution to the linear equations gives an exact solution.

This solution can be considered more accurate then the one given by the ellipsoidal approximation. Indeed, this case corresponds exactly to a faucet of the (real) "ellipsoid", where the approximation is less precise.

### Q3

In this case, if we assume a rotation around $P_4$ the friction limit is exceed, and the assumption thus wrong. Furthermore, even if we consider the COR to be in $P_3$ we find that such assumption is still wrong.

This means that the COR is not one of the contact point, and the general solution with the ellipsoid gives us the solution.