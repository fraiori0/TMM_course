# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg.linalg import _solve_dispatcher
from scipy.optimize import root_scalar, newton
from pprint import pprint

np.set_printoptions(precision=3)

filename = "Week8 manipulation with sensing/wrench-sequence.txt"

# Read file to retrieve the generated sensor readings
wrenches = np.loadtxt(filename, delimiter="\t", encoding='bytes')

# Following the formulation illustrated in the class, we can compute the position of the contact (given the sensed wrench)
# Normalized forces
f_v = wrenches[:, :3]/np.linalg.norm((wrenches[:, :3]), axis=1)[:, None]
# Normalized torques
tau_v = wrenches[:, 3:]/np.linalg.norm((wrenches[:, 3:]), axis=1)[:, None]
# Normalized lever arm
h_v = np.cross(f_v, tau_v, axis=1)
# Norm of lever arm, h_norm = f_norm/tau_norm
h_norm = np.linalg.norm((wrenches[:, 3:]), axis=1)[
    :, None] / np.linalg.norm((wrenches[:, :3]), axis=1)[:, None]

# sanity check, they should have all norm = ~1.0
# pprint(np.linalg.norm(h_v, axis=1))

# solve the equation r = h - alpha*f_v s.t. ||r||=R
R = 1.0


# cast the problem as finding the roots of f(alpha) = r(alpha)**2 - R**2
def r_alpha(alpha, i):
    return h_v[i] * h_norm[i] - alpha * f_v[i]


def f_root(alpha, i):
    return r_alpha(alpha, i).dot(r_alpha(alpha, i)) - R**2


# solve root-finding problem for each point
alpha_sol = []
r_sol = []
for i in range(wrenches.shape[0]):
    # start from both side (alpha less or greater than zero)
    sol1 = newton(f_root, args=(i,), x0=-1)
    sol2 = newton(f_root, args=(i,), x0=1)
    # select the right solution using the value of the scalar product
    r1 = r_alpha(sol1, i)
    if r1.dot(f_v[i]) <= 0:
        alpha_sol.append(sol1)
        r_sol.append(r1)
    else:
        r2 = r_alpha(sol2, i)
        alpha_sol.append(sol2)
        r_sol.append(r2)

r_sol = np.round(np.array(r_sol), decimals=4)
# same result as from the generator, in the noiseless case
pprint(r_sol)
