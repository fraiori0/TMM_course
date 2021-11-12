import numpy as np
from scipy.linalg.decomp_svd import null_space
from WrenchUtils import *
from math import pi
import array_to_latex as a2l
from numpy.linalg import matrix_rank
import scipy.linalg as linalg
import sys
from scipy.optimize import linprog


# debugging
def deb():
    sys.exit('Yeeet')


np.set_printoptions(precision=1)

# Contact location and local reference frame rotation
tp1 = np.array([1, 1, 0])
th1 = np.array([-pi/2, 0, 0])

tp2 = np.array([-1, 1, 0])
th2 = np.array([-pi/2, 0, 0])

tp3 = np.array([0, -1, 0])
th3 = np.array([pi/2, 0, 0])

tp = np.array([tp1, tp2, tp3])
th = np.array([th1, th2, th3])

# Contact forces in local reference frame
wc1 = np.array([0, 0, -1, 0, 0, 0])
wc2 = np.array([0, 0, -1, 0, 0, 0])
wc3 = np.array([0, 0, -1, 0, 0, 0])
wc = np.array([wc1, wc2, wc3])

# Wn matrix
Wn = np.empty((6, wc.shape[0]))
for i, (tpi, thi, wci) in enumerate(zip(tp, th, wc)):
    Wn[:, i] = Cartesmap(*tpi, *thi).T.dot(wci)

print('Q2-------\nWn:')
print(np.round(Wn))
print('\n')

print('Q2-------\nWn:')
print('1. Rank:\t', matrix_rank(Wn))
print('\n')


wc1f = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0],
]).T
wc2f = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0],
]).T
wc3f = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0],
]).T
# np.concatenate((wc1f, wc2f, wc3f), axis=1)
wcf = np.array((wc1f, wc2f, wc3f))
print(wcf.shape)
Wnf = []
for i in range(wcf.shape[0]):
    Wnf.append(Cartesmap(*tp[i], *th[i]).T.dot(wcf[i]))
Wnf = np.concatenate(Wnf, axis=1)


print('Q4-------\nWnf:')

print(a2l.to_ltx(np.round(Wnf), frmt='{:1.0f}'))
print('\nMatrix rank w/ friction:')
print(matrix_rank(Wnf))
print('\n')


print('Q6-------\n Internal forces: ')

Wnf_int = np.array((
    np.concatenate((tp2-tp1, -(tp2-tp1), np.zeros(3))),
    np.concatenate((tp3-tp1, np.zeros(3), -(tp3-tp1))),
    np.concatenate((np.zeros(3), -(tp3-tp2), tp3-tp2))
))

print(Wnf_int.shape)
G = np.concatenate((Wnf, Wnf_int), axis=0)
G[np.isclose(G, 0)] = 0
print('\nG:')
print(a2l.to_ltx(np.round(G), frmt='{:1.0f}'))


print('Q7-------\n')
print('linprog:')
f = np.ones(9)  # weights
b = -np.ones(3)  # This is for the inequality constraints (unisense forces)
# b = -np.array((1, 1, 1, 0, 0, 0, 0, 0, 0))
A = -np.array((
    (0, 0, 1, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 1)
))
w_ext = np.array([0, 0, 5, 0, 0, 0])
lb = np.array([-10, -10, 8, -10, -10, 8, -10, -10, 8])
ub = np.array([10, 10, 10, 10, 10, 10, 10, 10, 10])
bounds = np.array((lb, ub))
# Minimize f'*k
# Want F*k = 0 with all k positive (unisense) for internal forces.
# f = weights, A*k <= b gives the unisense constraint
# b=-1 is convenient choice and makes k >= 1
# F is the homogeneous grasp solution: Aeq=F so Aeq*k = beq = 0
# result = linprog(f, A_ub=A, b_ub=b, A_eq=G[:6], b_eq=w_ext)
result = linprog(f, A_ub=None, b_ub=None,
                 A_eq=G[:6], b_eq=w_ext, bounds=bounds.T)
print(result)
sol = result.x
sol[np.isclose(sol, 0)] = 0
print('Q7-------\nk =', sol)

print('Q8-----------\nG.dot(k):')
f_b = G.dot(result.x)
f_b[np.isclose(f_b, 0, atol=1e-4, rtol=1e-4)] = 0
print(f_b)


# print('Q5-------\n')

# print('Q5-------\n')

# print('Q5-------\n')
