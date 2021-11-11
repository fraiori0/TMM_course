# Week 5

Francesco Iori

### Q1

```python
# We leave the derivation of the [T] and [M] matrices as an exercise
# for the reader :-) You will similarly use eq (8) and (9)
from sympy import diag
xuu = simplify(diff(xu, u))
xuv = simplify(diff(xu, v))

Tb = Matrix([list(xuu/normpu), list(xuv/normpvpos)]).transpose()

T = yu.transpose() * Tb
T = simplify(T)
print("\nT = ")
pprint(T)

M = diag(normpu, normpvpos)
print("\nM = ")
pprint(M)
```

which results in

```console
T = 
⎡   -tan(u) ⎤
⎢0  ────────⎥
⎣      R    ⎦

M = 
⎡R     0    ⎤
⎢           ⎥
⎣0  R⋅cos(u)⎦
```

same as in Montana’s paper.

### Q2

