# Week 3

Francesco Iori

### Q1

Considering a jacobian of the right finger as 
$$
J_{\theta_2}\left[\begin{matrix}0 & - link & - link\\link & 0 & 0\\0 & - link & 0\\0 & 0 & 0\\0 & 1 & 1\\1 & 0 & 0\end{matrix}\right]
$$
and
$$
{}^P_BJ_2 = 
\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & w\\0 & 0 & 1 & 0 & - w & 0\\1 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & 0 & 1 & 0 & 0\end{matrix}\right]
$$
we get

```
--- Kb2 ---
Matrix([
[kb + kc,  -kc,     0, 0,       0,   -kc*w],
[    -kc,   kc,     0, 0,       0,    kc*w],
[      0,    0,    ka, 0,   -ka*w,       0],
[      0,    0,     0, 0,       0,       0],
[      0,    0, -ka*w, 0, ka*w**2,       0],
[  -kc*w, kc*w,     0, 0,       0, kc*w**2]])

--- Kbtotal ---
Matrix([
[2*kb + 2*kc,    0,    0, 0,         0,   -2*kc*w],
[          0, 2*kc,    0, 0,         0,         0],
[          0,    0, 2*ka, 0,         0,         0],
[          0,    0,    0, 0,         0,         0],
[          0,    0,    0, 0, 2*ka*w**2,         0],
[    -2*kc*w,    0,    0, 0,         0, 2*kc*w**2]])
```

which matches the result from the paper
$$
K_b = 
\left[\begin{matrix}2 k_b + 2 k_c & 0 & 0 & 0 & 0 & - 2 k_c w\\0 & 2 k_c & 0 & 0 & 0 & 0\\0 & 0 & 2 k_a & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 2 k_a w^{2} & 0\\- 2 k_c w & 0 & 0 & 0 & 0 & 2 k_c w^{2}\end{matrix}\right]
$$

### Q2

To consider a soft finger contact model, we modify $H$ to include torque along the (local) $z$-axis
$$
H = 
\left[\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 1\end{matrix}\right]
$$
Adding also the $C_{tip}$ matrix to the total compliance matrix, we obtain
$$
K_b = \left[\begin{matrix}4 k_a & 0 & 0 & 0 & 0 & - 2 k_a w\\0 & 2 k_a & 0 & 0 & 0 & 0\\0 & 0 & 2 k_a + 2 k_q & - 2 k_q & 0 & 0\\0 & 0 & - 2 k_q & 2 k_q & 0 & 0\\0 & 0 & 0 & 0 & 2 w^{2} \left(k_a + k_q\right) & 0\\- 2 k_a w & 0 & 0 & 0 & 0 & 2 k_a w^{2}\end{matrix}\right]
$$

### Q3

$$
K_j = 
\left[\begin{matrix}0 & 0 & 0 & 0 & 0 & 0\\0 & - 2 f_n & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & - 2 f_n w & 0\\0 & 0 & 0 & 0 & 0 & - 2 f_n w \left(w + 1\right)\end{matrix}\right]
$$

In particular, we obtain a tilting moment, due to geometry, that would tend to unstabilize the object (but the final effect depends also on the value of $K_b$).

### Q4

By adding the translation (in the local $\{l_i,m_i,n_i\}$ reference frame) produced by the infinitesimal rotation in `diff_i`(considered also to be in the local reference frame), we obtain a new $\Delta J$ .

This result in a new geometric stiffness
$$
K_j = 
\left[\begin{matrix}0 & 0 & 0 & 0 & 0 & 0\\0 & - 2 f_n & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 2 f_n \left(R - w\right) & 0\\0 & 0 & 0 & 0 & 0 & 2 f_n \left(R - w\right) \left(w + 1\right)\end{matrix}\right]
$$
If we consider the torque produced by this geometric effect around the body's $z$-axis,
$$
\tau_{b,z}^j = 2 f_n \left(R - w\right) \left(w + 1\right) \, \delta \theta_z^2
$$
we can see that it is positive (i.e. restoring torque) for $R>w$. So, for an object that is comparatively thin compared to the radius of the fingertips, the geometric effect is also stabilizing.