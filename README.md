# Delta-v calculation

## Quickstart

```bash
git clone https://github.com/Nughm3/delta-v

python3 -m venv .venv
pip install -r requirements.txt

./fetch-tle.sh
python3 delta_v.py
```

## Calculation details

The selected transfer strategy between two orbits exploits the $J_2$ perturbation to correct the RAAN of the chaser during a drift orbit (represented by $d$).

Optimal $\Delta v_{ijkm}$ for transfer between $i$-th, $j$-th debris starting at time $k$ and lasting for duration $m$ can be solved with nonlinear programming (NLP):

```math
\min_{a_d, i_d} \Delta v_{ijkm}
```

such that

```math
\Omega_c(t_{k+m}) = \Omega_j(t_{k+m})
```

where

```math
\Omega_c(t_{k+m}) = \Omega_c(t_0) + \dot \Omega_i(t_k - t-0) + \dot \Omega_d(t_{k+m} - t_k)
```

```math
\Omega_j(t_{k+m}) = \Omega(t_0) + \dot \Omega_j(t_{k+m} - t_0)
```

and $\dot \Omega_h$ in radians per second for any $h \in \{i, j, d\}$ is determined by the $J_2$ nodal precession formula:

```math
\dot \Omega_h = -\frac{3}{2} J_2 \sqrt{GM} R^2 a_h^{-\frac{7}{2}} \cos i_h
```

$\Delta v$ can be represented as a function of $a_d$ and $i_d$:

```math
\Delta v_{ijkm}(a_d, i_d) = \Delta v_{H_1,p} + \Delta v_{H_1,a} + \Delta v_{H_2,p} + \Delta v_{H_2,a}
```

```math
\Delta v_{H_1,p} = \sqrt{v_i^2 + v_{p,H_1}^2 - 2 v_i v_{p,H_1} \cos (\Delta i_i s_i)}
```

```math
\Delta v_{H_1,a} = \sqrt{v_d^2 + v_{a,H_1}^2 - 2 v_d v_{a,H_1} \cos (\Delta i_i (1 - s_i))}
```

```math
\Delta v_{H_2,p} = \sqrt{v_j^2 + v_{p,H_2}^2 - 2 v_j v_{p,H_2} \cos (\Delta i_j s_i)}
```

```math
\Delta v_{H_2,a} = \sqrt{v_d^2 + v_{a,H_2}^2 - 2 v_d v_{a,H_2} \cos (\Delta i_j (1 - s_i))}
```

where

```math
v_h = \sqrt{\frac{GM}{a_h}}
```

```math
\Delta i_i = |i_i - i_d|
```

```math
\Delta i_j = |i_j - i_d|
```

```math
s_i = \frac{1}{\Delta i_i} \arctan \left( \frac{\sin \Delta i_i}{\sqrt{a_d^3/a_i^3} + \cos \Delta i_i} \right)
```

```math
s_j = \frac{1}{\Delta i_j} \arctan \left( \frac{\sin \Delta i_j}{\sqrt{a_d^3/a_j^3} + \cos \Delta i_j} \right)
```
