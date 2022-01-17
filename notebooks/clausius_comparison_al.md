---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Clausius-Clapeyron Solver

This is a notebook that adapts the rootfind notebook to solve equation (8) of the Clausius-Clapeyron equation notes:

$$
l_v = (\phi^*_v - \phi_l)T\tag{PA8}
$$

For the saturation vapor pressure esat
and compare with Thompkins 2.13 or 2.15 for 10 temperatures between 0 and 20 degrees C. (Hint: in this case, the temperature T is fixed, the variable you are trying to guess is the the vapor pressure esat, which is buried inside $\phi^*_v$.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
```

*From PA's Clausius-Clapeyron notes, we get expressions for $\phi^*_v$ and $\phi_l$*

$$
\phi^*_v = c_{pv}log\frac{T}{T_p} - R_vlog\frac{e_{sat}}{e_{s0}} + \phi_0\tag{PA10}
$$

$$
\phi_l = c_llog\frac{T}{T_p}\tag{PA11}
$$

*Where:*

$c_{pv} = 1870 \text{ J/kgK}$ is the constant-pressure heat capacity of water vapor

$T_p = 273.16 \text{ K}$ is the triple point temperature

$R_v=461.5 \text{ J/kgK}$ is the ideal gas constant for water vapor

$e_{s0}=611Pa$ is the saturation vapor pressure at the triple point

$\phi_0 = \frac{l_{v0}}{T_p}$ is the entropy per unit mass required to evaporate/condense water at the triple point

$l_v = 2.5e6 \text{ J/kg}$ is the latent heat of vaporization

$c_l = 4187 \text{ J/kgK}$ is the heat capacity of liquid water

*and $e_{sat}$ is the quantity we want to excavate.* 

```{code-cell} ipython3
# constants
conts = dict(
    cpv=1870,  # J/kgK
    Tp=273.16,  # K
    Rv=461.5,  # J/K
    es0=611,  # Pa
    lv=2.5e6,  # J/kg
    cl=4187,  # J/kgK
)
#T = np.linspace(0, 20, 10) + 273.15  # K
```

```{code-cell} ipython3
def do_phi_star_v(esat, T, cpv, Tp, Rv, es0, lv, cl):
    """
    returns the entropy per unit mass of the vapor phase as per Austin eq10
    """
    phi_star_v = cpv * np.log(T / Tp) - Rv * np.log(esat / es0) + lv / Tp
    return phi_star_v

def do_phi_lq(T, cpv, Tp, Rv, es0, lv, cl):
    """
    returns the entropy per unit mass of the liquid phase as per Austin eq11
    """
    phi_lq = cl * np.log(T / Tp)
    return phi_lq

def do_eq_8(esat, *args):
    """
    solves Austin eq 8, to be fed into rootfinder to get esat
    """
    T = args[0]
    lv = args[1]["lv"]
    return (do_phi_star_v(esat, **args[1]) - do_phi_lq(**args[1])) * T - lv
```

```{code-cell} ipython3
# find esat using the brentq rootfinder

# set search limits on esat
esat_min = 1e-4
esat_max = 1e5
temps = np.linspace(0, 20, 11) + 273.15

# find the roots and show results
print("  Temp    |    e_sat \n", "-"*20)
for temp in temps:
    esat = optimize.zeros.brentq(do_eq_8, esat_min, esat_max, args=(temp, consts))
    print(temp, "K  | ", round(esat / 10), "HPa" )
```

### Compare with Thompkins 2.13 and 2.15

*Approximate form of the CC equation:*

$$
e_s = e_{s0}exp\left[\frac{L_v}{R_v}\left(\frac{1}{T_0} - \frac{1}{T}\right)\right]\tag{AT2.13}
$$

```{code-cell} ipython3
def do_T213(T, cpv, Tp, Rv, es0, lv, cl):
    """
    Solves Thomkins 2.13 for saturation vapor pressure
    """
    return es0 * np.exp((lv / Rv) * (Tp ** -1 - T ** -1))


for temp in temps:
    print(do_T213(temp, **conts))
```

*Or Bolton's formula:*

$$
e_s = 611.2exp\left(\frac{17.67\cdot T}{T + 243.5}\right)\tag{AT2.15}
$$

```{code-cell} ipython3
def do_T215(T):
    """
    calculates esat using Bolton's formula (Thompkins 2.15)
    """
    Tc = T - 273.15
    return 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))

for temp in temps:
    print(do_T215(temp))
```

```{code-cell} ipython3

```
