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

# Day 25 Problems - Andrew's Solution

## Problem 1

Given the critical supersaturation from the kohler notes:

$$
SS=S^∗−1=\left(\frac{4a^3}{27b}\right)^{1/2}
$$

show that this implies, for (NH4)2SO4, density $\rho_{aer} = 1775 kgm^{−3}$, van hoft $i=3$, that:

$$
S^∗−1 \approx 1.54\cdot10^{−12} m^{−0.5}_{aer}
$$

where $m_{aer}$ is the ammonium sulphate aerosol mass in kg.

Note that this is why a cloud chamber can get the aerosol mass distribution from a series of saturation and light scattering measurements as smaller and smaller aerosols are pushed over their critical supersaturation and activated.

---

+++

Substitute expressions for $a$ and $b$:

$$
a = \frac{2\sigma}{\rho_l R_v T}
$$

$$
b = \frac{imM_w}{4/3\pi \rho_s M_s}
$$

```{code-cell} ipython3
import json
from pathlib import Path
import numpy as np
from a405.thermo.constants import constants as c
```

```{code-cell} ipython3
# load the properties of ammonium sulphate from JSON file
ammonium_sulphate = json.load(open(Path("../src/a405/data/ammonium_sulphate.json"), "r"))
print(ammonium_sulphate)
```

```{code-cell} ipython3
def get_SS(T, aerosol):
    """
    finds supersturation given T (int/float), aerosol species (dict), aerosol mass (float)
    """
    a = (2 * aerosol["Sigma"]) / (c.rhol * c.Rv * T)
    b = (
        (aerosol["vanHoff"] * aerosol["Mw"])
        / ((4 / 3) * np.pi * aerosol["rho"] * aerosol["Ms"])
        * aerosol["mass"]
    )
    SS = (4 * a ** 3 / (27 * b)) ** 0.5
    return SS
```

```{code-cell} ipython3
SS = get_SS(280, ammonium_sulphate)
## my answer differs from Phil in b term. I used rho_s, he used rho_l
SS / ammonium_sulphate["mass"] ** -0.5
```

## Problem 2

Show that the expression for second derivative of the thermodynamic potential derived in the kohler stability notes:

$$
\frac{\delta^2G}{\delta r^2} = -4\pi R_v T\rho_l\left[2a - r^2\left(1 + \frac{b}{r^3}\right)\frac{3b}{r^4}\right] + 8\pi\sigma
$$

Changes sign from stable (positive) to unstable (negative) at $r_{crit}$

$$
r_{crit} = \left(\frac{3b}{a}\right)^{1/2}
$$


Hint – first show that the second derivative is zero at the critical radius. Then show that the third derivative is negative above and below the critical radius, which means that there has to be a sign change from + to -.

---

+++

As in the Kohler stability notes, drop the $r^5$ term since its small:

$$
\frac{\delta^2G}{\delta r^2} = -4\pi R_v T \rho_l \left[2a - \frac{3b}{r^2}\right] + 8\pi\sigma
$$

Plug in $r = r_{crit}$

$$
\frac{\delta^2G}{\delta r^2} = -4\pi R_v T \rho_l \left[2a - \frac{3b}{3b}a\right] + 8\pi\sigma
$$

$$
= -4\pi R_v T \rho_la + 8\pi\sigma
$$

Now plug in $a = \frac{2\sigma}{\rho_l R_v T}$

$$
-4\pi R_v T \rho_l \cdot\frac{2\sigma}{\rho_l R_v T} + 8\pi\sigma
$$

$$
= -8\pi\sigma + 8\pi\sigma
$$

$$
\boxed{=0}
$$

Calculate the third derivative:

$$
\frac{\delta^3 G}{\delta r^3} = -4\pi R_v T \rho_l\left[\frac{6b}{r^3}\right]
$$

substitute $b = \frac{imM_w}{4/3\pi \rho_s M_s}$

$$
\frac{\delta^3 G}{\delta r^3} = -4\pi R_v T \rho_l\left[\frac{6\left(\frac{imM_w}{4/3\pi \rho_s M_s}\right)}{r^3}\right]
$$

Realizing that all physical values of $r, im, \rho_l, T, M_w, M_s,\rho_s$ are $\ge 0$, we can conclude that the third deriviative  $\frac{\delta^3 G}{\delta r^3} \le 0 \space\forall\space r\square$
