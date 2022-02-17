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

# Kohler Problem - Andrew's Solution

```{code-cell} ipython3
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from a405.thermo.constants import constants as c
```

## Part 1

Plot the Kohler curve (Thompkins equation 4.15) for the aerosol in Phil's Kohler notes eq. 6 and 7

$$
\frac{e_s^r(sol)}{e_s^\infty} \approx 1 + \frac{a}{rT} - \frac{b}{r^3}\tag{AT 4.15}
$$

Question: Where does $T$ come from in the second term? Thompkins typo?
From page 3 of the Kohler notes:

$$
a = \frac{2\sigma}{\rho_l R_v T}
$$

$$
b = \frac{imM_w}{4/3\pi \rho_s M_s}
$$

```{code-cell} ipython3
# load the properties of ammonium sulphate from JSON file
ammonium_sulphate = json.load(open(Path("../src/a405/data/ammonium_sulphate.json"), "r"))
print(ammonium_sulphate)
```

```{code-cell} ipython3
def find_S(r, aerosol, T):
    """
    calculates supersaturation S given an aerosol (dict, keys Ms, Mw, Sigma, vanHoff, density),
    temperature T (float, K), and droplet radius r (float, m)
    uses Thompkins 4.15
    """
    # get coefficients
    a = (2 * aerosol["Sigma"]) / (c.rhol * c.Rv * T)
    #b = (aerosol["vanHoff"] * aerosol["Mw"]) / (4 / 3 * np.pi * aerosol["density"] * aerosol["Ms"])
    b=(aerosol["vanHoff"]*aerosol["Mw"])/((4./3.)*np.pi*c.rhol*aerosol["Ms"])* 1e-19 
    
    # apply AT 4.15
    S = 1 + a / r - b / r ** 3
    return S
```

```{code-cell} ipython3
# make a kohler curve from 0.1 to 0.5 microns, T = 300K
rad = np.linspace(0.1, 0.5) * 1e-6
S = find_S(rad, ammonium_sulphate, 300)
```

```{code-cell} ipython3
plt.plot(rad * 1e6, S)
plt.title("Kohler Curve for Ammonium Sulphate, 300K")
plt.xlabel("Droplet Radius ($\mu$m)")
plt.ylabel("Saturation");
```

## Part 2

Use our rootfinder to find the equilibrium radius for a haze particle at a relative humidity of 90% and a temperature of 15 deg C

```{code-cell} ipython3
def rootfind_r(r, *args):
    aerosol = args[0]
    T = args[1]
    S_target = args[2]
    S = find_S(r, aerosol, T)
    return S - S_target


RH = 0.9
T = 289  # K
r90 = optimize.zeros.brentq(
    rootfind_r, 0.03 * 1e-6, 0.5 * 1e-6, args=(ammonium_sulphate, T, RH)
)
print(f"Equilibrium radius for a haze particle at {RH * 100}% RH and {T} K: {round(r90 * 1e6,3)} um")
```
