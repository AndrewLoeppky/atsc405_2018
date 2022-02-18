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

+++ {"toc": true}

<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Wednesday-problems" data-toc-modified-id="Wednesday-problems-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Wednesday problems</a></span><ul class="toc-item"><li><span><a href="#Problem:-Taylor-series" data-toc-modified-id="Problem:-Taylor-series-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Problem: Taylor series</a></span></li><li><span><a href="#Problem:-Term-comparison" data-toc-modified-id="Problem:-Term-comparison-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Problem: Term comparison</a></span></li><li><span><a href="#Problem-surface-energy:" data-toc-modified-id="Problem-surface-energy:-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Problem surface energy:</a></span></li></ul></li></ul></div>

+++

# Wednesday problems

+++

## Problem: Taylor series

Show using two Taylor series expansions (one for $a_w$ and one for
$\exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ]$) show how to get  Thompkins 4.15:

$S= a_w \exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] \approx \left ( 1 + \frac{a}{r} - \frac{b}{r^3} \right )$

---

+++

Assume the argument $\frac{2\sigma}{\rho_lR_vTr}$ of the exponential is small, and use a first order Maclaurin series for $e^x$

$$
e^{\frac{2\sigma}{\rho_lR_vTr}} \approx 1 + \frac{2\sigma}{\rho_lR_vTr}
$$

Assume there is much more water than salt (or the amount of *aqueous* salt, which is what matters here):

$$
a_w = \frac{n_w}{n_w + n_s} = \frac{1}{1 + n_s/n_w}
$$

Expand about $n_s/n_w = 0$ to first order:

$$
\frac{1}{1 + n_s/n_w} \approx 1 - \frac{n_s}{n_w}
$$

Substitute into the original expression for $S$ and do some algebra:

$$
S = \left(1 - \frac{n_s}{n_w}\right)\left(1 + \frac{2\sigma}{\rho_lR_vTr}\right)
$$

$$
= 1 + \frac{n_s}{n_w} - \frac{2\sigma}{\rho_lR_vTr} - \left(\frac{n_s}{n_w}\frac{2\sigma}{\rho_lR_vTr}\right)
$$

Drop the last term - since we assumed that both $n_s/n_w$ and $\frac{2\sigma}{\rho_lR_vTr}$ are small, then their product must be negligible. Also substitute:

$$
a = \frac{2\sigma}{\rho_lR_vT}
$$

$$
b = \frac{imM_w}{4/3\pi\rho_sM_s}
$$

realizing that the ratio of salt to water molecules is related to $r$ by:

$$
\frac{b}{r^3} = \frac{n_s}{n_w}
$$

plug this all in to obtain:

$$
S \approx 1 + \frac{a}{r} - \frac{b}{r^3}\square
$$

+++

## Problem: Term comparison  

For the aerosol defined in the kohler.ipynb notebook ($10^{-19}$ kg, ammonium sulphate), inside a droplet of radius $r=1\ \mu m$
      find the size of the smallest term you've kept (either $\frac{a}{r}$ or $\frac{b}{r^3}$ and compare it to
      the size of the largest term you've dropped.   Repeat this for a droplet of radius   $r=0.1\ \mu m$

```{code-cell} ipython3
from pathlib import Path
import numpy as np
import json
from src.a405.thermo.constants import constants as c
# load the properties of ammonium sulphate from JSON file
ammonium_sulphate = json.load(open(Path("../src/a405/data/ammonium_sulphate.json"), "r"))
print(ammonium_sulphate)
```

```{code-cell} ipython3
a = 2 * ammonium_sulphate["Sigma"] / (c.rhol * c.Rv * 300)
b = (
    (ammonium_sulphate["vanHoff"] * ammonium_sulphate["Mw"])
    / ((4.0 / 3.0) * np.pi * c.rhol * ammonium_sulphate["Ms"])
    * 1e-19
)

# calculate the size of each term for a 1um droplet
r = 1e-6
print("1 micron droplet\n----------------")
print(f"a / r = {a / r}")
print(f"b / r3 = {b / r ** 3}")

# and for the 0.1um droplet
r = 1e-7
print("\n0.1 micron droplet\n----------------")
print(f"a / r = {a / r}")
print(f"b / r3 = {b / r ** 3}")
```

## Problem: Surface energy

Suppose you have $r_l$ =1 g/kg of liquid water spread among N spherical drops of radius 10 $\mu m$.  Compare the surface energy of this
      population (which we've been neglecting) with the enthalpy required to vaporize it:  $l_v r_l$.  Is it negligible in comparison?

```{code-cell} ipython3
# solve for N
rl = 1e-3 # kg/kg
vwat = rl / c.rhol # total volume of water in 1kg of air
droprad = 10e-6  # m
dropvol = 4 / 3 * np.pi * droprad ** 3 # volume of 1 droplet
dropsurf = 4 * np.pi * droprad ** 2 # surface area of each droplet
N = vwat / dropvol # number of droplets in 1kg of air
sigma = ammonium_sulphate["Sigma"]

total_surface_energy = sigma * dropsurf * N
print(f"Total surface energy: {total_surface_energy} J/kg")

# latent heat to create the bulk droplets
qvap = rl * c.lv0
print(f"Total latent heat: {qvap} J/kg")

# is it negligible?
perc = total_surface_energy / (qvap + total_surface_energy) * 100
print(f"Droplet surfaces account for {round(perc, 3)}% of the total energy")
```
