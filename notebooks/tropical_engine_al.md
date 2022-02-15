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
<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Carnot-cycle-problem:" data-toc-modified-id="Carnot-cycle-problem:-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Carnot cycle problem:</a></span><ul class="toc-item"><li><span><a href="#A.-get-the-surface-enthalpy-in-the-tropics-(point-A)" data-toc-modified-id="A.-get-the-surface-enthalpy-in-the-tropics-(point-A)-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>A. get the surface enthalpy in the tropics (point A)</a></span></li><li><span><a href="#B.-Lift-to-400-hPa-and-remove-80%-of-the-liquied-water" data-toc-modified-id="B.-Lift-to-400-hPa-and-remove-80%-of-the-liquied-water-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>B. Lift to 400 hPa and remove 80% of the liquied water</a></span></li><li><span><a href="#D.-descend-adiabatically-to-surface" data-toc-modified-id="D.-descend-adiabatically-to-surface-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>D. descend adiabatically to surface</a></span></li></ul></li></ul></div>

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units

import a405.thermo.thermlib as tl
import a405.thermo.constants as c
```

# Carnot cycle problem:

  Consider the following heat engine:

  1.  pseudo-adiabatic ascent in the tropical Pacific, surface pressure = 1000 hPa, where the SST is 300 K and
      the relative humidity is 80%.  Air rises to 400 hPa and loses 80% of its liquid water

  2.  Cooling at a constant pressure of 400 hPa by 20 K, with no change in total water

  3.  Descent to 1000 hPa, conserving water

  Do the following in a notebook:

  1\. Use makeSkewWet to  draw this cycle on a tephigram, using carnot.py and the midterm solution notebook as guidance.

  2\. Calculate the change in enthalpy at the surface and at 400 hPa, including both sensible and latent heat terms.

  3\. Find the total work done by the engine, and its efficiency.

  4\. Find the percentage of the total heat change that is due to addition and remove of water
     during the cycle.

+++

### Pseudo-adiabatic ascent in the tropical Pacific 

surface pressure = 1000 hPa, where the SST is 300 K and
      the relative humidity is 80%.  Air rises to 400 hPa and loses 80% of its liquid water

+++

*First, get the water content in terms of $T_d$, $r_v$ and $r_l$:*

```{code-cell} ipython3
# initial parcel state
T0 = 300 * units.kelvin
p0 = 1.0e5  * units.Pa
RH0 = 0.8

# calcs
rv0 = mpcalc.mixing_ratio_from_relative_humidity(p0, T0, RH0)  #tl.find_rsat(T0, p0) * RH0
rl0 = 0  # parcel is unsaturated, so no liquid water
rT0 = rv0 + rl0
Td0 = mpcalc.dewpoint_from_relative_humidity(T0, RH0).to("kelvin") # tl.find_Td(rv0, p0)

# print out the initial state
print(f"{T0 = }")
print(f"{rv0 = }")
print(f"{rl0 = }")
print(f"{rT0 = }")
```

```{code-cell} ipython3
# lift the parcel to 400hPa
p1 = 4.e4 * units.Pa
pres = np.linspace(p0, p1, 20) 
#prof01 = mpcalc.parcel_profile(pres, T0, Td0)

thetaet = tl.find_thetaet(Td0.magnitude, rv0.magnitude, T0.magnitude, p0.magnitude)
prof01 = np.array([tl.tinvert_thetae(thetaet, rT0, p) for p in pres.magnitude]) * units.kelvin

# rain out 80% of the liquid water
T1 = prof01[-1,0] # grab the temperature at the top of the first trajectory
rv1, rl1 = tl.find_rvrl(T1.magnitude, rv0.magnitude, p1.magnitude)
rl1 *= 0.2
rT1 = rv1 + rl1

# print out state at (1)
print(f"{T1 = }")
print(f"{rv1 = }")
print(f"{rl1 = }")
print(f"{rT1 = }")
```

```{code-cell} ipython3
# make a skewT
fig = plt.figure(figsize=(9, 9))

skew = SkewT(fig, rotation=45)
skew.ax.set_ylim(1020, 300)
skew.ax.set_xlim(-10, 40)
skew.plot_dry_adiabats(t0=np.arange(-40, 200, 10) * units.degC)
skew.plot_moist_adiabats()
skew.plot_mixing_lines(pressure=np.linspace(1020,100) * units.hPa)

# plot (1)
skew.plot(pres, prof01, "k", linewidth=3);
```

###  Cooling at a constant pressure of 400 hPa by 20 K, with no change in total water

```{code-cell} ipython3
T2 = T1 - (20 * units.kelvin)
prof12 = np.array([T1.magnitude, T2.magnitude]) * units.kelvin # gross....
p12 = np.ones(2) * pres[-1]

# find the vapor and liquid mixing ratios
rT2 = rT1 # no change in total water
rv2, rl2 = tl.find_rvrl(T2.magnitude, rT2, p1.magnitude)
rT2 = rv2 + rl2

# print out state at (2)
print(f"{T2 = }")
print(f"{rv2 = }")
print(f"{rl2 = }")
print(f"{rT2 = }")

# add trajectory to the skewT
skew.plot(p12, prof12, "k", linewidth=3)
display(fig)
```

### Descent to 1000 hPa, conserving water

```{code-cell} ipython3
## metpy doesnt let us descend parcels...
thetaes = tl.find_thetaes(T2.magnitude, p1.magnitude)
rT3 = rT2 # conserve water
prof23 = np.array([tl.tinvert_thetae(thetaes, rT3, p) for p in pres.magnitude]) * units.kelvin

T3 = prof23[0,0]
rv3, rl3 = tl.find_rvrl(T3.magnitude, rT3, p0.magnitude)
rT3 = rv3 + rl3

# print out state at (3)
print(f"{T3 = }")
print(f"{rv3 = }")
print(f"{rl3 = }")
print(f"{rT3 = }")

# plot the last line on the diagram
skew.plot(pres, prof23, "k", linewidth=3)
display(fig)
```
