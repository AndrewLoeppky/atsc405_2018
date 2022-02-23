---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Collision Coalescence - Andrew's Solution

## Question 1

Assuming that cloud condensation nuclei (CCN) are removed from the atmosphere by first serving as the centers on which cloud droplets form, and the droplets subsequently grow to form precipitation particles, estimate the residence time of a CCN in a column extending from the surface of the Earth to an altitude of 5 km. Assume that the annual rainfall is 100 cm and the cloud liquid water content is 0.30 g/kg . 

*Hint: Assume that all drops in the cloud have a radii of 10 microns and that every droplet contains exactly 1 CCN. How many CCN are in 1 kg of air? About how many kg of air are there in a 5 km column? About how many CCN are taken out by a rain rate of 1 m/year? Can you encode this information in a time constant of the form $1/N$ $dN/dt = 1/\tau$ ?*

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np
import matplotlib.pyplot as plt
```

```{code-cell} ipython3
rl = 0.3e-3  # kg/kg
RR = 1.  # meters of precip per year
droprad = 10e-6  # m

# How many CCN in 1kg of air?
kgvol = rl / c.rhol  # m3 of water per kg of air
kgN = kgvol / (4 / 3 * np.pi * droprad ** 3)

# How many kg of air in a 5km column?
# 5km roughly equals 500mb, which means half the atmosphere lies below (and above). So divide p0 by 2g
# units p0 ~ [kgm/s2 * m2], g ~[m/s2], units cancel to give kg of air in a 1m2 column
colmass = c.p0 / 2 / c.g0

# How many droplets in a column?
Ncol = colmass * kgN

# How many droplets rain out in a year?
Nyear = RR / (4 / 3 * np.pi * droprad ** 3)

# Calculate the time constant tau = N / (dn/dT)
tau = Ncol / Nyear * 365 # convert to days
print(f"The residence time for a CCN in the lower troposphere is {round(tau,2)} days")
```

## Question 2

A drop with an initial radius of $100 \mu m$ falls through a cloud containing 100 droplets per cubic centimeter that it collects in a continuous manner with a collection efficiency of 0.800. If all the cloud droplets have a radius of $10 \mu m$, how long will it take for the drop to reach a radius of 1 mm? You may assume that for the drops of the size considered in this problem the terminal fall speed $v$ (in $ms^{−1}$) of a drop of radius $r$ (in meters) is given by $v=8x10^3r$. Assume that the cloud droplets are stationary and that the updraft velocity in the cloud is negligible. 

*Hint: Integrate Thompkins equation 4.28*

---

$$
\frac{dR}{dt} = \frac{LV}{4\rho_L}\tag{AT 4.28}
$$

```{code-cell} ipython3
initrad = 100e-6  # m
dropdens = 100e6  # droplets per cubic m
coleff = 0.8  # collection efficiency
L = 4 / 3 * np.pi * 10e-6 ** 3 * dropdens * c.rhol # kg water/cubic m

# set up numerical integration
delt = 1  # s
t = 0
time = []
droprad = []

# do the integration
while initrad < 1e-3:
    v = 8e3 * initrad
    initrad += coleff * L * v / (4 * c.rhol)
    t += delt
    droprad.append(initrad)
    time.append(t)

plt.plot(np.array(time) / 60, np.array(droprad) * 1e6)
plt.xlabel("Elapsed Time (min)")
plt.ylabel("Droplet Radius ($\mu m$)");
```
