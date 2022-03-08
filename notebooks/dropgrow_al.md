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

# Droplet Growth - Andrew's Solution

<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Update:---Bjørn's-reported-aero.Mw-error-fixed-(see-#pha-2018/3/19)" data-toc-modified-id="Update:---Bjørn's-reported-aero.Mw-error-fixed-(see-#pha-2018/3/19)-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Update:   Bjørn's reported aero.Mw error fixed (see #pha 2018/3/19)</a></span></li><li><span><a href="#Use-a405.dropgrow.aerolib-to-find-equilibrium-drop-radii" data-toc-modified-id="Use-a405.dropgrow.aerolib-to-find-equilibrium-drop-radii-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Use a405.dropgrow.aerolib to find equilibrium drop radii</a></span></li><li><span><a href="#Read-in-the-inital-conditions-from-a-json-file" data-toc-modified-id="Read-in-the-inital-conditions-from-a-json-file-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Read in the inital conditions from a json file</a></span></li><li><span><a href="#Calculate-the-lognormal-aerosol-mass-distribution-and-get-the-number-concentration-in-each-of-30-bins" data-toc-modified-id="Calculate-the-lognormal-aerosol-mass-distribution-and-get-the-number-concentration-in-each-of-30-bins-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Calculate the lognormal aerosol mass distribution and get the number concentration in each of 30 bins</a></span><ul class="toc-item"><li><span><a href="#Find-the-number-of-aerosols-in-each-bin" data-toc-modified-id="Find-the-number-of-aerosols-in-each-bin-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>Find the number of aerosols in each bin</a></span></li><li><span><a href="#Find-the-equilibrium-radius-for-each-of-the-30-aerosol-masses" data-toc-modified-id="Find-the-equilibrium-radius-for-each-of-the-30-aerosol-masses-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Find the equilibrium radius for each of the 30 aerosol masses</a></span><ul class="toc-item"><li><span><a href="#Calculating-the-equilibrium-size-distribution-for-unactivated-aerosols" data-toc-modified-id="Calculating-the-equilibrium-size-distribution-for-unactivated-aerosols-4.2.1"><span class="toc-item-num">4.2.1&nbsp;&nbsp;</span>Calculating the equilibrium size distribution for unactivated aerosols</a></span></li></ul></li></ul></li></ul></div>

+++

## Use `a405.dropgrow.aerolib` to find equilibrium drop radii

```{code-cell} ipython3
import a405.dropgrow.aerolib
#
# new library for aerosol functions
#
from a405.dropgrow.aerolib import lognormal,create_koehler,find_koehler_coeffs
import numpy as np
from a405.utils.helper_funs import make_tuple, find_centers
from a405.thermo.rootfinder import find_interval, fzero
import json
import importlib_resources as ir  #pip install importlib_resources
import pandas as pd
#
# load pprint to print the yaml input
#
import pprint
pp = pprint.PrettyPrinter(indent=4)
```

## Read in the inital conditions from a json file

```{code-cell} ipython3
#
# use the importlib_resources module to 
# find the path to the data folder.  The file
# is located at a405/data/dropgrow.json
#
with ir.open_text('a405.data','dropgrow.json') as f:
    input_dict=json.load(f)
pp.pprint(input_dict)
```

## Calculate the lognormal aerosol mass distribution and get the number concentration in each of 30 bins

(code borrowed from aero.ipynb)

## Find the number of aerosols in each bin

```{code-cell} ipython3
mass_vals = np.linspace(-20, -16, 30)
mass_vals = 10 ** mass_vals
mu = input_dict["aerosol"]["themean"]
sigma = input_dict["aerosol"]["sd"]
totmass = input_dict["aerosol"]["totmass"]
mdist = totmass * lognormal(mass_vals, np.log(mu), np.log(sigma))
mdist = find_centers(mdist) * np.diff(mass_vals)
center_mass = find_centers(mass_vals)
ndist = mdist / center_mass
```

## Find the equilibrium radius for each of the 30 aerosol masses

(code borrowed from koehler.ipynb)

First make two debugging print functions ('info' and 'g') to easily print various floating point values.

```{code-cell} ipython3
def make_format(format_string="{:8.3g}"):
    """
    returns a function that formats with format_string
    """
    def inner_fun(value):
        return format_string.format(value)
    return inner_fun

#Now get closures from make_format and use it:

g = make_format()
info = make_format(format_string="debugging {}")

a=10
b=1.546e-23
print(info(a), g(b))
```

### Calculating the equilibrium size distribution for unactivated aerosols

Below we use the rootfinder to search on log(radius) to find the equilibrium
drop size for unactivaed aerosols on the left side of their Koehler curves

```{code-cell} ipython3
aero=make_tuple(input_dict['aerosol'])
parcel=make_tuple(input_dict['initial_conditions'])

a, b = find_koehler_coeffs(aero,parcel)

#
# sanity check
#
m=1.e-18
Scrit=(4.*a**3./(27.*b*m))**0.5;
rcrit = (3.*m*b/a)**0.5
print(("for aerosol with mass = {} kg, "
       "SScrit,rcrit are {:8.3g}, {:8.3g} microns")
        .format(m,Scrit,rcrit*1.e6))
```

```{code-cell} ipython3
koehler_fun = create_koehler(aero,parcel)

def find_diff(logr,S_target,m):
    """
    zero function for rootfinder
    """
    r = np.exp(logr)
    return S_target - koehler_fun(r,m)

S_target = parcel.Sinit
logr_start = np.log(0.1e-6)

initial_radius = []
for mass in center_mass:
    brackets = np.array(find_interval(find_diff,logr_start,S_target,mass))
    left_bracket, right_bracket = np.exp(brackets)*1.e6  #get brackets in microns for printing
    equil_rad = np.exp(fzero(find_diff,brackets,S_target,mass))
    
    Scrit=(4.*a**3./(27.*b*mass))**0.5
    
    initial_radius.append(equil_rad)
    print((f'mass = {mass:6.3g} kg\n'
           f'left bracket = {left_bracket:8.3e} microns\n'
           f'right bracket={right_bracket:8.3e} microns\n'
           f'critical supersaturation: {Scrit:6.3g}'))
    print(f'equlibrium radius at S={S_target} is {equil_rad*1.e6:5.3f} microns\n')
   
```

## Part 1

Create a plot of the ratio $n_s/n_w$ as a function of equilibrium radius for the 30 aerosols in the mass distribution. Do this twice, **A** once assuming that:

$$
n_w = \frac{4}{3}\pi r^3 \rho_l/M_w
$$

and **B**, using an exact calculation for $n_w$ – i.e. find the part of the volume that is due to the aerosol (just its mass/(aerosol density)) and subtract that off from the drop volume before calculating nw (assuming that the aerosol density is $1775 kg/m^3$).

```{code-cell} ipython3
from a405.thermo.constants import constants as c
from matplotlib import pyplot as plt
```

```{code-cell} ipython3
rad = np.array(initial_radius)
aer = input_dict["aerosol"]

ns = center_mass / aer["Ms"]

# Simplified version for A
nw = 4 / 3 * np.pi * rad ** 3 * c.rhol / aer["Mw"]
A_solution = ns / nw

# Exact solution B
aero_vol = center_mass / aer["rhoaero"]
nw_ext = (4 / 3 * np.pi * rad ** 3 - aero_vol) * c.rhol / aer["Mw"]
B_solution = ns / (nw_ext)

fig, ax = plt.subplots()
ax.plot(np.array(initial_radius) * 1e6, A_solution, label="Small Aerosol Appx")
ax.plot(np.array(initial_radius) * 1e6, B_solution, label="Exact Solution")
ax.set_xlim(0.02,0.2)
ax.set_xlabel("Droplet Radius ($\mu m$)")
ax.legend();
```

*Unsurprisingly, when there is relatively little water compared to aerosol, neglecting the volume occupied by the aerosol introduces a larger bias. This approximation is fine for large droplets or tiny aerosol particles*
