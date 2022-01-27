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

# Thermodiagrams -- DIY and with Metpy

Andrew Loeppky <br>
ATSC 405 - 2022

exposition

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units
```

## Part 1 - Building a Thermodiagram from Scratch in Python

### Task 1 - Plot Pressure Values With Log Spacing

Plot horizontal lines at pressures 1000 - 100 kPa, incrementing by 100 kPa, on a log scale. This gives our plot some realistic scaling, as this is how pressure levels are actually distributed in the atmposphere. See the hydrostatic approximation notes for further detail.

```{code-cell} ipython3
pres_lvls = np.arange(10000, 110000, 10000)  # Pa
pres = np.arange(10000, 100100, 100) # Pa

# initialize figure
fig, ax = plt.subplots(figsize=(9,9)) 
plt.gca().invert_yaxis()
ax.set_yscale("symlog")
ax.set_yticks(pres_lvls, pres_lvls // 100) # convert to kPa here
ax.set_ylim(104000, 10000)

# plot pressures
[ax.axhline(p, color="k", alpha=0.3) for p in pres_lvls];
```

### Task 2 - Plot Isotherms

Put temperatures ranging from -40 to +50 $^o$C along the x axis to create a $T, p$ grid. 

**Strong Opinion:**

It is a VERY good idea to always pass values around your program in SI units (Pa, K, kg/kg, etc.), and only convert when absolutely necessary. Also, when plotting, I like to leave *everything* in SI units and only change *labels* to C and kPa as needed. 

```{code-cell} ipython3
# Import a list of constants useful for meteorology and thermodynamics. 
# Ensures consistency and prevents us from always having to look them up ourselves
from a405.thermo.constants import constants as c
print(f"{c.p0 = } Pa") # how to call pressure p0
```

```{code-cell} ipython3
temps = np.arange(3, 500, 10) 

[ax.axvline(T, color="k", alpha=0.3) for T in temps]
ax.set_xticks(temps, np.round(temps - c.Tc)) # only change LABELS to Celsius
ax.set_xlim(-40 + c.Tc, 50 + c.Tc)

display(fig)
```

### Task 3 - Overlay an Upper Air Sounding

We now have the simplest possible thermodiagram with a $p$ and $T$ axis. Let's plot a sounding! The code below downloads and saves atmospheric soundings from the University of Wyoming's website, then reads them into Python as Pandas Dataframes. Your task is to select the sounding for July 17th and plot the temperature and dewpoint on our basic thermodiagram.

```{code-cell} ipython3
# todo for AL: rewrite this lib with pooch
from a405.soundings.wyominglib import write_soundings, read_soundings
values=dict(region='naconf',year='2012',month='7',start='0100',stop='3000',station='72340')
write_soundings(values, 'littlerock')
soundings= read_soundings('littlerock')
```

```{code-cell} ipython3
my_date = (2012, 7, 17, 0)
raw_sounding = soundings["sounding_dict"][my_date]
```

```{code-cell} ipython3
# put this in a library later (AL)
def UWYO_sounding_to_SI(sounding):
    """
    Takes in a University of Wyoming sounding (pandas dataframe) and returns
    a dataframe of the same sounding converted to SI units. The new sounding
    has units:
    
    pres: Pa, pressure
    hght: m, height
    temp: K, temperature
    dwpt: K, dew point temperature
    rehl: decimal, relative humidity (ie saturation = 1.0)
    mixr: kg/kg, water vapor mixing ratio
    drct: degrees, wind direction
    sknt -> wind: m/s (change column name so it makes sense)
    thta: K, potential temperature
    thte: K, equivalent potential temperature
    thtv: K, virtual potential temperature
    
    """
    out = pd.DataFrame()
    out["pres"] = sounding["pres"] * 100
    out["hght"] = sounding["hght"]
    out["temp"] = sounding["temp"] + c.Tc
    out["dwpt"] = sounding["dwpt"] + c.Tc
    out["relh"] = sounding["relh"] / 100
    out["mixr"] = sounding["mixr"] / 1000
    out["drct"] = sounding["drct"]
    out["wind"] = sounding["sknt"] * 0.5144447
    out["thta"] = sounding["thta"]
    out["thte"] = sounding["thte"]
    out["thtv"] = sounding["thtv"]
    

    return out
```

```{code-cell} ipython3
my_date = (2012, 7, 17, 0) 

# grab the sounding from the downloaded list and convert to SI units
sounding = UWYO_sounding_to_SI(soundings["sounding_dict"][my_date])
sounding
```

```{code-cell} ipython3
# plot the sounding on our simple thermodiagram
ax.plot(sounding.temp, sounding.pres, color="k", linewidth=3)
ax.plot(sounding.dwpt, sounding.pres, color="k", linewidth=3)
display(fig)
```

### Task 4 - The Skew-T Log-P Diagram

What we see here is true of pretty much all soundings; $T$ and $T_d$ decrease with height and run off the left side of the diagram. We can improve the readability of our chart by applying a coordinate transform and creating a *Skew-T log-P* diagram. Take every point $(T, p)$ on the diagram, and shift it by a parameter $\alpha$ with:

$$
T_{skew} = T + \left(\alpha \cdot log\frac{p_0}{p}\right)
$$

Write a function that implements this equation on input arrays for $T$ and $p$, and returns the skewed temperature $T_{skew}$. Re-initialize your plot and re-do parts 2 and 3 (isotherms and the sounding), using a skew parameter $\alpha= 20$

```{code-cell} ipython3
def T_to_skewT(T, p, alpha):
    """
    skews T by parameter alpha
    """
    Tskew = T + alpha * (np.log(c.p0 / p))
    return Tskew
```

```{code-cell} ipython3
# re-initialize plot (copy-paste step 1)
fig, ax = plt.subplots(figsize=(9,9)) 
plt.gca().invert_yaxis()
ax.set_yscale("symlog")
ax.set_yticks(pres_lvls, pres_lvls // 100) # convert to kPa here
ax.set_ylim(104000, 10000)
ax.set_xlim()

# plot pressures
[ax.axhline(p, color="k", alpha=0.3) for p in pres_lvls];

# create skewed isotherms
my_temps = np.arange(3.15,503.15,10)
alpha = 20 # define my skew parameter
for mt in my_temps:
    temp_vec = np.ones_like(pres) * mt
    skew_temp = T_to_skewT(temp_vec, pres, alpha)
    ax.plot(skew_temp, pres, color="k", alpha=0.3)

ax.set_xticks(temps, np.round(temps - c.Tc))
ax.set_xlim(-40 + c.Tc, 50 + c.Tc)

# plot the sounding in skewed coords
ax.plot(T_to_skewT(sounding.temp, sounding.pres, alpha), sounding.pres, color="k", linewidth=3)
ax.plot(T_to_skewT(sounding.dwpt, sounding.pres, alpha), sounding.pres, color="k", linewidth=3);
```

### Task 5

Do it again with MetPy! See if you can get the same plot using `metpy.plots.SkewT`. The default plot should look very similar to what we just created in task 4

```{code-cell} ipython3
#### do it with metpy ####
plt.rcParams['figure.figsize'] = (9, 9)
skew = SkewT()
skew.plot(sounding.pres / 100, sounding.temp - c.Tc, color="k", linewidth=3)
skew.plot(sounding.pres / 100, sounding.dwpt - c.Tc, color="k", linewidth=3)
# note, this took 4 lines to create instead of 100+. go metpy
```

### Task 6 - Dry Adiabats

Plot dry adiabats at 10K horizontal spacing that cover the whole plot. To convert from temperature to potential temperature, use Thompkins 1.38:

$$
\theta = T\left(\frac{p_0}{p}\right)^{\frac{R_d}{c_p}}\tag{AT 1.38}
$$

Plot your adiabats in red, dashed lines with an alpha value (transparancy) of 0.4

```{code-cell} ipython3
def make_dry_adiabat(theta, press, rv=0):
    """
    Creates a dry adiabat starting at a specified potential temperature
    along a pressure sounding
    """
    power = c.Rd / c.cpd * (1. - 0.24 * rv)
    tempOut = theta * (press / c.p0) ** power
    return tempOut
```

```{code-cell} ipython3
# put our dry adiabats on the figure
for mt in my_temps:
    adiabat = make_dry_adiabat(mt, pres)
    skew_adbt = T_to_skewT(adiabat, pres, alpha)
    ax.plot(skew_adbt, pres, color="r", alpha=0.4, linestyle="--")

display(fig)
```

### Task 7 - Saturated Adiabats

Add saturated adiabats to your plot with surface temperatures of -30$^o$C to +40$^o$C in 5$^o$C increments. The calculation for $\theta_e$ is less straightforward than $\theta$, feel free to use functions from `thermlib.py` to implement them in python.

```{code-cell} ipython3
import a405.thermo.thermlib as tl
np.seterr(all='ignore') # turn off divide by zero warnings, turns out that's fine

my_temps = np.arange(233.15, 313.15,5)
for mt in my_temps:
    rsat = tl.find_rsat(mt, c.p0)
    thetaes = tl.find_thetaet(mt, rsat, mt, c.p0)
    T_sat = [tl.find_Tmoist(thetaes, p, use_theta=True) for p in pres]
    T_sat_skew = T_to_skewT(T_sat, pres, alpha) 
    ax.plot(T_sat_skew, pres, color="b", alpha=0.4, linestyle="--")
    
display(fig)
```

### Task 8 - Mixing Lines

Write a function to calculate lines of constant mixing ratio, or use one from thermlib. Make 20 contours, log spaced between $10^{-5}$ and $10^{-1}$ kg/kg. Plot them in green dashed lines on your thermodiagram.

```{code-cell} ipython3
my_r = np.logspace(-5,-1,20) # mixing ratios in kg/kg
for mr in my_r: 
    T_r = [tl.find_Td(mr, p) for p in pres]
    T_r_skew = T_to_skewT(T_r, pres, alpha) 
    ax.plot(T_r_skew, pres, color="g", alpha=0.4, linestyle="--")
    

display(fig)
```

## Part 2 - Buoyancy, CAPE and Cloud Tops

We now have a pretty respectable looking thermodiagram. Let's use it to make some predictions about how clouds should develop on the day of our sounding. 

### Task 9 - Buoyancy 

To predict vertical motion in our sounding, we need to know the force/mass supplied by the difference in density, which we express as the *virtual temperature* $T_v$:

$$
<\text{equation for virtual temperature}>
$$


$$
B = g\left(\frac{T_{v,a}-T_{v,env}}{T_{v,env}}\right)\tag{AT 1.61*}
$$


+++

## Task xx

We can do this entire thing in 50 lines of code with metpy

```{code-cell} ipython3
#### do it with metpy ####

# this is all from the tutorial at:
# https://unidata.github.io/MetPy/latest/examples/Advanced_Sounding.html#sphx-glr-examples-advanced-sounding-py

p = sounding['pres'].values * units.Pa
T = sounding['temp'].values * units.kelvin
Td = sounding['dwpt'].values * units.kelvin
wind_speed = sounding['wind'].values * units('m/s')
wind_dir = sounding['drct'].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir)

fig2 = plt.figure(figsize=(9, 9))
add_metpy_logo(fig2, 115, 100)
skew = SkewT(fig2, rotation=45)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot.
skew.plot(p, T, 'r')
skew.plot(p, Td, 'g')
skew.plot_barbs(p, u, v)
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 60)

# Calculate LCL height and plot as black dot. Because `p`'s first value is
# ~1000 mb and its last value is ~250 mb, the `0` index is selected for
# `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
# i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
# should be selected.
lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

# Calculate full parcel profile and add to plot as black line
prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
skew.plot(p, prof, 'k', linewidth=2)

# Shade areas of CAPE and CIN
skew.shade_cin(p, T, prof, Td)
skew.shade_cape(p, T, prof)

# Add the relevant special lines
skew.plot_dry_adiabats(t0=my_temps * units.kelvin)
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

# Show the plot
plt.show()
```
