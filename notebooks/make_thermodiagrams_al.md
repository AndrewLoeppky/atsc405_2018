---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Thermodiagrams -- DIY and with Metpy

Andrew Loeppky <br>
ATSC 405 - 2022

exposition

## Learning Goals

- Make a sweet thermodiagram from scratch using matplotlib

- Create a thunderstorm forecast for the day by calculating CAPE, CIN, cloud top

- Learn about wrapper functions (metpy)

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore") # turn off warnings. DANGER!
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
my_date = (2012, 7, 6, 0) 

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

Write a function that implements this equation on input arrays for $T$ and $p$, and returns the skewed temperature $T_{skew}$. Re-initialize your plot and re-do parts 2 and 3 (isotherms and the sounding), using a skew parameter $\alpha= 30$. The pressure values/plot positions should remain unchanged.

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
alpha = 30 # define my skew parameter
for mt in my_temps:
    temp_vec = np.ones_like(pres) * mt
    skew_temp = T_to_skewT(temp_vec, pres, alpha)
    ax.plot(skew_temp, pres, color="k", alpha=0.3)

ax.set_xticks(temps, np.round(temps - c.Tc))
ax.set_xlim(-40 + c.Tc, 50 + c.Tc)

# plot the sounding in skewed coords
skew_Tsound = T_to_skewT(sounding.temp, sounding.pres, alpha)
skew_Tdsound = T_to_skewT(sounding.dwpt, sounding.pres, alpha)
ax.plot(skew_Tsound, sounding.pres, color="k", linewidth=3)
ax.plot(skew_Tdsound, sounding.pres, color="k", linewidth=3);
```

### Task 5 - Dry Adiabats

Plot dry adiabats at 10K horizontal spacing that cover the whole plot. To convert from temperature to potential temperature, use Thompkins 1.38:

$$
\theta = T\left(\frac{p_0}{p}\right)^{\frac{R_d}{c_p}}\tag{AT 1.38}
$$

Plot your adiabats in *red*, dashed lines with an alpha value (transparancy) of 0.4

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

### Task 6 - Saturated Adiabats

Add saturated adiabats to your plot with surface temperatures of -30$^o$C to +40$^o$C in 5$^o$C increments. The calculation for $\theta_e$ is less straightforward than $\theta$, feel free to use functions from `thermlib.py` to implement them in python. Plot these in *blue* with the same alpha=0.4.

```{code-cell} ipython3
import a405.thermo.thermlib as tl

my_temps = np.arange(233.15, 313.15,5)
for mt in my_temps:
    rsat = tl.find_rsat(mt, c.p0)
    thetaes = tl.find_thetaet(mt, rsat, mt, c.p0)
    T_sat = [tl.find_Tmoist(thetaes, p, use_theta=True) for p in pres]
    T_sat_skew = T_to_skewT(T_sat, pres, alpha) 
    ax.plot(T_sat_skew, pres, color="b", alpha=0.4, linestyle="--")
    
display(fig)
```

### Task 7 - Mixing Lines

Write a function to calculate lines of constant mixing ratio, or use one from thermlib. Make 20 contours, log spaced between $10^{-5}$ and $10^{-1}$ kg/kg. Plot them in *green* on your thermodiagram.

```{code-cell} ipython3
my_r = np.logspace(-5, -1, 20)  # mixing ratios in kg/kg
for r in my_r:
    T_r = [tl.find_Td(r, p) for p in pres]
    T_r_skew = T_to_skewT(T_r, pres, alpha)
    ax.plot(T_r_skew, pres, color="g", alpha=0.4, linestyle="--")

display(fig)
```

## Part 2 - Buoyancy, CAPE and Cloud Tops

We now have a pretty respectable looking thermodiagram. Let's use it to make some predictions about how clouds and possibly thunderstorms might develop on the day of our sounding. 

### Task 8 - Buoyancy 

To predict vertical motion in our sounding, we need to know the force/mass supplied by the difference in density, which we express as the *density temperature* $T_\rho$, or equivalently, *virtual temperature* $T_v$:

$$
T_\rho \approx T\frac{1 - \epsilon}{\epsilon}r_v - r_l - r_i\tag{AT 2.38}
$$

Then write the buoyant force per unit mass ($N/kg$, which reduces to $m/s^2$) is a function of the difference between the virtual temperature of the ascending parcel and the surrounding environment:

$$
B = g\left(\frac{T_{v,a}-T_{v,env}}{T_{v,env}}\right)\tag{AT 1.61*}
$$

Calculate the buoyant force as a function of pressure on an ascending parcel which originates at the surface with initial state $T_{surf}, T_{d,surf}, p_0$. Assume *pseudo-adiabatic ascent*, where all liquid and ice phase water immediately fall out of the parcel as precipitation. Add the parcel trajectory to your sounding, then calculate the buoyancy profile and plot it in a separate figure.

*HINT: See Thompkins figure 3.18 for more on parcel trajectories. The ascending air will either follow a dry adiabat or a moist adiabat depending on saturation.*

```{code-cell} ipython3
# Calculate and plot buoyancy as a function of pressure
def get_buoy(T, Td, pres):
    """
    Given an atmospheric sounding (Temperature and dew point profiles),
    calculate the buoyant force per unit mass (m/s2) of an air parcel
    originating from the surface as as it ascends. Thompkins Eqn 1.61
    """
    # initialize arrays
    T_parcel = np.zeros_like(np.asarray(pres))
    # buoy = np.zeros_like(np.asarray(pres))

    ## 1) get the trajectory of the ascending parcel ##
    Tsurf = T[0]
    Tdsurf = Td[0]
    psurf = pres[0]
    theta_surf = tl.find_theta(Tsurf, c.p0)
    Tlcl, plcl = tl.find_lcl(Tdsurf, Tsurf, psurf)

    # trajectory below LCL follows dry adiabat
    T_parcel[pres > plcl] = tl.make_dry_adiabat(theta_surf, pres[pres > plcl])

    # trajectory above LCL follows moist adiabat
    thetaes = tl.find_thetaes(Tlcl, plcl)
    T_parcel[pres <= plcl] = [
        tl.find_Tmoist(thetaes, p, use_theta=True) for p in pres[pres <= plcl]
    ]

    ## 2) calculate the buoyant force per unit mass
    # get the virtual temperatures of the sounding and pseudo-adiabatic ascent
    r_parcel = np.zeros_like(T_parcel)
    r_parcel[pres > plcl] = tl.find_rsat(Tdsurf, psurf)  # unsaturated ascent
    r_parcel[pres <= plcl] = tl.find_rsat(
        T_parcel[pres <= plcl], pres[pres <= plcl]
    )  # all water precipitates out
    Tv_parcel = tl.find_Tv(T_parcel, r_parcel)

    r_env = tl.find_rsat(Td, pres)
    Tv_env = tl.find_Tv(T, r_env)

    # Thompkins 1.61
    buoy = c.g0 * (Tv_parcel - Tv_env) / Tv_env
    return T_parcel, buoy
```

```{code-cell} ipython3
# get the buoyancy of our sounding and plot it in a new figure
T_parcel, buoy = get_buoy(sounding.temp, sounding.dwpt, sounding.pres)

# add subplot to show buoyant force/mass with height
cutoff = sounding.pres > 1000  # (Pa) virtual temp calculations fall apart at very 
                               # low moisture/pressure. Pick a suitable range
fig2, ax2 = plt.subplots()
plt.gca().invert_yaxis()
ax2.set_yscale("symlog")
ax2.set_yticks(pres_lvls, pres_lvls // 100) # convert to kPa here
ax2.set_ylim(104000, 10000)
ax2.set_xlim(-1,1)

ax2.plot(buoy[cutoff], sounding.pres[cutoff])
ax2.axvline(0, color="k", alpha=0.5, linestyle="--")
ax2.set_xlabel("Buoyancy (m/s$^2$)");
```

```{code-cell} ipython3
# Add parcel trajectory to our skewT diagram
skew_parcel = T_to_skewT(T_parcel, sounding.pres, alpha)
ax.plot(skew_parcel, sounding.pres, color="k", linewidth=2)

display(fig)
```

#### Analysis Questions

- Do you expect deep convection on the day of your sounding? How do you know?

    *Yes! (July 6/2012) As the surface heats up in the afternoon, parcels of air will be able to penetrate through the negatively buoyant zone between 900kPa to 850kPa. Above this level, parcels are positively buoyant until they reach 200kPa. This is definitely what we would call deep convection, so we expect thunderstorms later in the afternoon based on this sounding.*


- Near the top of the plot, the buoyancy becomes wildly negative. Explain why this is.

    *The tropopause is defined by a large isothermal layer (potential temperature increasing) followed by a temperature inversion. Air from the surface will not penetrate through the entire layer, confining deep convection to the troposphere.*

+++

### Task 9 - Key Pressure Levels

From the bouyancy profile, we can obtain a more useful variables often used for forecasting convection: *Convective available potential energy (CAPE)* and *convective inhibition (CIN)*. First, we need to identify key pressure levels:

- The surface: *Should be easy to find.*


- The level of free convection (LFC): *Whe point above which a parcel will spontaneously rise until it hits the tropopause (Not always super clearly defined, some soundings will alternate between positive and negative buoyancy throughout the profile. If this is the case, try a new sounding in July or August with an obvious negatively buoyant zone near the surface and a positive zone from $\approx$ 850 kPa to the tropopause. (July 6, 2012 at station 72340 works well for me).*


- The level of neutral bouoyancy (LNB): *The pressure level near the tropopause where your buoyancy curve crosses zero and the parcel begins to decelerate. Note this is sometimes also called the equilibrium level (EL)*


Write a function to identify each of these layers and mark them on your thermodiagram. It's super hard to automate this so that it works well for any general sounding, my approach is to look for buoyancy zero crossings within a pressure range that you prescribe, then pick out each variable by hand.

```{code-cell} ipython3
def get_buoy_crossing(sounding, pmin, pmax):
    """
    Returns the pressure level of a slice of a sounding nearest to neutral
    buoyancy within a specified pressure range (Pa)
    """
    # define buoy profile
    _, buoy = get_buoy(sounding.temp, sounding.dwpt, sounding.pres)
    
    # restrict profiles to within the range
    buoy_range = buoy[sounding.pres < pmax][sounding.pres > pmin]
    pres_range = sounding[sounding.pres < pmax][sounding.pres > pmin].reset_index()
    
    neutral = pres_range.pres[np.argmin(abs(buoy_range))]
    return neutral
```

```{code-cell} ipython3
# find the surface pressure
SFC = sounding.pres[0]
print(f"{SFC = } Pa")

# search for the LFC between 70000 and 80000 Pa
LFC = get_buoy_crossing(sounding, 70000, 80000)
print(f"{LFC = } Pa")

# search for the LNB between 15000 and 25000 Pa
LNB = get_buoy_crossing(sounding, 15000, 25000)
print(f"{LNB = } Pa")
```

```{code-cell} ipython3
# Mark them on our plot
ax.axhline(LNB, color="r", linestyle=":", linewidth=3)
ax.axhline(LFC, color="r", linestyle=":", linewidth=3)

display(fig)
```

### Task 10 - CAPE and CIN

CIN is defined as the energy that a parcel must be supplied to rise through the negative buoyancy zone and reach the level of free concection (paraphrased from Thompkins, see p48). We can calculate this energy by integrating buoyancy (force) over the distance separating the surface from the LFC (work = force x distance).

$$
CIN = \int_{SFC}^{LFC} B\cdot dh
$$

CAPE is the amount of energy gained by the parcel between the LFC and LNB, which we can obtain similarly to CIN, but integrating buoyancy between the LFC and LNB. 

$$
CAPE = \int_{LFC}^{LNB} B\cdot dh
$$

**Source:** Mish-mash of Thompkins 3.2 and Hobbs 1977 as referenced in the [MetPy docs](https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.cape_cin.html)

Calculate the CIN and CAPE for your sounding by numerically integrating. Note the $dh$ term is the *height difference between levels,* so you will have to find a way to interpolate either $h$ or $B$ so that you have the same number of $B$'s and $dh$'s to multiply together.

```{code-cell} ipython3
def get_CIN_CAPE(sounding, pmin, pmax):
    """
    calculates CAPE or CIN between two pressure levels
    by integrating the buoyancy profile
    """
    # get the buoyancy profile
    _, buoy = get_buoy(sounding.temp, sounding.dwpt, sounding.pres)
    
    # integration bounds
    buoy = np.asarray(buoy[pmin < sounding.pres][sounding.pres < pmax])
    ht = np.asarray(sounding.hght[pmin < sounding.pres][sounding.pres < pmax])
    
    # average the buoyancy in each layer so we have the same number ht's as buoy's
    lyr_buoy = [(buoy[i] + buoy[i + 1]) / 2 for i in range(len(buoy) - 1)]
    layer_ht = np.diff(ht)
    
    # do the integration
    CAPE_elems = lyr_buoy * layer_ht
    CAPE_CIN = np.sum(CAPE_elems)
    
    return CAPE_CIN
```

```{code-cell} ipython3
# calculate CIN and CAPE for our sounding
CIN = get_CIN_CAPE(sounding, LFC, SFC)
CAPE = get_CIN_CAPE(sounding, LNB, LFC)

print(f"{CIN  = } J/kg")
print(f"{CAPE = } J/kg")
```

```{code-cell} ipython3
# shade in CIN and CAPE on the thermodiagram
is_CIN = sounding.pres >= LFC
is_CAPE = (sounding.pres <= LFC) & (sounding.pres >= LNB)
ax.fill_betweenx(sounding.pres[is_CIN], skew_parcel[is_CIN], skew_Tsound[is_CIN], color="blue", alpha=0.3)
ax.fill_betweenx(sounding.pres[is_CAPE], skew_parcel[is_CAPE], skew_Tsound[is_CAPE], color="red", alpha=0.3)

display(fig)
```

### Task 11 - Updraft Speed and Cloud Top

Last task - we want to know the upper limit of how high we can expect convective cells do develop, i.e. how deep into the tropopause the updrafts can penetrate. To do this, lets assume that 100% of CAPE is converted into kinetic energy of the updrafts (This is a bold assumption. In reality entrainment, mixing, and aerodynamic drag all slow the parcel down). The kinetic energy per unit mass of a parcel is:

$$
E_k = \frac{1}{2}w^2
$$

Rearrange to get Thompkins 3.4:

$$
w_{max} = \sqrt{2\cdot CAPE}\tag{AT 3.4}
$$

Real updrafts in real deep convective cells have velocities something like 30m/s. Calculate the updraft speed given the CAPE of your sounding. How realistic is our zero entrainment assumption?

```{code-cell} ipython3
wmax = np.sqrt(2 * CAPE)
print(f"{wmax = } m/s")
```

Now let's find the highest possible cloud top height. Once an updraft passes the LNB, it begins to decelerate until it's kinetic energy reaches zero and the parcel stops. Use your CAPE function to find the pressure at which the parcel's CAPE is completely dissipated into the tropopause. Mark this point on your plot.

```{code-cell} ipython3
from scipy import optimize

def rootfind_CAPE(pres, *args):
    """
    finds CAPE, set up for rootfinding function to detect cloud top
    
    args = [sounding, LFC]
    """
    sounding = args[0]
    pmax = args[1]
    return get_CIN_CAPE(sounding, pres, pmax)
```

```{code-cell} ipython3
# use scipy rootfinder
CTOP = optimize.zeros.brentq(rootfind_CAPE, 1000, LNB, args=(sounding, LFC))
```

```{code-cell} ipython3
# plot cloud top
ax.axhline(CTOP, color="r", linestyle=":", linewidth=3)
display(fig)
```

## Part 2 - Repeat with MetPy

Now that we know how to painstakingly derive and plot everything by hand, we can repeat this entire analysis in like 50 lines of code with Metpy. Metpy is a python package for reading, visualizing, and performing calculations with weather data, and has a whole bunch of well-documented, maintained functions that wrap matplotlib's regular plotting package with meteorology applications in mind - like built in methods for handling unit conversions between C and K, Pa and hPa, etc. So handy!

[Metpy Home Page](https://unidata.github.io/MetPy/latest/index.html)


Install from the command line with:
```
$ mamba install metpy
```

Here's the analysis again, high speed:

```{code-cell} ipython3
#### do it with metpy ####
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units

# this is all lightly modified from the tutorial at:
# https://unidata.github.io/MetPy/latest/examples/Advanced_Sounding.html#sphx-glr-examples-advanced-sounding-py

# Task 1/2: initialize figure
fig3 = plt.figure(figsize=(9, 9))
#add_metpy_logo(fig3, 115, 100)
skew = SkewT(fig3, rotation=45)

# Task 3: read in sounding vars
p = sounding["pres"].values * units.Pa
T = sounding["temp"].values * units.kelvin
Td = sounding["dwpt"].values * units.kelvin

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot.
skew.plot(p, T, "k", linewidth=3)
skew.plot(p, Td, "k", linewidth=3)
skew.ax.set_ylim(1020, 100)
skew.ax.set_xlim(-40, 50)

# Task 5: dry adiabats
skew.plot_dry_adiabats(t0=np.arange(-40, 200, 10) * units.degC)

# Task 6: moist adiabats
skew.plot_moist_adiabats()

# Task 7: mixing lines (this one doesnt work as well as I'd hoped)
skew.plot_mixing_lines()

# Task 8/9: find key levels
# Calculate full parcel profile and add to plot as black line
prof = mpcalc.parcel_profile(p, T[0], Td[0]).to("degC")
skew.plot(p, prof, "k", linewidth=2)

# level of free convection
LFC = mpcalc.lfc(p, T, Td, prof)
skew.ax.axhline(LFC[0], color="r", linestyle=":", linewidth=3)

# LNB (metpy calls it EL, equilibrium level)
EL = mpcalc.el(p, T, Td, parcel_temperature_profile=prof, which='top')
skew.ax.axhline(EL[0], color="r", linestyle=":", linewidth=3)

# Task 10: CAPE and CIN
CAPE, CIN = mpcalc.cape_cin(p, T, Td, prof, which_lfc='bottom', which_el='top')
print(f"{CAPE = }\n{CIN = }")
skew.shade_cin(p[p > LFC[0]], T[p > LFC[0]], prof[p > LFC[0]]) # shade CIN below lfc
skew.shade_cape(p, T, prof) # shade all CAPE

# Task 11: cloud top
# (Metpy has no feature for wmax or cloud top)

# Show the plot
plt.show()
```
