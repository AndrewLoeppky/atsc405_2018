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

# Constant Mixing Ratio - Andrew's Solution

+++ {"toc": true}

## Table of Contents

<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#plot-your-sounding" data-toc-modified-id="plot-your-sounding-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>plot your sounding</a></span><ul class="toc-item"><li><span><a href="#set-sounding_folder-to-the-name-of-your-downloaded-sounding-folder" data-toc-modified-id="set-sounding_folder-to-the-name-of-your-downloaded-sounding-folder-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>set sounding_folder to the name of your downloaded sounding folder</a></span></li></ul></li><li><span><a href="#For-Monday-9am" data-toc-modified-id="For-Monday-9am-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>For Monday 9am</a></span></li></ul></div>

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
#!conda install -y pytz
```

## Plot your sounding

+++

Set sounding_folder to the name of your downloaded sounding folder

```{code-cell} ipython3
from a405.soundings.wyominglib import read_soundings
from a405.skewT.skewlib import makeSkewDry
from a405.thermo.thermlib import convertTempToSkew
import datetime
import pytz

sounding_folder = 'delrio'
soundings= read_soundings(sounding_folder)
print(soundings.keys())
print(soundings['sounding_dict'].keys())
```

```{code-cell} ipython3
the_date=(2017,5,1,0)
the_sounding=soundings['sounding_dict'][the_date]
attributes=soundings['attributes']
#print(attributes)
fig,ax =plt.subplots(1,1,figsize=(8,8))
ax,skew = makeSkewDry(ax)
temp=the_sounding['temp']
press = the_sounding['pres']
tdew = the_sounding['dwpt']
temp_skew = convertTempToSkew(temp,press,skew)
tdew_skew = convertTempToSkew(tdew,press,skew)
ax.plot(temp_skew,press)
ax.plot(tdew_skew,press)
the_date=datetime.datetime(*the_date,tzinfo=pytz.utc)
central=pytz.timezone('US/Central')
the_date_central=the_date.astimezone(central)
title=f'Dodge City KS sounding: {str(the_date_central)}'
ax.set_title(title);
#help(convertTempToSkew)
```

```{code-cell} ipython3

```

## For Monday 9am

Check in a notebook that puts your sounding on the tephigram and draws a line of constant saturation mixing ratio 
$r_s$ = 10 g/kg between 1000 and  400 hPa.  

Hint -- you want to rootfind the temperature that satisfies Thompkins (2.20):

$$r_s = \frac{\epsilon e_s(T)}{p - e_s(T)} = 0.01\ kg/kg$$

for a range of pressures then convert the temperatures to skew coordinates.

Here is the top part of my rootfinding cell.  What's missing is

1) code that calls find_esat and calculates rsat for the guess temperature temp at pressure press and subtracts it from rsat to get the residual

2) code that passes that to the rootfinder to find the temperature such that rsat(temp,press) = rsat_target

3) code that does that for a range of pressures and converts the resulting temp,press line to skew coords and adds them to your sounding

Give it a couple of hours, and if you're stuck try to give as clear a bug report as you can on your partially completed notebook

```{code-cell} ipython3
from a405.thermo.thermlib import find_esat
from scipy import optimize
```

```{code-cell} ipython3
# code that calls find_esat and calculates rsat for the guess temperature temp 
# at pressure press and subtracts it from rsat to get the residual
def find_rsat(temp, pres):
    """
    Calculates saturation mixing ration as per Thompkins 2.20
    
    Parameters
    ----------
    
    temp : float or array_like
           Temperature of parcel (K).
    
    Returns
    -------
    
    rs : float or list
        Saturation mixing ratio (kg/kg).
    """
    epsilon = 6.22 # from Stull, Practical Meteorology 2016
    rsat = epsilon * find_esat(temp) / (pres - find_esat(temp)) # AT 2.20
    return rsat

def residual(temp, *args):
    """
    returns the residual bw target and calculated saturation ratios,
    formatted to be fed into scipy.optimize.brentq minimizer
    
    Parameters
    ----------
    rsat_target : float or array_like
                  target saturation ratio
    temp : float or array_like
           Temperature of parcel (K).
    pres : float or array_like
           pressure of parcel (Pa)
    
    Returns
    -------
    
    residual : float or array_like
               target_rsat - rsat
    """
    rsat_target = args[0]
    pres = args[1]
    return rsat_target - find_rsat(temp, pres)
```

```{code-cell} ipython3
# code that passes that to the rootfinder to find the temperature such that rsat(temp,press) = rsat_target

# do conversions and assign the target rsat
press_pa = np.asarray(press) * 1e3
rsat_target = 0.01  # kg/kg

# use the rootfinder in a list comp to find the temperature coordinate which rsat occurs
# at each pressure in the sounding (guess temperatures bw 100 and 350 K)
T_rsat = [optimize.zeros.brentq(residual, 100, 350, args=(rsat_target, p)) for p in press_pa]
```

```{code-cell} ipython3
# code that does that for a range of pressures and converts the resulting temp,press
# line to skew coords and adds them to your sounding
T_rsat_c = np.asarray(T_rsat) - 273
T_rsat_skew = convertTempToSkew(T_rsat_c ,press,skew)
ax.plot(T_rsat_skew, press)
display(fig)
```
