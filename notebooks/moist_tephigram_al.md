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
<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#plot-a-sounding" data-toc-modified-id="plot-a-sounding-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>plot a sounding</a></span></li></ul></div>

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
from a405.soundings.wyominglib import read_soundings, write_soundings
from a405.skewT.skewlib import makeSkewDry
from a405.thermo.thermlib import convertTempToSkew, find_Tmoist, find_thetaet
import datetime
import pytz
#!conda install -y pytz
```

```{code-cell} ipython3
# pick a date to plot
my_month="7"
day="12"
```

```{code-cell} ipython3
# get a sounding (the station code for quilayute, WA is 72797)
values=dict(region='naconf',year='2017',month=my_month,start=day+"00",stop=day+'00',station='72797')
write_soundings(values, 'quillayute')
```

# plot a sounding

```{code-cell} ipython3
soundings= read_soundings('quillayute')
print(soundings.keys())
print(soundings['sounding_dict'].keys())
```

```{code-cell} ipython3
# grab selected sounding from the dictionary
the_sounding = soundings['sounding_dict'][(2017, int(my_month), int(day), 0)]
the_sounding
```

```{code-cell} ipython3
import a405.thermo.thermlib as tl
CtoK = 273.15  # convert kelvin to celsius and vv

# grab surface values
Tsurf = the_sounding.temp[1] + CtoK
rtsurf = the_sounding.mixr[0] / 1000
Tdsurf = the_sounding.dwpt[0] + CtoK
psurf = the_sounding.pres[0] * 100

# Find the lifting condensation level
Tlcl, plcl = tl.find_lcl(Tdsurf, Tsurf, psurf)
Tlcl -= CtoK
plcl /= 100

# Find the equivalent potential temperature at the LCL
thetaet = tl.find_thetaet(Tdsurf, rtsurf, Tsurf, psurf)

# Draw a dry adiabat from the surface T to the LCL
dabt = tl.make_dry_adiabat(Tsurf, the_sounding.pres * 100)
dabt -= CtoK

# Draw a line of constant mixing ratio from the surface Td
T_rs = np.array([tl.tinvert_rsat(Tdsurf + CtoK, rtsurf, pres) for pres in np.asarray(the_sounding.pres)[:10]])
T_rs -= CtoK

# Draw a moist adiabat from the LCL upwards
mabt = np.array([tl.find_Tmoist(thetaet, pres, use_theta=True) for pres in np.asarray(the_sounding.pres) * 100])
mabt -= CtoK
```

```{code-cell} ipython3
T_rs
```

mabt

```{code-cell} ipython3
from a405.skewT.fullskew import makeSkewWet, find_corners, make_default_labels


def label_fun():
    """
    override the default rs labels with a tighter mesh
    """
    from numpy import arange

    #
    # get the default labels
    #
    tempLabels, rsLabels, thetaLabels, thetaeLabels = make_default_labels()
    #
    # change the temperature and rs grids
    #
    tempLabels = range(-40, 50, 2)
    rsLabels = [0.1, 0.25, 0.5, 1, 2, 3] + list(np.arange(4, 28, 0.5))
    return tempLabels, rsLabels, thetaLabels, thetaeLabels


fig, ax = plt.subplots(1, 1, figsize=(12, 8))
corners = [-5, 25]
ax, skew = makeSkewWet(ax, corners=corners, skew=35, label_fun=label_fun)

# add temp and dewpoint to the diagram
skewtemp = convertTempToSkew(the_sounding["temp"], the_sounding["pres"], skew)
skewdwpt = convertTempToSkew(the_sounding["dwpt"], the_sounding["pres"], skew)
ax.plot(skewtemp, the_sounding["pres"], color="k", linewidth=3)
ax.plot(skewdwpt, the_sounding["pres"], color="k", linestyle="--", linewidth=3)

# add the extra stuff
skewlcl = convertTempToSkew(Tlcl, plcl, skew)  # the LCL
ax.scatter(skewlcl, plcl, color="r", linewidth=10)
skewmabt = convertTempToSkew(mabt, the_sounding["pres"], skew)  # moist adiabat
ax.plot(skewmabt, the_sounding["pres"], color="r", linestyle=":", linewidth=4)
skewdabt = convertTempToSkew(dabt, the_sounding["pres"], skew) # dry adiabat
ax.plot(skewdabt, the_sounding["pres"], color="r", linestyle=":", linewidth=4)
skewT_rs = convertTempToSkew(T_rs, np.asarray(the_sounding["pres"])[:10], skew) # constant r
ax.plot(skewT_rs, np.asarray(the_sounding["pres"])[:10], color="r", linestyle=":", linewidth=4)


ax.set_title("override")
xcorners = find_corners(corners, skew=skew)
ax.set(xlim=xcorners, ylim=[1050, 700]);
```

```{code-cell} ipython3
skewtemp
```

```{code-cell} ipython3
np.asarray(the_sounding.pres) * 10
```
