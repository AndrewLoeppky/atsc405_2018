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
# get a sounding (the station code for quilayute, WA is 72797)
values=dict(region='naconf',year='2021',month='7',start='1000',stop='1000',station='72797')
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
the_sounding = soundings['sounding_dict'][(2021, 7, 10, 0)]
the_sounding
```

```{code-cell} ipython3
# use find_Tmoist to get the temps along a moist adiabat
Tsurf = the_sounding.temp[0] + 273.15
Tdsurf = the_sounding.dwpt[0] + 273.15
rtsurf = the_sounding.mixr[0]
psurf = the_sounding.pres[0] * 100
thetaet = find_thetaet(Tdsurf, rtsurf, Tsurf, psurf)
#Tmoist = [find_Tmoist(thetaet, pres) for pres in the_sounding.pres * 100]
#TmoistC = np.array(Tmoist) - 273.15
thetaet
```

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

# add the other stuff
#moist_adiabat = convertTempToSkew(TmoistC, the_sounding.pres, skew)

ax.plot(skewtemp, the_sounding["pres"], color="k", linewidth=3)
ax.plot(skewdwpt, the_sounding["pres"], color="k", linewidth=3)
#ax.plot(moist_adiabat, the_sounding.pres, color="r", linestyle=":")
ax.set_title("override")
xcorners = find_corners(corners, skew=skew)
#ax.set(xlim=xcorners, ylim=[1000, 800])
```

```{code-cell} ipython3

```
