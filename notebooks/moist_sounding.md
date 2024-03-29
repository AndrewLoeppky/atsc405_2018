---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
#!conda install -y pytz
```

# plot a sounding

```{code-cell} ipython3
from a405.soundings.wyominglib import read_soundings
from a405.skewT.skewlib import makeSkewDry
from a405.thermo.thermlib import convertTempToSkew
import importlib_resources as ir
import datetime
import pytz

sounding_folder = ir.files('a405.data') / 'delrio'
soundings= read_soundings(sounding_folder)
print(soundings.keys())
print(soundings['sounding_dict'].keys())
```

```{code-cell} ipython3
the_sounding=soundings['sounding_dict'][(2017,5,10,0)]
```

```{code-cell} ipython3
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels

def label_fun():
    """
    override the default rs labels with a tighter mesh
    """
    from numpy import arange
    #
    # get the default labels
    #
    tempLabels,rsLabels, thetaLabels, thetaeLabels = make_default_labels()
    #
    # change the temperature and rs grids
    #
    tempLabels = range(-40, 50, 2)
    rsLabels = [0.1, 0.25, 0.5, 1, 2, 3] + list(np.arange(4, 28, 2)) 
    return tempLabels,rsLabels, thetaLabels, thetaeLabels

skew=35.
temp=the_sounding['temp']
press = the_sounding['pres']
tdew = the_sounding['dwpt']
temp_skew = convertTempToSkew(temp,press,skew)
tdew_skew = convertTempToSkew(tdew,press,skew)

fig,ax =plt.subplots(1,1,figsize=(11,11))
corners = [-5, 35]
ax, skew = makeSkewWet(ax, corners=corners, skew=skew,label_fun=label_fun)
ax.set_title('Del Rio Texas sounding')
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,200])
ax.plot(temp_skew,press,linewidth=5)
ax.plot(tdew_skew,press,linewidth=5);
fig.savefig('delrio.pdf')
```

```{code-cell} ipython3

```
