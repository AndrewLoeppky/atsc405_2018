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

# Cape Part 1 - Solution

+++ {"toc": true}

## Table of Contents

<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Supress-autoscrolling" data-toc-modified-id="Supress-autoscrolling-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Supress autoscrolling</a></span></li><li><span><a href="#Draw-a-moist-adiabat-through-the-LFC" data-toc-modified-id="Draw-a-moist-adiabat-through-the-LFC-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Draw a moist adiabat through the LFC</a></span><ul class="toc-item"><li><span><a href="#Grab-a-Little-Rock-sounding" data-toc-modified-id="Grab-a-Little-Rock-sounding-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Grab a Little Rock sounding</a></span></li><li><span><a href="#Select-one-sounding" data-toc-modified-id="Select-one-sounding-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Select one sounding</a></span></li><li><span><a href="#Save-the-metadata-for-plotting" data-toc-modified-id="Save-the-metadata-for-plotting-2.3"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>Save the metadata for plotting</a></span></li><li><span><a href="#Convert-temperature-and-dewpoint-to-skew-coords" data-toc-modified-id="Convert-temperature-and-dewpoint-to-skew-coords-2.4"><span class="toc-item-num">2.4&nbsp;&nbsp;</span>Convert temperature and dewpoint to skew coords</a></span></li><li><span><a href="#Plot-the-sounding,-making-the-sounding-lines-thicker" data-toc-modified-id="Plot-the-sounding,-making-the-sounding-lines-thicker-2.5"><span class="toc-item-num">2.5&nbsp;&nbsp;</span>Plot the sounding, making the sounding lines thicker</a></span></li><li><span><a href="#turn-off-log(0)-warning" data-toc-modified-id="turn-off-log(0)-warning-2.6"><span class="toc-item-num">2.6&nbsp;&nbsp;</span>turn off log(0) warning</a></span></li><li><span><a href="#find-the-$\theta_{es}$-of-the--LCL" data-toc-modified-id="find-the-$\theta_{es}$-of-the--LCL-2.7"><span class="toc-item-num">2.7&nbsp;&nbsp;</span>find the $\theta_{es}$ of the  LCL</a></span><ul class="toc-item"><li><span><a href="#What-is-the-LCL-of-this-air?" data-toc-modified-id="What-is-the-LCL-of-this-air?-2.7.1"><span class="toc-item-num">2.7.1&nbsp;&nbsp;</span>What is the LCL of this air?</a></span></li></ul></li><li><span><a href="#get-the-limits-of-integration" data-toc-modified-id="get-the-limits-of-integration-2.8"><span class="toc-item-num">2.8&nbsp;&nbsp;</span>get the limits of integration</a></span></li></ul></li></ul></div>

+++

## Supress autoscrolling

```{code-cell} ipython3
%%javascript
IPython.OutputArea.prototype._should_scroll = function(lines) {
    return false;
}
```

## Draw a moist adiabat through the LFC

```{code-cell} ipython3
import numpy as np
import pandas as pd
from pprint import pformat

from a405.thermo.constants import constants as c
from a405.thermo.thermlib import convertSkewToTemp, convertTempToSkew
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels
```

```{code-cell} ipython3
from a405.soundings.wyominglib import write_soundings, read_soundings
from matplotlib import pyplot as plt
```

## Grab a Little Rock sounding

```{code-cell} ipython3
values=dict(region='naconf',year='2012',month='7',start='0100',stop='3000',station='72340')
write_soundings(values, 'littlerock')
soundings= read_soundings('littlerock')
```

```{code-cell} ipython3
soundings['sounding_dict'].keys()
```

## Select one sounding

```{code-cell} ipython3
the_time=(2012,7,17,0)
sounding=soundings['sounding_dict'][the_time]
sounding.columns
```

## Save the metadata for plotting

```{code-cell} ipython3
title_string=soundings['attributes']['header']
index=title_string.find(' Observations at')
location=title_string[:index]
print(f'location: {location}')

units=soundings['attributes']['units'].split(';')
units_dict={}
for count,var in enumerate(sounding.columns[1:]):
    units_dict[var]=units[count]
#
# use the pretty printer to print the dictionary
#
print(f'units: {pformat(units_dict)}')
```

## Convert temperature and dewpoint to skew coords

```{code-cell} ipython3
skew=35.
triplets=zip(sounding['temp'],sounding['dwpt'],sounding['pres'])
xcoord_T=[]
xcoord_Td=[]
for a_temp,a_dew,a_pres in triplets:
    xcoord_T.append(convertTempToSkew(a_temp,a_pres,skew))
    xcoord_Td.append(convertTempToSkew(a_dew,a_pres,skew))
```

## Plot the sounding, making the sounding lines thicker

```{code-cell} ipython3
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
```

```{code-cell} ipython3
fig,ax =plt.subplots(1,1,figsize=(8,8))
corners = [10, 35]
ax, skew = makeSkewWet(ax, corners=corners, skew=skew,label_fun=label_fun)
#ax,skew = makeSkewWet(ax,corners=corners,skew=skew)
out=ax.set(title=title_string)
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,400]);
l1,=ax.plot(xcoord_T,sounding['pres'],color='k',label='temp')
l2,=ax.plot(xcoord_Td,sounding['pres'],color='g',label='dew')
[line.set(linewidth=3) for line in [l1,l2]];
```

## turn off log(0) warning

```{code-cell} ipython3
np.seterr(all='ignore');
```

## find the $\theta_{es}$ of the  LCL

```{code-cell} ipython3
from a405.thermo.thermlib import find_Tmoist,find_thetaep,find_rsat,find_Tv,find_lcl,find_thetaes,find_thetaet
#
# find thetae of the surface air, at index 0
#
sfc_press,sfc_temp,sfc_td =[sounding[key][0] for key in ['pres','temp','dwpt']]
#
sfc_press,sfc_temp,sfc_td = sfc_press*100.,sfc_temp+c.Tc,sfc_td+c.Tc
```

### What is the LCL of this air?

```{code-cell} ipython3
Tlcl, plcl=find_lcl(sfc_td, sfc_temp,sfc_press)
```

```{code-cell} ipython3
print(f'found Tlcl={Tlcl} K, plcl={plcl} Pa')
```

```{code-cell} ipython3
#  convert to mks and find surface rv and thetae
#
sfc_rvap = find_rsat(sfc_td,sfc_press)
lcl_rvap = find_rsat(Tlcl,plcl)
sfc_thetae=find_thetaes(Tlcl,plcl)
press=sounding['pres'].values*100.
#
# find the index for 200 hPa pressure -- searchsorted requires
# the pressure array to be increasing, so flip it for the search,
# then flip the index.  Above 200 hPa thetae goes bananas, so
# so trim so we only have good values
#
toplim=len(press) - np.searchsorted(press[::-1],2.e4)
presstrun=press[:toplim]
#
# find temps along that adiabat
#
adia_temps= np.array([find_Tmoist(sfc_thetae,the_press) 
                      for the_press in presstrun])
adia_rvaps = find_rsat(adia_temps,presstrun)
adia_rls = sfc_rvap - adia_rvaps
env_temps = (sounding['temp'].values + c.Tc)[:toplim]
env_Td = (sounding['dwpt'].values + c.Tc)[:toplim]
height = sounding['hght'].values[:toplim]
pairs = zip(env_Td,presstrun)
env_rvaps= np.array([find_rsat(td,the_press) for td,the_press in pairs])
env_Tv = find_Tv(env_temps,env_rvaps)
adia_Tv = find_Tv(adia_temps,adia_rvaps,adia_rls)
xcoord_thetae=[]
press_hPa = presstrun*1.e-2
#
# convert the adiabatic thetae sounding to skewT coords
#
for a_temp,a_press in zip(adia_temps - c.Tc,press_hPa):
    out=convertTempToSkew(a_temp,a_press,skew)
    xcoord_thetae.append(out)
ax.plot(xcoord_thetae,press_hPa,color='r',label='thetae',linewidth=3.)
xlcl=convertTempToSkew(Tlcl - c.Tc,plcl*0.01,skew)
ax.plot(xlcl,plcl*0.01,'bo',markersize=8)
display(fig)
```

```{code-cell} ipython3
def find_buoy(adia_Tv,env_Tv):
  buoy=c.g0*(adia_Tv - env_Tv)/env_Tv
  return buoy

#
# moved find_buoy to library
#
from a405.thermo.thermlib import find_buoy
    
fig,ax =plt.subplots(1,1)
buoy=find_buoy(adia_Tv,env_Tv)
ax.plot(buoy,presstrun*1.e-2)
ax.invert_yaxis()
out=ax.set(ylabel='press (hPa)',xlabel='buoyancy (m/s^2)')
ax.grid(which='both')
```

## Get the limits of integration

```{code-cell} ipython3
#np.searchsorted finds the first crossing
zerobot=np.searchsorted(buoy,0)
#flip the array and search backwards for top crossing
zerotop=len(buoy) - np.searchsorted(buoy[::-1],0)
print(zerotop, toplim)
print('pressure levels for crossing: {} hPa {} hPa'\
      .format(press[zerobot]*0.01,press[zerotop]*0.01))
```

```{code-cell} ipython3
clipped_buoy = buoy[zerobot:zerotop]
clipped_height = height[zerobot:zerotop]
#average the levels to get layer buoyancy
layer_buoy = (clipped_buoy[:-1] + clipped_buoy[1:])/2.
cape = np.sum(layer_buoy*np.diff(clipped_height))
print('cape is {:6.2f} J/kg'.format(cape))
```

```{code-cell} ipython3
clipped_buoy = buoy[zerobot:]
clipped_height = height[zerobot:]
new_press= clipped_press[zerobot:]
#average the levels to get layer buoyancy
layer_buoy = (clipped_buoy[:-1] + clipped_buoy[1:])/2.
layer_height = (clipped_height[:-1] + clipped_height[:-1])/2.
layer_press = (new_press[:-1] + new_press[1:])/2.
cape = np.cumsum(layer_buoy*np.diff(clipped_height))
fig,ax=plt.subplots(1,1)
ax.plot(cape,layer_press*1.e-2)
ax.invert_yaxis()
```
