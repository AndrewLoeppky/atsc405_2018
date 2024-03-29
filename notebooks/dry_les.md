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

# Dry LES

+++ {"toc": true}

## Table of Contents

<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#The-Dry-LES-dataset" data-toc-modified-id="The-Dry-LES-dataset-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>The Dry LES dataset</a></span></li><li><span><a href="#Intro-to-netcdf" data-toc-modified-id="Intro-to-netcdf-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Intro to netcdf</a></span></li><li><span><a href="#Intro-to-python-packages" data-toc-modified-id="Intro-to-python-packages-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Intro to python packages</a></span></li><li><span><a href="#Dumping-the-netcdf-metadata" data-toc-modified-id="Dumping-the-netcdf-metadata-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Dumping the netcdf metadata</a></span><ul class="toc-item"><li><span><a href="#Plot-$\theta$-profile-for-every-third-timestep-(i.e.-every-30-minutes)" data-toc-modified-id="Plot-$\theta$-profile-for-every-third-timestep-(i.e.-every-30-minutes)-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>Plot $\theta$ profile for every third timestep (i.e. every 30 minutes)</a></span></li><li><span><a href="#Color-contour-plot-of-one-level-for-realization-c1,-last-timestep" data-toc-modified-id="Color-contour-plot-of-one-level-for-realization-c1,-last-timestep-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Color contour plot of one level for realization c1, last timestep</a></span></li></ul></li></ul></div>

```{code-cell} ipython3
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
```

## The Dry LES dataset

This notebook looks at a portion of a dataset that was generated by running a large eddy simulation 10 different times with identical conditions.  The 10 realizations of temperature and pressure are stored as a single netcdf file

+++

## Intro to netcdf

See the function descriptions and tutorial at http://unidata.github.io/netcdf4-python/

+++

## Intro to python packages

a. do the following to install the course python code using [pip][1]:
   
     cd atsc405
     git fetch origin
     git reset --hard origin/master
     pip install -e .
    
   (this is called an "editable install", for reasons I'll explain in class)
   
   [1]: https://en.wikipedia.org/wiki/Pip_(package_manager)
 
b. Check the install by executing the cell below:

   If it succeeds, you should see:
   
       download case_60_10.nc: size is    499.3 Mbytes

```{code-cell} ipython3
from  a405.utils.data_read import download
the_root="http://clouds.eos.ubc.ca/~phil/docs/atsc500/data/"
the_file='case_60_10.nc'
out=download(the_file,root=the_root)
```

## Dumping the netcdf metadata

+++

Netcdf file layout:  10 groups corresponding to 10 different ensemble members.  Small slice of larger domain of LES run with surface heat flux of 60 W/m^2 and stable layer with dT/dz = 10 K/km.  Snapshots every 10 minutes for 8 hours.

We can read the metdata using

```{code-cell} ipython3
!pyncdump case_60_10.nc
```

### Plot $\theta$ profile for every third timestep (i.e. every 30 minutes)

```{code-cell} ipython3
%matplotlib inline

def make_theta(temp,press):
    """
      temp in K
      press in Pa
      returns theta in K
    """
    p0=1.e5
    Rd=287.  #J/kg/K
    cpd=1004.  #J/kg/K
    theta=temp*(p0/press)**(Rd/cpd)
    return theta

case_name='case_60_10.nc'
#
#  look at the first ensemble member
#
ensemble='c1'
with Dataset(case_name,'r') as ncin:
    #
    # grab the group variables
    #
    group = ncin.groups['c1']
    temp=group.variables['TABS'][...]
    press=ncin.variables['press'][...]
    z=ncin.variables['z'][...]
mean_temp=temp.mean(axis=(3,2))

fig,ax=plt.subplots(1,1,figsize=(10,8))
for i in np.arange(0,temp.shape[0],3):
    theta = make_theta(mean_temp[i,:],press)
    ax.plot(theta,z)
out=ax.set(xlabel=r'$\overline{\theta}$ (K)',ylabel='height (m)',
       title='LES dry run for realization 1:  surface flux=60 $W\,m^{-2}$, $\Gamma$=10 K/km')
```

```{code-cell} ipython3
temp.shape
```

### Color contour plot of one level for realization c1, last timestep

+++

1. Find the index of the level closest to 400 meters
2. Retrieve the horizontal temperature field for this realization at the last timestep

```{code-cell} ipython3
index=np.searchsorted(z,400.)
temp_400=temp[-1,index,:,:]
```

```{code-cell} ipython3
temp_diff=temp_400 - temp_400.mean(axis=(0,1))
fig,ax=plt.subplots(1,1,figsize=(10,8))
with Dataset(case_name,'r') as ncin:
    x=ncin.variables['x'][...]
    y=ncin.variables['y'][...]
cs=ax.pcolormesh(x,y,temp_diff)
cb=fig.colorbar(cs)
cb.set_label('temperature perturbation (K)',rotation=-90)
```

```{code-cell} ipython3
!pwd
```

```{code-cell} ipython3

```
