
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#Working-with-a-3D-large-eddy-simulation-of-shallow-convection" data-toc-modified-id="Working-with-a-3D-large-eddy-simulation-of-shallow-convection-0.1"><span class="toc-item-num">0.1&nbsp;&nbsp;</span>Working with a 3D large eddy simulation of shallow convection</a></span><ul class="toc-item"><li><span><a href="#The-simulation" data-toc-modified-id="The-simulation-0.1.1"><span class="toc-item-num">0.1.1&nbsp;&nbsp;</span>The simulation</a></span></li><li><span><a href="#The-dataset-----netccdf" data-toc-modified-id="The-dataset-----netccdf-0.1.2"><span class="toc-item-num">0.1.2&nbsp;&nbsp;</span>The dataset  -- netccdf</a></span></li><li><span><a href="#liquid-water-cross-section-at-1-km" data-toc-modified-id="liquid-water-cross-section-at-1-km-0.1.3"><span class="toc-item-num">0.1.3&nbsp;&nbsp;</span>liquid water cross section at 1 km</a></span></li><li><span><a href="#zoom-in-on--the-top-left-corner" data-toc-modified-id="zoom-in-on--the-top-left-corner-0.1.4"><span class="toc-item-num">0.1.4&nbsp;&nbsp;</span>zoom in on  the top left corner</a></span></li><li><span><a href="#Get-a-vertical-cross-section-along-y-=-2km" data-toc-modified-id="Get-a-vertical-cross-section-along-y-=-2km-0.1.5"><span class="toc-item-num">0.1.5&nbsp;&nbsp;</span>Get a vertical cross section along y = 2km</a></span></li><li><span><a href="#Find-the-vapor-mixing-ratio-along-this-cross-section" data-toc-modified-id="Find-the-vapor-mixing-ratio-along-this-cross-section-0.1.6"><span class="toc-item-num">0.1.6&nbsp;&nbsp;</span>Find the vapor mixing ratio along this cross section</a></span></li><li><span><a href="#For-Wednesday" data-toc-modified-id="For-Wednesday-0.1.7"><span class="toc-item-num">0.1.7&nbsp;&nbsp;</span>For Wednesday</a></span></li></ul></li></ul></li><li><span><a href="#vertical-cross-section-of-temperature" data-toc-modified-id="vertical-cross-section-of-temperature-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>vertical cross section of temperature</a></span></li><li><span><a href="#Vertical-cross-section-of-relative-humidity" data-toc-modified-id="Vertical-cross-section-of-relative-humidity-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Vertical cross section of relative humidity</a></span><ul class="toc-item"><li><span><a href="#note-the-bimodal-RH-distribution----boundary-layer-air-is-different-from-free-atmos" data-toc-modified-id="note-the-bimodal-RH-distribution----boundary-layer-air-is-different-from-free-atmos-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>note the bimodal RH distribution -- boundary layer air is different from free atmos</a></span></li></ul></li><li><span><a href="#make-a-palette-that-limits-the-range-of-colors" data-toc-modified-id="make-a-palette-that-limits-the-range-of-colors-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>make a palette that limits the range of colors</a></span></li><li><span><a href="#Vertical-profiles-of-mean-and-standard-deviation-for-rh" data-toc-modified-id="Vertical-profiles-of-mean-and-standard-deviation-for-rh-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Vertical profiles of mean and standard deviation for rh</a></span></li></ul></div>

# ## Working with a 3D large eddy simulation of shallow convection

# ### The simulation
# 
# * Objective: compare a single column of a GCM with large eddy simlations for three different cloud types (stratus, stratocumulus, trade cumulus).  Basic approach is to run the same case studies twice.  Once using the physics from a single column of several global climate models:
# 
# [GCM paper](https://doi-org.ezproxy.library.ubc.ca/10.1002/2013MS000246)
# 
# and also with several large eddy models:
# 
# [LES paper](https://agupubs-onlinelibrary-wiley-com.ezproxy.library.ubc.ca/doi/abs/10.1002/jame.20025)
# 
# * We started with the trade cumulus simulation, then perturbed it by raising the temperature to 300 K and 301 K.
# 
# http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/cgils_ctl_s6_synthetic_albedo.mp4
# 
# http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/cgils_sst_300K_synthetic_albedo.mp4
# 
# http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/cgils_sst_301K_synthetic_albedo.mp4

# ###  The dataset  -- netccdf
# 
# [An example of reading a netCDF4 file ](http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html)
# 

# In[1]:


import glob
from netCDF4 import Dataset
import numpy as np
from a405.utils.ncdump import ncdump
from a405.utils.data_read import download

do_download = True
if do_download:
    root = 'https://clouds.eos.ubc.ca/phaustin/a405'
    the_file = 'ENT_CGILS_CTL_S6_3D_384x384x194_25m_1s_96_0000014160.nc'
    out = download(the_file, root=root)
    
the_file = glob.glob("*CTL*")[0]
with Dataset(the_file,'r') as ncin:
    ncdump(ncin)


# ### liquid water cross section at 1 km
# 
# 

# In[2]:


def get_var(the_file,varname):
    with Dataset(the_file,'r') as ncin:
         out=ncin.variables[varname][...]
         x = ncin.variables['x'][...]
         y = ncin.variables['y'][...]
         z = ncin.variables['z'][...]
         out = out.squeeze()  #remove the time dimension, since we only have one timestep
    return x,y,z,out
x,y,z,qn = get_var(the_file, 'QN')
print(the_file)
print(qn.max())


# In[3]:


#
#  find the index for z = 1000 meters
#

level = np.searchsorted(z, 1000)
print(level)


# In[4]:


#
# get the cloud liquid water at 1000 m
#
horiz_cross_sec = qn[level,:,:]
#
# find the cross section cloud fraction
#
cloud_frac=np.sum(horiz_cross_sec > 0)/horiz_cross_sec.size
print('cloud fraction: {:5.3f}'.format(cloud_frac))


# In[5]:


from matplotlib import pyplot as plt
plt.close('all')
fig,ax =plt.subplots(1,1,figsize=(10,10))
whole_scene=ax.imshow(horiz_cross_sec)
cax=plt.colorbar(whole_scene,ax=ax)
cax.set_label('liquid water content (g/kg)')
title = 'horizontal qn cross section at z=1000 m'
ax.set_title(title)


# ### zoom in on  the top left corner
# 
# Switch from [imshow](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow) to 
# [pcolormesh](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.pcolormesh) so we can orient the axes along model x,y, and z coordinates.  Note that if y is north/south (north up), then imshow plots the image upside down.

# In[6]:


#
# it helps in checking your orientation to make the selection
# have different numbers of rows and columns
#
end_col = 200
end_row = 180
fig,ax =plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],y[:end_row],horiz_cross_sec[:end_row,:end_col])
ax.set(xlabel='distance east',ylabel='distance north')
cax = plt.colorbar(image,ax=ax)
cax.set_label('liquid water content (g/kg)')
ax.set_title('zoomed horiz qn cross section at z=1000 m')


# ### Get a vertical cross section along y = 2km

# In[7]:


row_number_2000 = np.searchsorted(y,2000)  #(y index of 80)
print(f'y= 2km occurs at index {row_number_2000}')


# In[8]:


vert_cross_sec = qn[:,row_number_2000,:end_col]
print(f'shape of cross section along y = 2000 m: {vert_cross_sec.shape}')
fig,ax = plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec[:,:end_col])
cax = plt.colorbar(image,ax=ax)
cax.set_label('liquid water mixing ratio qn (g/kg)')
ax.set_title('vertical qn cross section along y=2 km')
ax.set(xlabel='distance east (m)',ylabel='height (m)');


# ### Find the vapor mixing ratio along this cross section

# In[10]:


plt.close('all')
x,y,z,qv = get_var(the_file, 'QV')
vert_cross_sec = qv[:,row_number_2000,:end_col]
fig,ax = plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec[:,:end_col])
cax = plt.colorbar(image,ax=ax)
cax.set_label('water vapor mixing ratio qv (g/kg)')
ax.set_title('vertical qv cross section along y=2 km');


# ### For Wednesday
# 
# 1\.  Read Thompkins chapter 4 parameterization notes through section 4.7.1
# 
# 
# 2\.  Hand in a notebook that adds cells to cgilsI.ipynb to:
# 
#      * use pcolormesh to plot a vertical cross section of the relative humidity for along y=2 km, x= 0-5 km
#      
#      * use plot to plot a vertical profile of the horizontal mean RH in for this cross section as a function of height
#      
#      * use plot to plot a vertical profile of the horizontal standard deviation of RH as a function of height

# # vertical cross section of temperature

# In[11]:


from a405.thermo.thermlib import find_rsat, find_esat
from a405.thermo.constants import constants as c
x,y,z,temp = get_var(the_file,'TABS')
x,y,z,press = get_var(the_file,'p')
vert_cross_sec_temp = temp[:,row_number_2000,:end_col]
plt.close('all')
fig,ax = plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec_temp)
cax = plt.colorbar(image,ax=ax)
cax.set_label('temperature (K)')
ax.set_title('vertical temp cross section along y=2 km');


# # Vertical cross section of relative humidity
# 
# use [broadcasting](https://eli.thegreenplace.net/2015/broadcasting-arrays-in-numpy/) 
# so I don't need to write a loop

# In[12]:


esat = find_esat(temp)
print(f'esat shape: {esat.shape}')
#
# broadcasting: add dimensions to the 1D pressure vectory so
# that it can be used in the denominator 
#
press = press[:,np.newaxis,np.newaxis]*100.  #convert to Pa
rsat = c.eps*esat/(press - esat)*1.e3  #convert to g/kg`
rh = qv/rsat

fig,ax = plt.subplots(1,1,figsize=(10,10))
vert_cross_sec_rh = rh[:,row_number_2000,:end_col]
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec_rh)
cax = plt.colorbar(image,ax=ax)
cax.set_label('relative humidity')
ax.set_title('relative humidity cross section along y=2 km');


# ## note the bimodal RH distribution -- boundary layer air is different from free atmos

# In[13]:


fig, ax = plt.subplots(1,1)
ax.hist(vert_cross_sec_rh.flat);


# # make a palette that limits the range of colors
# 
# below 0.81 is gray, between 0.81 and 0.85 is black, more than 0.98 is red
# 
# Some links about colors:
# 
# * [matplotlib palettes](https://matplotlib.org/examples/color/colormaps_reference.html)
# 
# * [xkcd color survey](https://blog.xkcd.com/2010/05/03/color-survey-results/)
# 
# * [xkcd colors from matplotlib](https://seaborn.pydata.org/generated/seaborn.xkcd_palette.html)
# 
# * [wikipedia article on RGB colors](https://en.wikipedia.org/wiki/RGB_color_model)
# 
# * use [a masked array](http://www.scipy-lectures.org/intro/numpy/numpy.html#masked-arrays) to eliminate some values
# 

# In[14]:


from matplotlib.colors import Normalize
fig,ax=plt.subplots(1,1,figsize=(10,10))
pal = plt.get_cmap('plasma')
pal.set_bad('0.75') #75% grey
pal.set_over('r')
pal.set_under('k')
vmin= 0.85
vmax= 0.98
#
#mask relative humidities > 1
#using a masked array
#
import numpy.ma as ma
mask = vert_cross_sec_rh < 0.81
ma_rh_vert = ma.array(vert_cross_sec_rh, mask = mask)
the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
image=ax.pcolormesh(x[:end_col],z,ma_rh_vert,cmap=pal,norm=the_norm)
cax=plt.colorbar(image,ax = ax,extend='both')
cax.set_label('RH (fraction)')
ax.set(xlabel='along x (meters)',ylabel='height (m)')
ax.set(ylim=[0,2500]);


    


# # Vertical profiles of mean and standard deviation for rh
# 
# focus on bottom 2500 meters

# In[15]:


end_z = np.searchsorted(z,2500)
rh_subset_bottom= rh[:end_z,:end_row,:end_col]


# the shape is (z,y,x) = axes (0,1,2)  so average over axes 1 and 2

# In[16]:


mean_rh_xy = rh_subset_bottom.mean(axis=(1,2))
sd_rh_xy = rh_subset_bottom.var(axis=(1,2))**0.5
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,10))
ax1.plot(mean_rh_xy*100,z[:end_z])
ax2.plot(sd_rh_xy*100.,z[:end_z])
ax1.set(ylabel='height (m)',xlabel='mean rh (percent)')
ax2.set(xlabel='rh standard deviation (per cent)');


# In[19]:


mask = rh_subset_bottom > 1
ma_rh_subset = ma.array(rh_subset_bottom,mask=mask)
#help(ma_rh_subset.mean)
rh_xyavg = ma_rh_subset.mean(axis = (1,2))


# In[20]:


rh_xyavg

