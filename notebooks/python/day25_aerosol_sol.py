
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Given-the-critical-supersaturation-from-the-kohler-notes:" data-toc-modified-id="Given-the-critical-supersaturation-from-the-kohler-notes:-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Given the critical supersaturation from the kohler notes:</a></span><ul class="toc-item"><li><span><a href="#specify-the-aerosol-properties" data-toc-modified-id="specify-the-aerosol-properties-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>specify the aerosol properties</a></span></li><li><span><a href="#Turn-the-aerosol-dictionary-into-a-named_tuple" data-toc-modified-id="Turn-the-aerosol-dictionary-into-a-named_tuple-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Turn the aerosol dictionary into a named_tuple</a></span></li><li><span><a href="#Define-a-function-to-calculate-a-and-b-for-any-aerosol" data-toc-modified-id="Define-a-function-to-calculate-a-and-b-for-any-aerosol-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Define a function to calculate a and b for any aerosol</a></span></li><li><span><a href="#Call-it-on-this-aerosol-and-calculate-SScrit" data-toc-modified-id="Call-it-on-this-aerosol-and-calculate-SScrit-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Call it on this aerosol and calculate SScrit</a></span></li></ul></li><li><span><a href="#save-these-aerosol-data-for-future-calculations" data-toc-modified-id="save-these-aerosol-data-for-future-calculations-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>save these aerosol data for future calculations</a></span></li><li><span><a href="#Read-the-json-files-back-in-and-work-with-them" data-toc-modified-id="Read-the-json-files-back-in-and-work-with-them-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Read the json files back in and work with them</a></span></li></ul></div>

# # Given the critical supersaturation from the kohler notes:
# 
# 
# $SS=S^* - 1= \left ( \frac{4 a^3}{27b} \right )^{1/2}$
# 
# show that this implies, for $(NH_4)_2 SO_4$, density $\rho_{aer}$ = 1775
# ${kg}\,{m^{-3}}$ , van hoft i=3, that:
#   
# $S^* -1 \approx 1.54 \times 10^{-12}~ m_{aer}^{-0.5}$
# 
# where $m_{aer}$ is the ammonium sulphate aerosol mass in kg.
# 
# Note that this is why a cloud chamber can get the aerosol mass distribution from a series of
# saturation and light scattering measurements as smaller and smaller aerosols are pushed over
# their critical supersaturation and activated.
# 
# 

# ## specify the aerosol properties

# In[73]:


aero_amonium_sulphate=dict([   ('Ms', 132),
   ('Mw', 18.0),
   ('Sigma', 0.075),
   ('vanHoff', 3.0),
   ('comments','ammonum sulfate (NH4)2SO4')])

aero_sodium_chloride=dict([   ('Ms', 58),
   ('Mw', 18.0),
   ('Sigma', 0.075),
   ('vanHoff', 2.0),
   ('comments','sodium chloride NaCl')])

#
# load pprint to print nested data structures
#
import pprint
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(aero_amonium_sulphate)


# ## Turn the aerosol dictionary into a named_tuple
# 
# I would rather write the vanHoff factor as:
# 
#      aerosol.vanHoff
#      
# instead of:
# 
#      aerosol['vanHoff']
#      
# in my formulae.
# 
# The utility function [a405.utils.helper_funs.make_tuple](https://github.com/phaustin/atsc405_2018/blob/master/a405/utils/helper_funs.py#L26-L47) turns a dictionary into a [namedtuple](https://docs.python.org/3/library/collections.html#namedtuple-factory-function-for-tuples-with-named-fields)

# In[75]:


from a405.utils.helper_funs import make_tuple
aero = make_tuple(aero_amonium_sulphate)
print(aero)


# ## Define a function to calculate a and b for any aerosol
# 
# Follow the hint and copy it from [day24](https://clouds.eos.ubc.ca/~phil/courses/atsc405/html/kohler.html)

# In[76]:


from a405.thermo.constants import constants as c
import numpy as np

def calc_kohler_coeffs(aero, Temp):
    a=(2.*aero.Sigma)/(c.Rv*Temp*c.rhol)  #curvature term
    b=(aero.vanHoff*aero.Mw)/((4./3.)*np.pi*c.rhol*aero.Ms)  #Raoult term without aerosol mass
    return a,b


# ## Call it on this aerosol and calculate SScrit

# In[77]:


a, b = calc_kohler_coeffs(aero,280.)
coeff = (4.*a**3./(27.*b))**0.5
print(f'The coefficient is {coeff:6.2e} kg**(0.5)')


# # save these aerosol data for future calculations
# 
# How do I save the two aerosol dictionaries so I can load them from other programs?
# 
# I have a data directory that I use for this kind of information called a405/data.  I'll write out the two dictionaries aero_amonium_sulpate and aero_sodium_chloride as [json](https://en.wikipedia.org/wiki/JSON) files.  I'll use the
# [pathlib module](https://docs.python.org/3/library/pathlib.html) to construct the path to the individual json files
# starting from the location of the a405.data folder.

# In[87]:


import a405.data
from pathlib import Path
import json
#
# locate the a405/data folder on this computer
#
dirpath=a405.data.__path__

the_path = Path(dirpath[0])
#
#  add the indivdual file names
#
nacl_file = the_path / Path('aero_nacl.json')
amsulphate_file = the_path / Path('ammonium_sulphate.json')
file_paths=[amsulphate_file,nacl_file]
dicts = [aero_amonium_sulphate, aero_sodium_chloride]
for the_path,the_dict in zip(file_paths,dicts):
    with open(the_path,'w') as f:
        json.dump(the_dict,f,indent=4)


# Check the a405/data folder to make sure these files have been written.

# # Read the json files back in and work with them
# 
# How do I find the files from another program/notebook once I've written them?   There is a module 
# called [importlib_resources](http://importlib-resources.readthedocs.io/en/latest/) that will be part of the next python release (3.7).  For now you can install with pip:
#          
#           pip install importlib_resources
#           
# As long as you know that the data folder is in the a405 package you can get its location on your computer like this:

# In[97]:


import importlib_resources as ir 
contents=list(ir.contents('a405.data'))
print(f'contents of the data folder={contents}\n')

with ir.open_text('a405.data','ammonium_sulphate.json') as f:
    aero_dict=json.load(f)

print('here is the sulfate aerosol info:')
pp.pprint(aero_dict)


# In[96]:


list(ir.contents('a405.data'))
dir(ir)

