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
<div class="toc"><ul class="toc-item"><li><span><a href="#Parcel-model-with-30-aerosol-masses,-lognormal-distribution" data-toc-modified-id="Parcel-model-with-30-aerosol-masses,-lognormal-distribution-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Parcel model with 30 aerosol masses, lognormal distribution</a></span><ul class="toc-item"><li><span><a href="#Read-in-the-json-file-and-set-the-koehler-function-for-this-aerosol" data-toc-modified-id="Read-in-the-json-file-and-set-the-koehler-function-for-this-aerosol-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Read in the json file and set the koehler function for this aerosol</a></span></li><li><span><a href="#initialize-the-lognormal-mass-and-number-distributions-for-30-bins" data-toc-modified-id="initialize-the-lognormal-mass-and-number-distributions-for-30-bins-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>initialize the lognormal mass and number distributions for 30 bins</a></span></li><li><span><a href="#find-the-equilibrium-radius-for-each-bin-at-saturation-Sinit" data-toc-modified-id="find-the-equilibrium-radius-for-each-bin-at-saturation-Sinit-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>find the equilibrium radius for each bin at saturation Sinit</a></span></li><li><span><a href="#now-add-the-intial-conditions-to-the-cloud_vars-dictionary-and-make-it-a-namedtuple" data-toc-modified-id="now-add-the-intial-conditions-to-the-cloud_vars-dictionary-and-make-it-a-namedtuple-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>now add the intial conditions to the cloud_vars dictionary and make it a namedtuple</a></span></li><li><span><a href="#use-odeint-to-integrate-the-variable-in-var_vec-from-tinit-to-tfin-with-outputs-every-dt-seconds" data-toc-modified-id="use-odeint-to-integrate-the-variable-in-var_vec-from-tinit-to-tfin-with-outputs-every-dt-seconds-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>use odeint to integrate the variable in var_vec from tinit to tfin with outputs every dt seconds</a></span></li><li><span><a href="#create-a-dataframe-with-33-columns-to-hold-the-data" data-toc-modified-id="create-a-dataframe-with-33-columns-to-hold-the-data-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>create a dataframe with 33 columns to hold the data</a></span></li><li><span><a href="#store-the-dataframe-in-an-csv-file,-including-a-copy-of-the-input-dictionary-for-future-reference" data-toc-modified-id="store-the-dataframe-in-an-csv-file,-including-a-copy-of-the-input-dictionary-for-future-reference-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>store the dataframe in an csv file, including a copy of the input dictionary for future reference</a></span></li><li><span><a href="#Matt's-solution" data-toc-modified-id="Matt's-solution-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Matt's solution</a></span><ul class="toc-item"><li><span><a href="#do-test-at-z-=-1010-meters" data-toc-modified-id="do-test-at-z-=-1010-meters-1.8.1"><span class="toc-item-num">1.8.1&nbsp;&nbsp;</span>do test at z = 1010 meters</a></span></li><li><span><a href="#find-the-differentials" data-toc-modified-id="find-the-differentials-1.8.2"><span class="toc-item-num">1.8.2&nbsp;&nbsp;</span>find the differentials</a></span></li><li><span><a href="#write-everything-out-as-a-tuple" data-toc-modified-id="write-everything-out-as-a-tuple-1.8.3"><span class="toc-item-num">1.8.3&nbsp;&nbsp;</span>write everything out as a tuple</a></span></li><li><span><a href="#return-$\Delta-SS$-given-the-tuple-t" data-toc-modified-id="return-$\Delta-SS$-given-the-tuple-t-1.8.4"><span class="toc-item-num">1.8.4&nbsp;&nbsp;</span>return $\Delta SS$ given the tuple t</a></span></li><li><span><a href="#Check-against-the-output-$\Delta-SS$" data-toc-modified-id="Check-against-the-output-$\Delta-SS$-1.8.5"><span class="toc-item-num">1.8.5&nbsp;&nbsp;</span>Check against the output $\Delta SS$</a></span></li></ul></li></ul></li></ul></div>

+++

# Parcel model with 30 aerosol masses, lognormal distribution

```{code-cell} ipython3
import json
import a405.utils
from pathlib import Path
import numpy as np
from a405.dropgrow.aerolib import lognormal,create_koehler
from a405.utils.helper_funs import make_tuple, find_centers
from collections import OrderedDict as od
from a405.thermo.thermlib import find_esat, find_lv
from a405.thermo.rootfinder import find_interval, fzero
from a405.dropgrow.drop_grow import find_diff, rlcalc, find_derivs, Scalc, rlderiv
from a405.thermo.constants import constants as c
from scipy.integrate import odeint
import pandas as pd
from matplotlib import pyplot as plt
import datetime
import importlib_resources as ir 
import pprint
pp = pprint.PrettyPrinter(indent=4)
```

```{code-cell} ipython3
%matplotlib inline
```

## Read in the json file and set the koehler function for this aerosol

```{code-cell} ipython3
with ir.open_text('a405.data','dropgrow.json') as f:
    input_dict=json.load(f)
pp.pprint(input_dict)

aero=make_tuple(input_dict['aerosol'])
parcel=make_tuple(input_dict['initial_conditions'])

koehler_fun = create_koehler(aero,parcel)
    
```

## initialize the lognormal mass and number distributions for 30 bins

```{code-cell} ipython3
#
#set the edges of the mass bins
#31 edges means we have 30 droplet bins
#
numrads = 30
mass_vals = np.linspace(-20,-16,numrads+1) 
mass_vals = 10**mass_vals  #aerosol mass in kg
mu=input_dict['aerosol']['themean']
sigma = input_dict['aerosol']['sd']
totmass = input_dict['aerosol']['totmass']
mdist = totmass*lognormal(mass_vals,np.log(mu),np.log(sigma))
mdist = find_centers(mdist)*np.diff(mass_vals)  #kg/m^3 of aerosol in each bin
center_mass = find_centers(mass_vals)
ndist = mdist/center_mass  #number/m^3 of aerosol in each bin
#save these in an ordered dictionary to pass to functions
cloud_vars = od()
cloud_vars['mdist'] = mdist
cloud_vars['ndist'] = ndist
cloud_vars['center_mass'] = center_mass
cloud_vars['koehler_fun'] = koehler_fun
```

## find the equilibrium radius for each bin at saturation Sinit

```{code-cell} ipython3
S_target = parcel.Sinit
logr_start = np.log(0.1e-6)

initial_radius = []
dry_radius = []
for mass in center_mass:
    brackets = np.array(find_interval(find_diff,logr_start,S_target,mass,koehler_fun))
    left_bracket, right_bracket = np.exp(brackets)*1.e6  #get brackets in microns for printing
    equil_rad = np.exp(fzero(find_diff,brackets,S_target,mass,koehler_fun))

    initial_radius.append(equil_rad)
    dry_rad = (mass/(4./3.*np.pi*aero.rhoaero))**(1./3.)
    dry_radius.append(dry_rad)

    print('mass = {mass:6.3g} kg'.format_map(locals()))
    print('equlibrium radius at S={} is {:5.3f} microns\n'.format(S_target,equil_rad*1.e6))
```

## now add the intial conditions to the cloud_vars dictionary and make it a namedtuple

the vector var_vec holds 30 droplet radii plus three extra variables at the
end of the vector: the temperature, pressure and height.

```{code-cell} ipython3
nvars =4
```

```{code-cell} ipython3
cloud_vars['initial_radiius'] = initial_radius
cloud_vars['dry_radius'] = dry_radius
cloud_vars['masses'] = center_mass
numrads = len(initial_radius)
var_vec = np.empty(numrads +  nvars)
for i in range(numrads):
    var_vec[i] = initial_radius[i]

#
# temp, press and height go at the end of the vector
#
var_vec[-4] = parcel.Sinit-1 # initialize SS
var_vec[-3] = parcel.Tinit
var_vec[-2] = parcel.Pinit
var_vec[-1] = parcel.Zinit

cloud_tup = make_tuple(cloud_vars)
#calculate the total water (kg/kg)
rl=rlcalc(var_vec,cloud_tup,nvars);
e=parcel.Sinit*find_esat(parcel.Tinit);
rv=c.eps*e/(parcel.Pinit - e)
#save total water
cloud_vars['rt'] = rv + rl
cloud_vars['wvel'] = parcel.wvel
cloud_vars['wvel'] = 7.
#
# pass this to the find_derivs function
#
cloud_tup= make_tuple(cloud_vars)

```

## use odeint to integrate the variable in var_vec from tinit to tfin with outputs every dt seconds

```{code-cell} ipython3
def find_derivs_SS(var_vec,the_time,cloud_tup):
    """
    calcuate derivatives of var_vec 

    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    the_time: float
       timestep 

    cloud_tup: namedtuple
           tuple of necessary coefficients
    

    Returns
    -------

    deriv_vec: vector(float)
         derivatives of each of var_vec
    
    """
    #print('inside: ',var_vec)
    SS,temp,press,height = var_vec[-4:]
    numrads = len(var_vec) - nvars
    dry_radius = cloud_tup.dry_radius
    rho=press/(c.Rd*temp)
    #
    # find the evironmental S by water balance
    #
    ## S=Scalc(var_vec,cloud_tup)
    S = SS + 1
    deriv_vec=np.zeros_like(var_vec)
    #dropgrow notes equaton 18 (W&H p. 170)
    for i in range(numrads):
        m=cloud_tup.masses[i]
        if var_vec[i] < dry_radius[i]:
            var_vec[i] = dry_radius[i]
        Seq=cloud_tup.koehler_fun(var_vec[i],m)  
        rhovr=(Seq*find_esat(temp))/(c.Rv*temp)
        rhovinf=S*find_esat(temp)/(c.Rv*temp)
        #day 25 drop_grow.pdf eqn. 18
        deriv_vec[i]=(c.D/(var_vec[i]*c.rhol))*(rhovinf - rhovr)

    #
    # moist adiabat day 25 equation 21a
    #
    deriv_vec[-3]=find_lv(temp)/c.cpd*rlderiv(var_vec,deriv_vec,cloud_tup,nvars) - c.g0/c.cpd*cloud_tup.wvel
    #
    # hydrostatic balance  dp/dt = -rho g dz/dt
    #
    deriv_vec[-2]= -1.*rho*c.g0*cloud_tup.wvel
    #
    # how far up have we traveled?
    #
    deriv_vec[-1] = cloud_tup.wvel
    
        #
    # dSS/dt using equil_super notes
    #
    es = find_esat(temp)
    dp_dt = deriv_vec[-2]
    dT_dt = deriv_vec[-3]
    des_dT = c.eps * find_lv(temp) * es / (c.Rd * (temp ** 2))
    dSS_bracket = (-c.eps * es / (press ** 2)) * dp_dt + (c.eps / press) * des_dT * dT_dt
    dSS_divide = c.eps * es / press
    drl_dt = rlderiv(var_vec,deriv_vec,cloud_tup, nvars)
    deriv_vec[-4]= (drl_dt - (1+SS) * dSS_bracket ) / dSS_divide
    
    return deriv_vec
```

```{code-cell} ipython3
var_out = []
time_out =[]

tinit=input_dict['integration']['dt']
dt = input_dict['integration']['dt']
tfin = input_dict['integration']['tend']

t = np.arange(0,tfin,dt)
sol = odeint(find_derivs_SS,var_vec, t, args=(cloud_tup,))
```

## create a dataframe with 33 columns to hold the data

```{code-cell} ipython3
colnames = ["r{}".format(item) for item in range(30)]
colnames.extend(['SS','temp','press','z'])
df_output = pd.DataFrame.from_records(sol,columns = colnames)
```

## store the dataframe in an csv file, including a copy of the input dictionary for future reference

```{code-cell} ipython3
if input_dict['dump_output']:
    outfile_name = f'{input_dict["output_file"]}.csv'
    with open(outfile_name,'w') as store:
       df_output.to_csv(store)

    metadata_name = f'{input_dict["output_file"]}.json'
    date=datetime.datetime.now().strftime('%Y-%M-%d')
    with open(metadata_name,'w') as meta:
        history ="file produced by drop_grow.py on {}".format(date)
        print('history: ',history)
        input_dict['history']=history
        json.dump(input_dict,meta,indent=4)
```

```{code-cell} ipython3
print(f"{df_output['SS'][0]+1} {parcel.Sinit}")
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1,figsize=[10,8])
for i in colnames[:-nvars]:
    ax.plot(df_output[i]*1.e6,df_output['z'],label=i)
out=ax.set(ylim=[1000,1040],xlim=[0,6],
       xlabel='radii (microns)',ylabel='height (m)',
              title='radii vs. height in a {} m/s updraft'.format(cloud_tup.wvel))
```

cool

```{code-cell} ipython3
#Svals = []
#for index,row in df_output.iterrows():
#    var_vec = row.values
#    Svals.append(Scalc(var_vec,cloud_tup))
fig,ax = plt.subplots(1,1,figsize=[10,8])
ax.plot(df_output['SS']+1,df_output['z'])
out=ax.set(ylim=[1000,1150],title='Saturation in a {} m/s updraft'.format(cloud_tup.wvel))
```

```{code-cell} ipython3

```
