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

# Marshall-Palmer Distribution

<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Marshall-Palmer" data-toc-modified-id="Marshall-Palmer-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Marshall Palmer</a></span></li><li><span><a href="#Rain-rate" data-toc-modified-id="Rain-rate-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Rain rate</a></span></li></ul></div>

```{code-cell} ipython3
### Rewrite the MP distribtuion as a function
```

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
%matplotlib inline
plt.style.use('ggplot')

def marshallpalmer(R):
    """
    marshall palmer size distribution
    given rainrate R in mm/hr, return
    n(D), the number concentration of drops with
    diameter D

    Parameters
    ----------
    R: float
        rainrate (mm/hr)

    Returns
    -------

    d: vector (float)
      drop diameters (cm)

    n: vector (float)
     the number distribution n(d) #m^{-3} mm^{-1}

    """
    D=np.arange(0,8,0.01)
    Dmm=D
    Dcm=D*0.1
    N0=0.08*1.e6*1.e-1 #m**{-3} mm^{-1}
    theLambda=41*R**(-0.21)
    n=N0*np.exp(-theLambda*Dcm)
    return Dcm,n

curve_dict={}
Rvals = [1,5,25]
for R in Rvals:
    diam,ndist = marshallpalmer(R)
    curve_dict[R] = ndist
fig, ax = plt.subplots(1,1,figsize=(10,8))
for R in Rvals:
    ax.semilogy(diam,curve_dict[R],label='{} mm/hr'.format(R))
ax.set_xlabel('Drop diameter (mm)')
ax.set_ylabel('$n(D) m^{-3} mm^{-1}$')
ax.set_title('Marshall Palmer distribution for three rain rates')
out=ax.legend(loc='best')
```

**Confirm that the mean diameter of the Marshall Palmer distribution is $1/\Lambda$**

Number of droplets between diameter $D$ and $D + dD$
$$
N(D) = N_0e^{-\Lambda D}\tag{AT 4.31}
$$

Integrate and divide by the number of drops to get the mean:

$$
\overline{N(D)} = \frac{1}{N_0}\int^\infty_0 N_0e^{-\Lambda D}
$$

$$
= -\frac{1}{\Lambda}e^{-\Lambda D}\bigg\rvert_{D=0}^\infty
$$

$$
= -\frac{1}{\Lambda}\left(e^{\infty} - e^{0}\right)
$$

$$
\boxed{\overline{N(D)} = \frac{1}{\Lambda}}
$$

+++

## Rain rate


- Calculate the precipitation flux (mm/hr) by integrating the total volume in the Marshall Palmer size distribution and with the fall speed of Villermaux and Bossa (2009): $w = - \sqrt{\rho_l/\rho_{air} * g *D}$ where
$D$ is the drop diameter and $\rho_l,\rho_{air}$ are the liquid and air densities.  Show
that you get about $R=15 mm/hour$ back from this calculation if you use the $\Lambda$ that's appropriate for $R=15 mm/hour$

+++

$$
-\Lambda (R) = 41R^{-0.21}\tag{AT 4.32}
$$

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np
```

```{code-cell} ipython3
def Lambda_from_R(R):
    """
    calculates Marshall palmer parameter Lambda using Thompkins 4.32
    
    inputs
    -------
    R: (float) rainfall rate (mm/hour)
    
    returns
    -------
    Lambda: (float) marshall-palmer parameter
    """
    return 41 * R ** -0.21

def N_from_Lambda(lam, D):
    """
    calculates the number of droplets per 
    """
```

```{code-cell} ipython3
RR = Lambda_from_R(15)
RR
```

```{code-cell} ipython3
N0 = 0.08 # cm^-4
rhoa = 1. # density of air
w = - np.sqrt(c.rhol/rhoa * c.g0 * D)
```

```{code-cell} ipython3
diam = np.linspace(0,1e-2,100) # droplet diameter (m)
vol = 4 / 3 * np.pi * (diam / 2) ** 3 # volume of each droplet (m3)
vol
```
