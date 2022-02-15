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

# Rootfinding - Solution

+++ {"toc": true}

## Table of Contents

<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#The-cell-below-shows-how-to-use-a-rootfinder-to-solve-cos(x)-=-0.75" data-toc-modified-id="The-cell-below-shows-how-to-use-a-rootfinder-to-solve-cos(x)-=-0.75-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>The cell below shows how to use a rootfinder to solve cos(x) = 0.75</a></span></li><li><span><a href="#Problem-for-Friday:--9am" data-toc-modified-id="Problem-for-Friday:--9am-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Problem for Friday:  9am</a></span><ul class="toc-item"><li><span><a href="#Solution" data-toc-modified-id="Solution-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Solution</a></span></li></ul></li><li><span><a href="#Bracket-finding" data-toc-modified-id="Bracket-finding-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Bracket finding</a></span><ul class="toc-item"><li><span><a href="#example----find-a-bracket-for-sin(x)=0-near-x=12-radians-~-700-degrees" data-toc-modified-id="example----find-a-bracket-for-sin(x)=0-near-x=12-radians-~-700-degrees-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>example -- find a bracket for sin(x)=0 near x=12 radians ~ 700 degrees</a></span></li></ul></li><li><span><a href="#now-use-the-fzero-wrapper-to-find-the-root-of-sin(x)=0--(720-degrees)" data-toc-modified-id="now-use-the-fzero-wrapper-to-find-the-root-of-sin(x)=0--(720-degrees)-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>now use the fzero wrapper to find the root of sin(x)=0  (720 degrees)</a></span></li><li><span><a href="#Redo-theta-example-with-find_interval" data-toc-modified-id="Redo-theta-example-with-find_interval-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Redo theta example with find_interval</a></span></li></ul></div>

+++

The cell below shows how to use a rootfinder to solve cos(x) = 0.75

```{code-cell} ipython3
%matplotlib inline
from scipy import optimize
from matplotlib import pyplot as plt
import numpy as np
plt.style.use('ggplot')


xvals=np.linspace(0,10.)
fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(xvals,np.cos(xvals))
straight_line=np.ones_like(xvals)
ax.plot(xvals,straight_line*0.75,'b-')

def root_function(x):
    """Function we want to find the root of
       input: x value
       output: y value -- should be zero when x is a root
    """
    return np.cos(x) - 0.75

root1 = optimize.zeros.brentq(root_function,0,2)
root2 = optimize.zeros.brentq(root_function,4,6)
root3 = optimize.zeros.brentq(root_function,6,8)
xvals=np.array([root1,root2,root3])
yvals=np.cos(xvals)
ax.plot(xvals,yvals,'ro',markersize=20)

print(optimize.zeros.__file__)
```

## Problem for Friday:  9am

Write a function:
    
    def temp_from_theta(theta,press)
    
that uses brentq to find the temperature for a theta of 280 K and pressures of press=[1.e5, 7.e4, 4.e4] Pascals

Email a working notebook as an attachment, or send a url to a dropbox-like repository

+++

### Solution

```{code-cell} ipython3
def root_function(Tguess,theta,press):
    """Function we want to find the root of
       input: Tguess (K), target theta (K), press (Pa)
       output: difference between guess and target -- should be zero when x is a root
    """
    p0=1.e5  #Pa
    Rd=287. #J/kg/K
    cp=1004.  #J/kg/K     
    theta_guess=Tguess*(p0/press)**(Rd/cp)
    return theta - theta_guess
```

```{code-cell} ipython3
def temp_from_theta(theta,press):
    """
       input: theta (K), press (Pa)
       output: temp (K) found by rootfinder
    """     
    left=10 #K
    right=1000 #K
    temp = optimize.zeros.brentq(root_function,left,right,args=(theta,press))
    return temp

for press in [1.e5,7.e4,4.e4]:
    print('Temp {:5.2f} (K) at pressure of {:5.2f} kPa'.format(temp_from_theta(280.,press),press*1e-2))
```

### Bracket finding

I've written a couple of convenience functions called rootfinder.find_interval and
rootfinder.fzero to make rootfinding a little easier.   The new module is 
[rootfinder.py](https://github.com/phaustin/atsc405_2018/blob/master/a405/thermo/rootfinder.py)

```{code-cell} ipython3
from a405.thermo import rootfinder as rf
```

```{code-cell} ipython3
help(rf.find_interval)
```

```{code-cell} ipython3
help(rf.fzero)
```

### example -- find a bracket for sin(x)=0 near x=12 radians ~ 700 degrees

```{code-cell} ipython3
brackets=rf.find_interval(np.sin,12)
brackets
```

### now use the fzero wrapper to find the root of sin(x)=0  (720 degrees)

```{code-cell} ipython3
print(rf.fzero(np.sin,brackets)*180/np.pi)
```

### Redo theta example with find_interval

```{code-cell} ipython3
import a405.thermo.rootfinder as rf

def theta_zero(Tguess,theta,press):
    """Function we want to find the root of
       input: Tguess (K), target theta (K), press (Pa)
       output: difference between guess and target -- should be zero when x is a root
    """
    p0=1.e5  #Pa
    Rd=287. #J/kg/K
    cp=1004.  #J/kg/K     
    theta_guess=Tguess*(p0/press)**(Rd/cp)
    return theta - theta_guess
```

```{code-cell} ipython3
def temp_from_theta(theta,press):
    """
       input: theta (K), press (Pa)
       output: temp (K) found by rootfinder
    """     
    #
    #  use theta as guess for bracket and pass theta,press to theta_zero
    #
    brackets=rf.find_interval(theta_zero,theta,theta,press)
    the_temp = rf.fzero(theta_zero,brackets,theta,press)
    return the_temp

for press in [1.e5,7.e4,4.e4]:
    print('Temp {:5.2f} (K) at pressure of {:5.2f} kPa'.format(temp_from_theta(280.,press),press*1e-2))
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```
