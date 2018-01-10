
# coding: utf-8

# ## The cell below shows how to use a rootfinder to solve cos(x) = 0.75

# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
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
    """
    This is the function we want to find the root of
    
    Parameters
    ----------
    
       x: float
          dependent variable
          
    Returns
    -------
    
    Result of function, should be zero when x is a root

    """
    return np.cos(x) - 0.75 

#
# Solve cos(x) - 0.75 = 0 for three different x intervals
#
#  0-2, 4-6, and 6-8
#
root1 = optimize.zeros.brentq(root_function,0,2)
root2 = optimize.zeros.brentq(root_function,4,6)
root3 = optimize.zeros.brentq(root_function,6,8)
xvals=np.array([root1,root2,root3])
yvals=np.cos(xvals)
ax.plot(xvals,yvals,'ro',markersize=20)

print(optimize.zeros.__file__)

## Problem for Friday:  9am

Write a function:
    
    def temp_from_theta(theta,press)
    
that uses brentq to find the temperature for a $\theta$ of 280 K and pressures of press=[1.e5, 7.e4, 4.e4] Pa

Share the google drive url to your notebook.