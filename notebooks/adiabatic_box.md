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
<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Notebook-practice" data-toc-modified-id="Notebook-practice-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Notebook practice</a></span></li><li><span><a href="#Adiabatic-box" data-toc-modified-id="Adiabatic-box-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Adiabatic box</a></span></li></ul></div>

+++

# Notebook practice

    
Write a function called "eqstate" that calculates the density of dry air.  Use it to find the dry air density
at a pressure of 80 kPa and temperatures of temp=[270, 280, 290] K

```{code-cell} ipython3
def eqstate(Temp,Press):
    dens = xxx
    return dens
```

# Adiabatic box

+++

The figure below represents an insulated box with two
compartments A and B, each containing dry air. They are separated by
an insulating and perfectly flexible wall, so that the pressure is
equal on the two sides. Initially each compartment measure one and
the gas is at a pressure of 100 kPa and a temperature of 273 K. Heat
is then supplied to the gas in box A using a resistor, until the
pressure in both compartments rises to 300 kPa. Calculate:

1. The final temperature  $T_B$ (K) in box B

2. The time integrated work ($J\,kg^{-1}$) performed on the air in B by the
   membrane.

3. The final temperature $T_A$ (K) in box A

4. The time-integrated heating rate ($J kg^{-1}$) of the gas in box A.

Hint: use the fact that entropy (and therefore $\theta$ ) is conserved for an adiabatic process and therefore can't change in compartment B, and
that the total volume of the combined compartments A and B has to stay constant at 2 $m^3$.

```{code-cell} ipython3
from IPython.display import Image
Image(filename="images/insulated_box.png")
```
