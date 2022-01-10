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

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
```

# Notebook practice

    
Write a function called "eqstate" that calculates the density of dry air.  Use it to find the dry air density
at a pressure of 80 kPa and temperatures of temp=[270, 280, 290] K

+++

*Solve Lohmann 2.4 for $\rho$ to obtain:*

$$
\rho_d = \frac{p}{R_dT}\tag{Lohmann 2.4*}
$$

*with $p$ in Pascals, $T$ in K, and $R_d = 287 J kg^{−1} K^{−1}$. Write the function in Python:*

```{code-cell} ipython3
def eqstate(Temp, Press):
    Rd = 287  # J/kgK
    Press *= 1e3  # convert to Pa
    dens = Press / (Rd * Temp)
    return dens
```

*Use the function to find the dry air density at a pressure of 80 kPa and temperatures of temp=[270, 280, 290] K*

```{code-cell} ipython3
p = 80  # kPa

print("Andrew's Soln\n" + "*" * 30)
for temp in [270, 280, 290]:
    print(
        f"""    P:    {p} kPa
    T:    {temp} K
    dens: {round(eqstate(temp, p),3)} kg/m^3
    """
    )
```

# Adiabatic box

+++

The figure below represents an insulated box with two
compartments A and B, each containing dry air. They are separated by
an insulating and perfectly flexible wall, so that the pressure is
equal on the two sides. Initially each compartment measures one cubic meter and
the gas is at a pressure of 100 kPa with a temperature of 273 K. Heat
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

***1) Final Temperature in Box B***

*Lohmann 2.37 equation for the quantity $\theta$, which is conserved in box B:*

$$
\theta = T\left(\frac{p_0}{p}\right)^\kappa\tag{Lohmann 2.37}
$$

*Initially, $p_0 = p$, so*

$$
\theta = T_0
$$

*$\kappa$ for dry air is equal to 0.286 (Lohmann pg 38). Initially, $\theta = T$, and we can solve for the final temperature $T_B$ in box B by rearranging and plugging in the final pressure $p=300$ kPa*

$$
T_B = T_0\left(\frac{p}{p_0}\right)^\kappa
$$

```{code-cell} ipython3
# initial state of each box (they're the same)
P_init = 100  # kPa
T_init = 273  # K
V_init = 1  # m3
kappa = 0.286

# final temperature of box B
P_final = 300  # kPa
TB = T_init * (P_final / P_init) ** kappa
print(f"TB = {round(TB,2)} K")
```

***2) Time-Integrated Work on Box B***

*Looking at figure 2.5, we see the work done on box B can be represented as the area underneath an adiabat connnecting the initial state $(P_i,V_i)$ to the final state $(P_f, V_f)$.* 

$$
W = \int p\cdot dv
$$

*Integrate and sketch the curve with matplotlib*

```{code-cell} ipython3
res = 1e-3  # numerical integration resolution
PB = np.arange(P_init, P_final, res)  # create an array of pressures
TB = T_init * (PB / P_init) ** kappa  # calculate corresponding temperatures
VB = eqstate(T_init, P_init) / eqstate(TB, PB)  # get corresponding volumes

# integrate to get total work in the box
ttl_work = sum(PB[1:] * np.diff(VB))

# convert to J/Kg as specified - divide by initial density
work_B = ttl_work / eqstate(T_init, P_init)
```

```{code-cell} ipython3
# plot the result
fig, ax = plt.subplots()
ax.plot(VB, PB, "k")
ax.scatter([VB[0], VB[-1]], [PB[0], PB[-1]], color="k")
ax.set_xlim(VB[-1] - 0.1, VB[0] + 0.1)
ax.set_ylim(0, max(PB) + 1e5)
ax.fill_between(VB, PB, color="k", alpha=0.2)
ax.annotate(f"W = {round(work_B)} J/kg", (0.6, 50000))
ax.set_xlabel("V (m$^3$)")
ax.set_ylabel("P (Pa)")
ax.set_title("Work Done on Box B");
```

***3) Final Temperature in Box A***

*Since we know the final volume of box B, along with the total volume of both boxes, and the pressure, we can get the final temperature from the equation of state:*

$$
PV = nKT
$$

$$
\frac{P_1V_1}{T_1} = \frac{P_2V_2}{T_2}
$$

$$
T_2 = \frac{P_2V_2}{P_1V_1}\cdot T_1
$$

```{code-cell} ipython3
VA = 2 - VB[-1] # final volume of box A (m3)
T2 = (P_final * VA * T_init) / (P_init * V_init)
print(f"T2 = {round(T2,2)} K")
```

**4) Total Heat Added to Box A**

The total work done on the system is the energy reqired to heat the air in box A + the energy required to compress box B. 

$$
q_T = c_v\Delta T_A + W_B
$$

with $c_v = 718 J/kg \cdot K$ (Lohmann 2.13)

```{code-cell} ipython3
cv = 718  # J/kgK
qT = cv * (T2 - T_init) + work_B
print(f"Total heat added to A: {round(qT,2)} J")
```
