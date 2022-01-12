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

# Two Compartment

+++

Consider two compartments, each with a volume of $1 (m^3)$, separated by a rigid, perfectly insulating membrane.
Initially $(T_A) = 300 K$ and $(p_A) = (10^5) Pa$, and $(T_B) = 100 K$ and $(p_B) = (10^3) Pa$. Suppose the membrane breaks. Find the final temperature and pressure in the $2 (m^3)$ compartment. Does the total entropy change? By how much?

```{code-cell} ipython3
import numpy as np

# constants
Rd = 287  # J/kgK
cv = 718  # J/kgK
cp = 1005  # J/kgK
kappa = 0.286  # L page 38
p0 = 1e5  # Pa (reference pressure)
```

*First calculate the total internal energy as the sum of $U$ in each compartment:*

$$
du = c_vdT\tag{L2.14}
$$

*Integrate from $(0,0)$ to $(T, u)$, then multiply by the mass (density $\times$ volume). This is the energy required to "assemble" each compartment from nothing.*

$$
\int_0^udu = c_v\int_0^TdT
$$

$$
u = c_vT \text{    (J/kg)}
$$

$$
U = c_vT \rho V\text{    (J)}\tag{1}
$$

*Total internal energy after mixing is the sum of $U$ in each compartment (energy is conserved):*

$$
U_T = U_A + U_B
$$

*Invert (1) to back out the final temperature*

$$
T_f = \frac{U_T}{c_v \rho_F V_f}
$$

*Final density $\rho_f$ is just the average of the initial densities (since the initial volumes are the same). $V_f = 2 m^3$*

$$
\rho_f = \frac{\rho_A + \rho_B}{2}
$$

*Now get the final pressure from the equation of state:*

$$
p = \rho_dR_d T \tag{L2.4}
$$

```{code-cell} ipython3
# initial state of compartment A
P_A = 1e5  # Pa
T_A = 300  # K
rho_A = P_A / (Rd * T_A)  # kg/m3
U_A = cv * T_A * rho_A  # J

# intiial state of compartment B
P_B = 1e3  # Pa
T_B = 100  # K
rho_B = P_B / (Rd * T_B)  # kg/m3
U_B = cv * T_B * rho_B  # J

# final state
U_T = U_A + U_B # total internal energy (J)
rho_F = (rho_A + rho_B) / 2  # final density (kg/m3)
T_F = U_T / (cv * rho_F * 2 )  # final temperature (K)
P_F = (P_A / T_A + P_B / T_B) * T_F / 2

print(f"Final temperature: {T_F} K\nFinal Pressure:    {P_F} Pa")
```

Does the total entropy change?

*Yes. This process is obviously irreversible (I cannot think of an easy way to put all the air back in compartment A)*

By how much?

*Starting at the entropy notes eq (6) (equiv to Lohmann 2.40):*

$$
s_1 - s_0 = c_p log\frac{\theta_1}{\theta_0}\tag{PA6}
$$

*Calculate the entropy change in each compartment separately using the initial and final potential temperatures (both compartments have the same final state). Set $p_0 = 10^5 Pa$.*

$$
\theta = T\left(\frac{p_0}{p}\right)^\kappa\tag{L2.37}
$$

$$
\theta_{Ai} = T_{A}
$$

$$
\theta_{Bi} = T_B\left(\frac{p_0}{p_B}\right)^\kappa
$$

$$
\theta_f = T_f\left(\frac{p_0}{p_f}\right)^\kappa
$$

$$
\Delta s_{total} = c_p log\frac{\theta_f}{\theta_{Ai}} + c_p log\frac{\theta_f}{\theta_{Bi}}
$$

$$
\Delta s_{total} = c_p log\left(\frac{\theta_f^2}{\theta_{Ai}\theta_{Bi}}\right)\tag{2}
$$

```{code-cell} ipython3
# get potential temperature in each state
theta_A = T_A * (p0 / P_A) ** kappa
theta_B = T_B * (p0 / P_B) ** kappa
theta_F = T_F * (p0 / P_F) ** kappa

# plug into (2)
del_s = cp * np.log(theta_F ** 2 / (theta_A * theta_B))

print(f"Total entropy change: {del_s} J/K")
```
