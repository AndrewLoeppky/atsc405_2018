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

# Day 10 Problem Set

Andrew Loeppky

ATSC 405

+++

### Question 1 

Derive equation 7 in my Scale height notes (feel free to use the linked indefinite integral).

+++

### Question 2

Show starting from the first law that the entropy of liquid water is

$$
\phi_l = c_llog(T)
$$

stating all your assumptions.

---

+++

*Start by writing out the first law and the definition of entropy:*

$$
du = dq - pdv\tag{First Law}
$$

$$
d\phi = \frac{dq}{T}\tag{AT's Entropy Def}
$$

*Liquid water is incompressible and the thermal expansion coefficient is small as far as atmospheric scientists are concerned, so $dv=0$, leaving:*

$$
du = dq
$$

*Assume all the internal energy goes to raising the temperature (not changing phase or something else):*

$$
du = dq = c_ldT
$$

*Divide through by $T$ to get the LHS equal to the entropy:*

$$
d\phi = \frac{dq}{T} = c_l\frac{dT}{T}
$$

*Integrate from some reference state $(T_0,\phi_0)$ to $(T,\phi)$:*

$$
\int_{\phi_0}^{\phi_l} d\phi = c_l\int_{T_0}^T \frac{dT'}{T'}
$$

$$
\phi_l = c_llog\left(\frac{T}{T_0}\right)\square
$$
