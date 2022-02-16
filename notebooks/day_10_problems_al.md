---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Day 10 Problem Set - Andrew's Solution

Andrew Loeppky

ATSC 405

+++

### Question 1 

Derive equation 7 in my Scale height notes (feel free to use the linked indefinite integral).

---

+++

*Expected value of a continuous random variable $x$ on $0 < x < \infty$:*

$$
\bar{x} = \int_0^\infty xp(x)dx
$$

*Exponential distribution:*

$$
p(x) = \lambda exp(-\lambda x)
$$

*Cheat by looking up the integral [on Wikipedia](https://en.wikipedia.org/wiki/List_of_integrals_of_exponential_functions)*

$$
\lambda\int xe^{-\lambda x}dx = \lambda e^{-\lambda x}\left(\frac{\lambda x - 1}{\lambda^2}\right)
$$

$$
= e^{-\lambda x} \left(\frac{\lambda x - 1}{\lambda}\right)
$$

$$
= e^{-\lambda x} \left(1 - \frac{1}{\lambda}\right)
$$

*Evaluate at the limits $x = 0$ and $x=\infty$:*

$$
\bar{x} = e^{-\infty}(\infty - 1/\lambda) - e^0(0 - 1/\lambda)
$$

$$
\bar{x} = \frac{1}{\lambda}\square
$$

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
