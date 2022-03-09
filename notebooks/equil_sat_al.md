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

# Equalibrium Saturation Derivation - Andrew's Solution

Derivation of (6) from Phil's Equilibrium Saturation notes, starting from:

+++

$$
r_v = \underbrace{(1 + SS)}_{A}\underbrace{\frac{\epsilon e_s}{p}}_{B}\tag{PA 5}
$$

Derive equation (6). Differentiating (PA 5) using the chain rule yields:

$$
\frac{dr_v}{dt} = \underbrace{(1 + SS)\left(\underbrace{\frac{d}{dt}\left[\frac{1}{p}\right]}_{I}\epsilon e_s + \frac{\epsilon}{p}\underbrace{\frac{de_s}{dt}}_{II} \right)}_{A\frac{dB}{dt}} + \underbrace{\frac{\epsilon e_s}{p}\underbrace{\frac{d}{dt}(1 + SS)}_{III}}_{B\frac{dA}{dt}}
$$

$(I)$:

$$
\frac{d}{dt}\left[\frac{1}{p}\right] = \frac{d(1/p)}{dp}\frac{dp}{dt} = \frac{-1}{p^2}\frac{dp}{dt} 
$$

From the notes, use equation 7b:

$$
\frac{dp}{dt} = -\rho g \frac{dz}{dt} = \frac{-gpV}{R_dT}\tag{PA 7b}
$$

$$
(I) = \frac{-\epsilon e_s}{p^2}\left(\frac{-gpV}{R_dT}\right)
$$

$(II)$:

$$
\frac{de_s}{dt} = \frac{de_s}{dT}\frac{dT}{dt}
$$

using equation 7a:

$$
\frac{de_s}{dT} = \frac{\epsilon Le_s}{R_dT^2}\tag{PA 7a}
$$

$$
(II) = \frac{\epsilon Le_s}{R_dT^2}\frac{dT}{dt}
$$

$(III)$:

$$
\frac{d}{dt}(1 + SS) = \frac{dSS}{dt}
$$

Combine $I, II, III$ to get (PA 6):

$$
\frac{dr_v}{dt} = (1 + SS)\left[\frac{-\epsilon e_s}{p^2}\left(\frac{-gpV}{R_dT}\right) + \frac{\epsilon}{p}\left(\frac{\epsilon Le_s}{R_dT^2}\right)\frac{dT}{dt}\right] + \frac{\epsilon e_s}{p}\frac{dSS}{dt}\tag{PA 6}\square
$$
