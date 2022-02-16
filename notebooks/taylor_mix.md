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

# Air Parcel Mixing Approximation

+++

On p. 57, Thompkins shows that the moist static energy mixies linearly. Show that that is
also roughly true for $\theta_e$ by doing the following:

1. Recall that, if we neglect the increase in entropy due to irreversible mixing, that the entropies
of two components will sum to give the entropy of the mixture:

$$
s_{mix} = fs_{env} + (1 - f)s_{cld}\tag{1}
$$

where $f$ and $(1 - f)$ are the mass fractions of the environment and cloudy air and $s_{env}$ and
$s_{cld}$ are the specifc entropies of each component in J/kg/K. We also know that, if we neglect the effect of liquid and vapor on the heat capacities:

$$
s = c_p log(\theta_e)\tag{2}
$$

so we can rewrite (1) as

$$
log \theta_{mix} \approx f log(\theta_{e,env}) + (1 - f) log(\theta_{e,cld})\tag{3}
$$

or equivalently, taking exp of each side:

$$
\theta_{e,mix} = \theta_{e,env}^f \theta_{e,cld}^{(1-f)}\tag{4}
$$


Use a Taylor series expansion to show that if $(\theta_{e,cld} - \theta_{e,env})/\theta_{e,env}= \delta \ll 1$ then (4) is
approximately equivalent to:

$$
\theta_{e,mix} = f\theta_{e,env} + (1 - f)\theta_{e,cld}\tag{5}
$$


*Hint: see if ou can rewrite $\theta_{e,cld}/\theta_{e,env}$ as $(\theta_{e,env} +(\theta_{e,cld} - \theta_{e,env}))/\theta_{e,env}$ and show by expanding in a Taylor series that:*

$$
(1 + \delta)^f \approx 1 + f\delta\tag{6}
$$

*and*

$$
f + (1 - f) = 1 \tag{7}
$$
