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
<div class="toc"><ul class="toc-item"><li><span><a href="#Demonstrate-the-sounding-retrieval-code" data-toc-modified-id="Demonstrate-the-sounding-retrieval-code-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Demonstrate the sounding retrieval code</a></span></li><li><span><a href="#Read-the-soundings-back-into-python" data-toc-modified-id="Read-the-soundings-back-into-python-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Read the soundings back into python</a></span></li><li><span><a href="#Examine-the-nested-dictionaries-inside-soundings" data-toc-modified-id="Examine-the-nested-dictionaries-inside-soundings-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Examine the nested dictionaries inside soundings</a></span></li><li><span><a href="#Get-the-first-sounding" data-toc-modified-id="Get-the-first-sounding-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Get the first sounding</a></span></li><li><span><a href="#Plot-it" data-toc-modified-id="Plot-it-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Plot it</a></span></li></ul></div>

+++

# Demonstrate the sounding retrieval code

```{code-cell} ipython3
from a405.soundings.wyominglib import write_soundings, read_soundings
from matplotlib import pyplot as plt
```

Ask for north american soundings between July 1, 2017 00Z and July 18, 2017 00Z for
Dodge City, Kansas  from http://weather.uwyo.edu/upperair/sounding.html

```{code-cell} ipython3
values=dict(region='naconf',year='2017',month='5',start='0100',stop='1800',station='72261')
```

Write the soundings into a folder called soundingdir

```{code-cell} ipython3
write_soundings(values, 'delrio')
```

# Read the soundings back into python

```{code-cell} ipython3
soundings= read_soundings('delrio')
```

# Examine the nested dictionaries inside soundings

```{code-cell} ipython3
print((f'soundings keys: {list(soundings.keys())}\n'),
      (f'soundings attributes: {list(soundings["attributes"])}\n'),
      (f'sounding_dict keys: {list(soundings["sounding_dict"].keys())}'))
```

# Get the first sounding

```{code-cell} ipython3
target_date=list(soundings['sounding_dict'].keys())[0]
the_sounding = soundings['sounding_dict'][target_date]
print(the_sounding.columns)
```

# Plot it

```{code-cell} ipython3
%matplotlib inline
m2km=1.e-3  #convert meters to km
fig,ax=plt.subplots(1,1,figsize=(8,10))
ax.plot(the_sounding['temp'],the_sounding['hght']*m2km,label='temp')
ax.plot(the_sounding['dwpt'],the_sounding['hght']*m2km,label='dewpoint')
ax.legend()
out=ax.set(xlabel="temperature (K)",ylabel="height (km)",
      title ="First sounding")
ax.grid(True)
```

```{code-cell} ipython3

```
