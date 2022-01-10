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
<div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Getting-started-with-vapor" data-toc-modified-id="Getting-started-with-vapor-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Getting started with vapor</a></span></li></ul></div>

+++

# Getting started with vapor

+++

1.  Fetch the netcdf file tropical.nc from my webserver (50 Mbytes)

```{code-cell} ipython3
!mv  tropical_clouds.zip test.zip
```

```{code-cell} ipython3
from a405.utils.data_read import download
the_root = "http://clouds.eos.ubc.ca/~phil/courses/atsc500/docs"
the_file =  "tropical_clouds.zip"
out=download(the_file,root=the_root)
```

2.  unzip tropical_clouds.zip in your notebook folder

+++

3. Install vapor per [these instructions](https://www.vapor.ucar.edu/docs/usage/getting-started-vapor)

+++

Start vapor, go to "Data -> Load dataset into current session" and point it at tropical_clouds/outname.vdf  See if you can get a 3-d visualization of 3-d visualization of TABS

```{code-cell} ipython3

```
