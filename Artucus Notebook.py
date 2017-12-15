
# coding: utf-8

# # Artucus Spectra
# 
# Infrared Arcturus Atlas (Hinkle+ 1995)
# These are currently not telluric corrected but you can find some that are
# Resolving power of 100000,

# In[1]:



import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import string
import numpy as np
from astropy.io import fits

from astro_scripts.plot_fits import get_wavelength


# In[2]:


files_1 = glob.glob("/home/jneal/Phd/data/artucus/1*.fits")
files_2 = glob.glob("/home/jneal/Phd/data/artucus/2*.fits")
for f in files_1:
    print(f)


# In[3]:


plt.figure(figsize=(18, 8))
for f in files_1:
    data, hdr = fits.getdata(f, header=True)
    wave = get_wavelength(hdr, convert=False) / 10
    plt.xlabel("Wavelenght (nm)")
    
    plt.plot(wave, data)
plt.title("Artucus at 1000nm")
plt.show()


# In[4]:


plt.figure(figsize=(18, 8))
for f in files_2:
    data, hdr = fits.getdata(f, header=True)
    wave = get_wavelength(hdr, convert=False) / 10
    
    plt.plot(wave, data)
plt.title("Artucus at 2000nm")
plt.xlabel("Wavelenght (nm)")
plt.show()


# In[5]:


print(fits.getheader(f))

