
# coding: utf-8

# # Comparing Phoenix Spectra

# Plot Phoenix aces and Phoenix BT-Settl for same parameters
# 
# Also AMES Dusty Cond.

# #### FORMAT OF THE SPECTRA OUTPUT FILES
# 
# You can find the  pre-computed grids, also accessible via links on
# the bottom part of the simulator presentation page, or using this link:
# http://phoenix.ens-lyon.fr/Grids/
# 
# The file names contain the main parameters of the models:
# lte{Teff/10}-{Logg}{[M/H]}a[alpha/H].GRIDNAME.7.spec.gz/bz2/xz
# is the synthetic spectrum for the requested effective temperature
# (Teff),surface gravity (Logg), metallicity by log10 number density with
# respect to solar values ([M/H]), and alpha element enhencement relative     
# to solar values [alpha/H]. The model grid is also mentioned in the name.
# 
# Spectra are provided in an ascii format (\*.7.gz):
# 
# column1: wavelength in Angstroem
# column2: 10\*\*(F_lam + DF) to convert to Ergs/sec/cm\*\*2/A
# column3: 10\*\*(B_lam + DF) i.e. the blackbody fluxes of same Teff in same units.
# 
# Additional columns, obtained systematically when computing spectra using the
# Phoenix simulator, give the information to identify atomic and molecular
# lines. This information is used by the idl scripts lineid.pro and plotid.pro 
# which are provided in the user result package.  
#    
# With the stacked ascii format (\*.spec.gz files ) we have rather:
# 
# line1: Teff logg [M/H] of the model
# line2: number of wavelengths
# line3: F_lam(n) X 10\*\*DF , n=1,number of wavelengths
# lineX: B_lam(n) X 10\*\*DF , n=1,number of wavelengths
# 

# In[1]:


import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import string
import numpy as np
from astropy.io import fits
from spectrum_overload import Spectrum
from loading_phoenix import load_phoenix_aces, load_Allard_Phoenix, align2model

get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


os.chdir("/home/jneal/Phd/2017/Compare_PHOENIX/")
os.getcwd()
files = glob.glob("data/*")
# print(files)


# In[ ]:



def plot_allard_phoenix(fname, band=None):
    wav, flux, bb_flux = load_Allard_Phoenix(fname)
    plt.figure(figsize=(10, 8))
    plt.plot(wav, flux, label="flux")
    plt.plot(wav, bb_flux, label="bb_flux")
    plt.legend()
    plt.title(fname)
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")


# In[ ]:


from PyAstronomy.pyasl.phoenixUtils.read import readUnit7, readDTable, decomposeFilename
#readUnit7 is in units of per cm not per A

decomposeFilename("lte4300-2.50-0.0a+0.0.BT-dusty-giant-2013.cf128.sc.spid.fits")




read7 = readUnit7("data/lte043.0-2.5-0.0a+0.0.BT-Settl.spec.7")
print(read7.shape)
print(read7)
plt.figure(figsize=(10, 8))
plt.plot(read7[:, 0]/10, 10**read7[:, 1], "--", label="read7")
plt.plot(read7[:, 0]/10, 10**read7[:, 2], "--", label="read7 bb ")
plt.legend()
plt.xlim([700, 5000])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Ergs/sec/cm**2/A")
plt.title("BT-Settle 4300-2.5-0.0")
plt.show()


# In[ ]:


wav, flux, bb_flux = load_Allard_Phoenix("data/lte043.0-2.5-0.0a+0.0.BT-Settl.spec.7")
read7 = readUnit7("data/lte043.0-2.5-0.0a+0.0.BT-Settl.spec.7")
plt.figure(figsize=(10, 8))
plt.plot(wav, flux, label="flux")
plt.plot(read7[0], read7[1], "--", label="read7")
plt.plot(wav, bb_flux, label="bb_flux")
plt.plot(read7[0], read7[1], "--", label="read7 bb ")
plt.legend()
plt.xlim([700, 5000])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Ergs/sec/cm**2/A")
plt.title("BT-Settle 4300-2.5-0.0")
plt.show()


# In[ ]:


plot_allard_phoenix("data/lte043.0-2.5-0.0a+0.0.BT-Settl.spec.7")
plt.show()


# In[ ]:


wav, flux, bb_flux = load_Allard_Phoenix("data/lte430-3.5-0.5a+0.2.BT-Cond.7")

plt.figure(figsize=(10, 8))
plt.plot(wav, flux, label="flux")
plt.plot(wav, bb_flux, label="bb_flux")
plt.legend()
plt.xlim([700, 5000])
plt.ylim([1e1, 1e10])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Ergs/sec/cm**2/A")
plt.title("BT-Cond 4300-3.5-0.5")
#plt.ylim([10, 18])
plt.show()


# In[ ]:



wav, flux, bb_flux = load_Allard_Phoenix("data/lte043-2.5-0.0.BT-Dusty.spec.7")

plt.figure(figsize=(10, 8))
plt.plot(wav, flux, label="flux")
plt.plot(wav, bb_flux, ".", label="bb_flux")
plt.legend()
plt.xlim([700, 5000])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Ergs/sec/cm**2/A")
plt.title("BT-Dusty 4300-2.5-0.0")
plt.legend()
plt.show()


# In[ ]:



wav, flux, bb_flux = load_Allard_Phoenix("data/lte043.0-2.5-0.0a+0.0.BT-Settl.spec.7")

plt.figure(figsize=(10, 8))
plt.plot(wav, flux, label="flux")
plt.plot(wav, bb_flux, label="bb_flux")
plt.legend()
plt.xlim([700, 5000])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Ergs/sec/cm**2/A")
plt.title("BT-Dusty 4300-2.5-0.0")
plt.legend()
plt.show()


# ### COMPARING MODELS

# In[ ]:


Artucus = [4300, 1.50, -0.5]
HD30501 = [5200, 4.5, 0.0]
ACES_bottom = [2300, 4.5, 0.0]
Sun = [5800, 4.5, 0.0]


# In[ ]:


# Teff 5800, 4.5, 0.0
# comparison_plot("Sun", *Sun)
from spectrum_overload import Spectrum


w_next, f_next, bb_next = load_Allard_Phoenix("data/lte580-4.5-0.0a+0.0.BT-NextGen.7")
next_spec = Spectrum(xaxis=w_next, flux=f_next) 
w_dusty_spec, f_dusty_spec, bb_dusty_spec = load_Allard_Phoenix("data/lte058-4.5-0.0.BT-Dusty.spec.7")
dusty_spec = Spectrum(xaxis=w_dusty_spec, flux=f_dusty_spec) 
w_settl, f_settl, bb_settl = load_Allard_Phoenix("data/lte058.0-4.5-0.0a+0.0.BT-Settl.spec.7")
settl_spec = Spectrum(xaxis=w_settl, flux=f_settl) 
w_cond, f_cond, bb_cond = load_Allard_Phoenix("data/lte580-4.5-0.0a+0.0.BT-Cond.7")
cond_spec = Spectrum(xaxis=w_cond, flux=f_cond) 
w_aces, f_aces = load_phoenix_aces("data/lte05800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")
aces_spec = Spectrum(xaxis=w_aces, flux=f_aces) 
w_dusty_fits, f_dusty_fits, bb_dusty_fits = load_Allard_Phoenix("data/lte5800-4.50-0.0a+0.0.BT-dusty-giant-2013.cf128.sc.spid.fits")
w_dusty_fits = w_dusty_fits*1000
dusty_fits_spec = Spectrum(xaxis=w_dusty_fits, flux=f_dusty_fits) 

plt.figure(figsize=(15, 10))
plt.plot(w_dusty_spec, f_dusty_spec/max(f_dusty_spec), label="Dusty spec")
plt.plot(w_settl, f_settl/max(f_settl), label="BT-SETTL")
plt.plot(w_cond, f_cond/max(f_cond), label="BT-COND")
plt.plot(w_aces, f_aces/max(f_aces), label="PHOENIX ACES")
plt.plot(w_dusty_fits, f_dusty_fits/max(f_dusty_fits), label="Dusty fits")
plt.plot(w_next, f_next/max(f_next), label="NEXTGEN")
plt.title("Sun - 5800K")
plt.xlabel("Wavelength(nm)")
plt.ylabel("Flux")
plt.legend()
plt.show()


# In[ ]:


limits = [2116, 2120]
next_spec.wav_select(*limits)
#print(next_spec.xaxis, next_spec.flux)
next_spec = next_spec.normalize("exponential")
dusty_spec.wav_select(*limits)
dusty_spec = dusty_spec.normalize("exponential")
settl_spec.wav_select(*limits)
settl_spec = settl_spec.normalize("exponential")
cond_spec.wav_select(*limits)
cond_spec = cond_spec.normalize("exponential")
aces_spec.wav_select(*limits)
aces_spec = aces_spec.normalize("exponential")
dusty_fits_spec.wav_select(*limits)
dusty_fits_spec = dusty_fits_spec.normalize("exponential")


# In[ ]:



next_spec.plot(label="nextgen")
dusty_spec.plot(label="dusty_spec")
cond_spec.plot(label="cond")
aces_spec.plot(label="aces")
dusty_fits_spec.plot(label="dusty_fits")
settl_spec.plot(label="settl")
plt.legend()
plt.show()


# In[ ]:



next_spec.plot(label="nextgen")
plt.legend()
plt.show()
dusty_spec.plot(label="dusty_spec")
plt.legend()
plt.show()
settl_spec.plot(label="settl")
plt.legend()
plt.show()
cond_spec.plot(label="cond")
plt.legend()
plt.show()
aces_spec.plot(label="aces")
plt.legend()
plt.show()
dusty_fits_spec.plot(label="dusty_fits")

plt.legend()
plt.show()

dusty_spec.plot(label="dusty_spec")
dusty_fits_spec.plot(linestyle="--",label="dusty_fits")
plt.legend()
plt.show()


# In[ ]:


aces_spec.plot(label="aces")
dusty_spec.plot(linestyle="-.", label="dusty_spec")
settl_spec.plot(linestyle=":", color="k", label="settl")
dusty_fits_spec.plot(linestyle="--", label="dusty_fits")

plt.legend()
plt.show()


# ## PHONEIX ACES Limit - Teff 2300, 4.5, 0.0

# In[ ]:



# comparison_plot("ACES_bottom", *ACES_bottom)

w_next, f_next, bb_next = load_Allard_Phoenix("data/lte023-5.0-0.0.BT-NextGen.7")
next_spec = Spectrum(xaxis=w_next, flux=f_next) 

w_dusty_spec, f_dusty_spec, bb_dusty_spec = load_Allard_Phoenix("data/lte023-4.5-0.0.BT-Dusty.spec.7")
dusty_spec = Spectrum(xaxis=w_dusty_spec, flux=f_dusty_spec) 

w_settl, f_settl, bb_settl = load_Allard_Phoenix("data/lte023.0-4.5-0.0a+0.0.BT-Settl.spec.7")
settl_spec = Spectrum(xaxis=w_settl, flux=f_settl) 

w_cond, f_cond, bb_cond = load_Allard_Phoenix("data/lte230-4.5-0.0a+0.0.BT-Cond.7")
cond_spec = Spectrum(xaxis=w_cond, flux=f_cond) 

w_aces, f_aces = load_phoenix_aces("data/lte02300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")
aces_spec = Spectrum(xaxis=w_aces, flux=f_aces) 
w_dusty_fits, f_dusty_fits, bb_dusty_fits = load_Allard_Phoenix("data/lte5800-4.50-0.0a+0.0.BT-dusty-giant-2013.cf128.sc.spid.fits")
w_dusty_fits = w_dusty_fits*1000
dusty_fits_spec = Spectrum(xaxis=w_dusty_fits, flux=f_dusty_fits) 


limits = [2100, 2200]
next_spec.wav_select(*limits)
# print(next_spec.xaxis, next_spec.flux)
next_spec = next_spec.normalize("exponential")
dusty_spec.wav_select(*limits)
dusty_spec = dusty_spec.normalize("exponential")
settl_spec.wav_select(*limits)
settl_spec = settl_spec.normalize("exponential")
cond_spec.wav_select(*limits)
cond_spec = cond_spec.normalize("exponential")
aces_spec.wav_select(*limits)
aces_spec = aces_spec.normalize("exponential")
dusty_fits_spec.wav_select(*limits)
dusty_fits_spec = dusty_fits_spec.normalize("exponential")



# In[ ]:



aces_spec.plot(label="aces")
dusty_spec.plot(linestyle="-.", label="dusty_spec")
settl_spec.plot(linestyle=":", color="k", label="settl")

plt.title("Spectra at 2300 K")
plt.legend()
plt.show()


# In[ ]:



aces_spec.plot(label="aces")
dusty_spec.plot(linestyle="-.", label="dusty_spec")
settl_spec.plot(linestyle=":", color="k", label="settl")

plt.title("Spectra at 2300 K")
plt.xlim([2120, 2120.5])
plt.legend()
plt.show()


# In[ ]:





# In[ ]:


aces_spec.plot(label="aces")
dusty_spec.plot(linestyle="-.", label="dusty_spec")
settl_spec.plot(linestyle=":", color="k", label="settl")
dusty_fits_spec.plot(linestyle="--", label="dusty_fits")

plt.legend()
plt.show()


# In[ ]:


# Teff 5200, 4.5, 0.0
# HD30501 Host

w_next, f_next, bb_next = load_Allard_Phoenix("data/lte520-4.5-0.0a+0.0.BT-NextGen.7")
next_spec = Spectrum(xaxis=w_next, flux=f_next) 

w_dusty, f_dusty, bb_dusty = load_Allard_Phoenix("data/lte052-4.5-0.0.BT-Dusty.spec.7")
dusty_spec = Spectrum(xaxis=w_dusty, flux=f_dusty) 

w_settl, f_settl, bb_settl = load_Allard_Phoenix("data/lte052.0-4.5-0.0a+0.0.BT-Settl.spec.7")
settl_spec = Spectrum(xaxis=w_settl, flux=f_settl) 

w_cond, f_cond, bb_cond = load_Allard_Phoenix("data/lte520-4.5-0.5a+0.0.BT-Cond.7")
cond_spec = Spectrum(xaxis=w_cond, flux=f_cond) 

w_aces, f_aces = load_phoenix_aces("data/lte05200-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")
aces_spec = Spectrum(xaxis=w_aces, flux=f_aces) 

w_dusty_fits, f_dusty_fits, bb_dusty_fits = load_Allard_Phoenix("data/lte5200-4.50-0.0a+0.0.BT-dusty-giant-2013.cf128.sc.spid.fits")
w_dusty_fits = w_dusty_fits*1000
dusty_fits_spec = Spectrum(xaxis=w_dusty_fits, flux=f_dusty_fits) 



# In[ ]:



limits = [2100, 2200]
next_spec.wav_select(*limits)
print(next_spec.xaxis, next_spec.flux)
next_spec = next_spec.normalize("exponential")
dusty_spec.wav_select(*limits)
dusty_spec = dusty_spec.normalize("exponential")
settl_spec.wav_select(*limits)
settl_spec = settl_spec.normalize("exponential")
cond_spec.wav_select(*limits)
cond_spec = cond_spec.normalize("exponential")
aces_spec.wav_select(*limits)
aces_spec = aces_spec.normalize("exponential")
dusty_fits_spec.wav_select(*limits)
dusty_fits_spec = dusty_fits_spec.normalize("exponential")


# In[ ]:


aces_spec.plot(label="aces")
dusty_spec.plot(linestyle="-.", label="dusty_spec")
settl_spec.plot(linestyle=":", color="k", label="settl")
dusty_fits_spec.plot(linestyle="--", label="dusty_fits")

plt.legend()
plt.show()


# # Comparing to ARTUCUS
# 
# it is 100000 Resolving power

# In[ ]:


files_1 = glob.glob("/home/jneal/Phd/data/artucus/1*.fits")
files_2 = glob.glob("/home/jneal/Phd/data/artucus/2*.fits")
for f in files_1:
    print(f)
for f in files_2:
    print(f)


# In[ ]:


from astro_scripts.plot_fits import get_wavelength


for f in files_1:
    data, hdr = fits.getdata(f, header=True)
    wave = get_wavelength(hdr, convert=False)
    
    plt.plot(wave, data)
plt.title("Artucus at 1000nm")
plt.show()

for f in files_2:
    data, hdr = fits.getdata(f, header=True)
    wave = get_wavelength(hdr, convert=False)
    
    plt.plot(wave, data)
plt.title("Artucus at 2000nm")
plt.show()


# In[ ]:



artucus_1 = "/home/jneal/Phd/data/artucus/10097-10155_s-obs.fits"
data, hdr = fits.getdata(artucus_1, header=True)

artucus_1 = Spectrum(xaxis=get_wavelength(hdr)/10, flux=data, header=hdr )

artucus_2 = "/home/jneal/Phd/data/artucus/21380-21518_s-obs.fits"
data, hdr = fits.getdata(artucus_2, header=True)

artucus_2 = Spectrum(xaxis=get_wavelength(hdr)/10, flux=data, header=hdr)


# # Comparing Artucus at 1um

# In[ ]:


artucus_1.plot()
plt.show()


# In[ ]:


artucus_2.plot()
plt.show()


# In[ ]:


# Teff 4300, 1.5 logg, -0.5 [Fe/H]
# comparison_plot("Artucus", *Artucus)

w_next, f_next, bb_next = load_Allard_Phoenix("data/lte043-2.5-0.0a+0.0.BT-NextGen.7")
next_spec = Spectrum(xaxis=w_next, flux=f_next) 

w_dusty, f_dusty, bb_dusty = load_Allard_Phoenix("data/lte043-2.5-0.0.BT-Dusty.spec.7")
dusty_spec = Spectrum(xaxis=w_dusty, flux=f_dusty) 

w_settl, f_settl, bb_settl = load_Allard_Phoenix("data/lte043.0-2.5-0.0a+0.0.BT-Settl.spec.7")
settl_spec = Spectrum(xaxis=w_settl, flux=f_settl) 

w_cond, f_cond, bb_cond = load_Allard_Phoenix("data/lte043-2.5-0.0a+0.0.BT-Cond.7")
cond_spec = Spectrum(xaxis=w_cond, flux=f_cond) 

w_aces, f_aces = load_phoenix_aces("data/lte04300-1.50-0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")
aces_spec = Spectrum(xaxis=w_aces, flux=f_aces) 

#w_dusty_fits, f_dusty_fits, bb_dusty_fits = load_Allard_Phoenix("data/lte4300-2.50-0.0a+0.0.BT-dusty-giant-2013.cf128.vo0.spid.fits")
#w_dusty_fits = w_dusty_fits*1000
#dusty_fits_spec = Spectrum(xaxis=w_dusty_fits, flux=f_dusty_fits) 


# In[ ]:


limits = [2100, 2200]
limits = [artucus_1.xaxis[0], artucus_1.xaxis[-1]]
next_spec.wav_select(*limits)
print(next_spec.xaxis, next_spec.flux)
next_spec = next_spec.normalize("exponential")
dusty_spec.wav_select(*limits)
dusty_spec = dusty_spec.normalize("exponential")
settl_spec.wav_select(*limits)
settl_spec = settl_spec.normalize("exponential")
cond_spec.wav_select(*limits)
cond_spec = cond_spec.normalize("exponential")
aces_spec.wav_select(*limits)
aces_spec = aces_spec.normalize("exponential")
dusty_fits_spec.wav_select(*limits)
dusty_fits_spec = dusty_fits_spec.normalize("exponential")


# In[ ]:


artucus_1.plot("Artucus")
aces_spec.plot(label="aces")
dusty_spec.plot(linestyle="-.", label="dusty_spec")
settl_spec.plot(linestyle=":", color="k", label="settl")
dusty_fits_spec.plot(linestyle="--", label="dusty_fits")

plt.legend()
plt.show()


# In[ ]:


# BTSETTL RESOLUTION
# The spectral resolution used to compute the *.spec.7.bz2 spectra is following or better:
# 0.02A in the optical, and 0.05A in the infrared for the most recent files.

# At 2000 nm this amounts to a Resolution of  R = lambda/delta lambda
wav = 20000  # Angstrom
dwav = 0.05
print("Resolution of {}A at {}A gives, resolving power R={}".format(dwav, wav, wav/dwav))

# Need to convovlve to a give resolution.


# 
#  BT-Cond Teff =2600K to 70000K
# 
# 
# from FORMAT
# 
# In the case of the most recent models, the Barber & Tennison (UCL) so-called BT2 water vapor line list
# has been used explaining why all those models bear names starting with 'BT-'. 
# BT-Dusty refers to dust in equilibrium with the gas phase (sedimentation is
# neglected), while BT-Cond includes dust condensation in equilibrium with the
# gas phase while neglecting their opacities in the radiative transfer. BT-Settl
# means that gravitational settling of sedimentation is accounted for in the
# frame of a detailed cloud model 
