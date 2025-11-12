import remove_planet as rp
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

df = pd.read_csv('data-sets/DS1_timeSeries.csv')

col_dict = {
    'Standard File Name' : 'file',
    'Time [eMJD]'        : 'time',
    'RV [m/s]'           : 'rv',
    'RV Err. [m/s]'      : 'e_rv',
    'Exp. Time [s]'      : 'exptime',
    'Airmass'            : 'airmass',
    'BERV [km/s]'        : 'berv',
    'Instrument'         : 'inst',
    'CCF FWHM [km/s]'    : 'fwhm',
    'CCF FWHM Err. [km/s]' : 'e_fwhm',
    'CCF Contrast'       : 'contrast',
    'CCF Contrast Err.'  : 'e_contrast',
    'BIS [m/s]'          : 'bis',
    'H-alpha Emission'   : 'ha',
    'CaII Emission'      : 'caii'
}

df = df.rename(columns=col_dict)

t = df.time
rv_obs = df.rv

rp.remove_planet(t, rv_obs, 2.91, planet='b', save_plots='teste')
