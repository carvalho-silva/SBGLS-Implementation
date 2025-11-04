import LS_like_periodograms as lslp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('data-sets/DS1_timeSeries.csv')

col_dict = {
    'Standard File Name' : 'file',
    'Time [eMJD]' : 'time',
    'RV [m/s]' : 'rv',
    'RV Err. [m/s]' : 'e_rv',
    'Exp. Time [s]' : 'exptime',
    'Airmass' : 'airmass',
    'BERV [km/s]' : 'berv',
    'Instrument' : 'inst',
    'CCF FWHM [km/s]' : 'fwhm',
    'CCF FWHM Err. [km/s]' : 'e_fwhm',
    'CCF Contrast' : 'contrast',
    'CCF Contrast Err.' : 'e_contrast',
    'BIS [m/s]' : 'bis',
    'H-alpha Emission' : 'ha',
    'CaII Emission' : 'caii'
}

renamed_df = df.rename(columns=col_dict) # rename the columns
renamed_df.head()

def sine_func(params, t, P):
	a, b, c = params
	omega = 2 * np.pi / P
	return a * np.sin(omega * t) + b * np.cos(omega * t) + c
	
def func_fit(t, rv, P):
	omega = 2 * np.pi / P
	t_phase = omega * t

	X = np.column_stack([
		np.sin(t_phase),
		np.cos(t_phase),
		np.ones_like(t_phase)
	])

	coef, _, _, _ = np.linalg.lstsq(X, rv, rcond=None)
	a, b, c = coef
	print("a, b, c =", a, b, c)

	
	rv_fit = sine_func(coef, t, P)
	
	residuals = rv - rv_fit
	
	return coef, rv_fit, residuals, omega

   
time = renamed_df['time'] 
rv = renamed_df['rv']

P = 2.92

coef, rv_fit, residuals, omega = func_fit(time, rv, P)

t_ = np.linspace(time.min(),time.max(),500)

phase = (time % P) / P


idx = np.argsort(phase)
phase_sorted = phase[idx]
rv_fit_sorted = rv_fit[idx]

plt.scatter(phase, rv) 
plt.plot(phase_sorted, rv_fit_sorted)
plt.show()

renamed_df['rv_fit'] = rv_fit 
renamed_df['residuals'] = residuals
renamed_df.to_csv('DS1_timeSeries_b.csv')
print(renamed_df)
