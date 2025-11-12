import remove_planet as rp
import pandas as pd



DSno_list = [str(i) for i in range(1, 10)]

for DSno in DSno_list:
	DSno_str = str(DSno)
	bestP_results = pd.read_csv(f'bestP_results_DS{DSno_str}_b.csv')
	bestP_results['DSno'] = bestP_results['DSno'].astype(str) 
	
	mask = (bestP_results['mode']=='random without replacement')&(bestP_results['ptype']=='BGLS')
	best_P = bestP_results[mask].loc[bestP_results.DSno == DSno,'best_P_SNR'].iloc[0]
	print('DS: ',DSno_str)
	print('best P = ',best_P)
	
	df = pd.read_csv(f'DS{DSno}_timeSeries_b.csv')

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
		'CaII Emission' : 'caii',
		'residuals' : 'rv_b'
	}

	renamed_df = df.rename(columns=col_dict)
	renamed_df.head()
	time = renamed_df.time

	RV = renamed_df.rv_b
	RV_error = renamed_df.e_rv


	rp.remove_planet(renamed_df, best_P, DS_no=DSno_str, planet='c',save_plots=True)



