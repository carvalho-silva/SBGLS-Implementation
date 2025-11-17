import pandas as pd
import LS_like_periodograms as aux



type_list = ['BGLS']

mode_list = ['chronological','random','random without replacement']

DSno_list = [str(i) for i in range(1, 10)]
#DSno_list = DSno_list = [str(2)]

for DSno in DSno_list:
	results = []
	DSno_str = str(DSno)
	df = pd.read_csv(f'data-sets/DS{DSno}_timeSeries.csv')

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
	time = renamed_df.time

	RV = renamed_df.rv
	RV_error = renamed_df.e_rv
	for ptype in type_list:
		periodogram_type = ptype
		for mode in mode_list:
			mode=mode


			data = aux.stacked_periodogram(time, RV, RV_error, N_min=55,
						                   periodogram_type=periodogram_type, p_min=1.,
						                   p_max=60., num_periods = 5000,
						                   mode=mode, exclude_periods=(0.98,1.02))

			aux.plot_stacked_periodogram_heatmap(data, norm='log', highlight_strong_signal = True,
			plot_SNR=True, save_plots=f'output/{periodogram_type}_DS{DSno}_{mode}')
			
			best_P = data['SNR']['best_P']
						                   
			results.append({
                'DSno'        : DSno_str,
                'ptype'       : periodogram_type,
                'mode'        : mode,
                'best_P_SNR'  : best_P
            })
	results_df = pd.DataFrame(results)
	results_df.to_csv(f'output/bestP_results_DS{DSno}.csv', index=False, sep=',')
