import pandas as pd
import LS_like_periodograms as aux
import os
import numpy as np


TYPE_LIST = ['GLS','BGLS']

MODE_LIST = ['chronological','random','random without replacement']

DSNO_LIST = [str(i) for i in range(1, 10)]

COL_DICT = {
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

PROXIES_LIST = ['rv','bis', 'fwhm', 'contrast', 'ha', 'caii']


for DSno in DSNO_LIST:

    path = f'data-sets/DS{DSno}_timeSeries.csv'
    if not os.path.exists(path):
        print(f'File {path} not found, skipping...')
        continue
    df = pd.read_csv(path)
    all_results = []
    print(f'Processing DS{DSno}...')

    renamed_df = df.rename(columns=COL_DICT)

    time = renamed_df['time'].values
    for proxy in PROXIES_LIST:
        if proxy not in renamed_df.columns:
            print(f"  '{proxy}' not found in DS{DSno}, skipping...")
            continue
        y = renamed_df[proxy].values
        if f'e_{proxy}' in renamed_df.columns:
            yerr = renamed_df[f'e_{proxy}'].values
        else:
            yerr = np.ones_like(y)
        for ptype in TYPE_LIST:
                print(f'    Processing proxy {proxy} in {ptype} periodogram...')
                for mode in MODE_LIST:
                    print(f'    Mode: {mode}')
                    data = aux.stacked_periodogram(
                                                    time,
                                                    y,
                                                    yerr,
                                                    N_min=55,
                                                    periodogram_type=ptype,
                                                    p_min=1.5,
                                                    p_max=60., 
                                                    num_periods = 5000,
                                                    mode=mode
                                                    )
                    
                    save_name = f'{proxy}_{ptype}_DS{DSno}_{mode}'

                    aux.plot_stacked_periodogram_heatmap(
                                                        data,
                                                        norm='log',
                                                        highlight_strong_signal = True,
                                                        plot_SNR=True,
                                                        save_plots=save_name
                                                        )
                    
                    best_P = data['SNR']['best_P']

                    print(f'    Best period found: {best_P:.4f} days\n')

                    all_results.append({
                                        'DSno'        : DSno,
                                        'series'      : proxy,         
                                        'ptype'       : ptype,
                                        'mode'        : mode,
                                        'best_P_SNR'  : best_P
                                        })
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(f'bestP_results_DS{DSno}_all_proxies.csv', index=False, sep=',')
        