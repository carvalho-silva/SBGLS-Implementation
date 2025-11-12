import pandas as pd
import numpy as np

import LS_like_periodograms as ls
import remove_planet as rp

import os

os.makedirs("output", exist_ok=True)

planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']


def get_bestP_with_stacked(time, rv, rv_err,
                           DSno,
                           label,
                           periodogram_type='GLS',
                           mode='random without replacement',
                           N_min=55,
                           p_min=1.0,
                           p_max=60.0,
                           num_periods=5000):
    results = ls.stacked_periodogram(
        time, rv, rv_err,
        N_min=N_min,
        periodogram_type=periodogram_type,
        p_min=p_min,
        p_max=p_max,
        num_periods=num_periods,
        mode=mode
    )

    ls.plot_stacked_periodogram_heatmap(
        results,
        norm='log',
        highlight_strong_signal=True,
        plot_SNR=True,
        save_plots=f'output/{periodogram_type}_DS{DSno}_{label}_{mode}'
    )

    best_P = results['SNR']['best_P']
    return best_P, results


def iterative_remove_planets_for_DS(DSno,
                                    max_planets=3,
                                    periodogram_type='GLS',
                                    mode='random without replacement',
                                    N_min=55,
                                    p_min=1.0,
                                    p_max=60.0,
                                    num_periods=5000):
    DSno_str = str(DSno)

    df = pd.read_csv(f'data-sets/DS{DSno_str}_timeSeries.csv')

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

    renamed_df = df.rename(columns=col_dict)

    current_df = renamed_df.copy()
    current_df['rv_orig'] = current_df['rv']

    if 'e_rv' in current_df.columns:
        rv_err = current_df['e_rv'].values
    else:
        rv_err = np.ones_like(current_df['rv'].values)

    summary = []

    for i in range(max_planets):
        planet = planets[i] if i < len(planets) else f'p{i+1}'
        print(f'\n>>> DS{DSno_str}: planet {planet}')

        time = current_df['time'].values
        rv   = current_df['rv'].values

        best_P, stacked_results = get_bestP_with_stacked(
            time, rv, rv_err,
            DSno=DSno_str,
            label=planet,
            periodogram_type=periodogram_type,
            mode=mode,
            N_min=N_min,
            p_min=p_min,
            p_max=p_max,
            num_periods=num_periods
        )

        print(f'  best period (SNR) = {best_P:.5f} days')

		
        residuals, params_fit = rp.remove_planet(
            current_df.time,
            current_df.rv,
            P_init=best_P,
            planet=planet,
            save_plots=f'DS{DSno_str}'  
        )

        summary.append({
            'DSno'   : DSno_str,
            'planet' : planet,
            'best_P' : best_P
        })


        current_df['residuals'] = residuals
        current_df['rv'] = residuals

    df_final = current_df
    return df_final, pd.DataFrame(summary)
    


if __name__ == "__main__":
    DSno_list = [str(i) for i in range(1, 10)]

    all_results = []

    for DSno in DSno_list:
        df_final, summary = iterative_remove_planets_for_DS(
            DSno,
            max_planets=2,             
            periodogram_type='BGLS',  
            mode='random without replacement',
            N_min=55,
            p_min=2.0,
            p_max=90.0,
            num_periods=5000
        )

        summary.to_csv(f'output/periods_DS{DSno}.csv', index=False)
        all_results.append(summary)

    if len(all_results) > 0:
        all_summary = pd.concat(all_results, ignore_index=True)
        all_summary.to_csv('output/periods_all_DS.csv', index=False)

