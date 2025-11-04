import pandas as pd
import LS_like_periodograms as aux

df = pd.read_csv('DS1_timeSeries_b.csv')
print(df)

time = df.time

RV = df.rv
RV_error = df.e_rv


data = aux.stacked_periodogram(time, RV, RV_error, N_min=50,
                               periodogram_type='BGLS', p_min=1,
                               p_max=5., num_periods = 5000,
                               mode='random')
                               
aux.plot_stacked_periodogram_heatmap(data, norm='linear', highlight_strong_signal = True, plot_SNR=True)
                               
print(data['SNR']['best_P'])
