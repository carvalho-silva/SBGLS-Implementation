import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def sine_func(params, t, P):
	a, b, c = params
	omega = 2 * np.pi / P
	return a * np.sin(omega * t) + b * np.cos(omega * t) + c
	
def fit_func(t, rv, P):
	'''
	Linear least-squares fit of a sine function to the radial velocity data.
	The linear fit can be simplified and written terms of X.
	'''
	omega = 2 * np.pi / P
	t_phase = omega * t

	X = np.column_stack([
		np.sin(t_phase),
		np.cos(t_phase),
		np.ones_like(t_phase)
	])

	coef, _, _, _ = np.linalg.lstsq(X, rv, rcond=None)
	a, b, c = coef
	print('a, b, c =', a, b, c)

	
	rv_fit = sine_func(coef, t, P)
	
	residuals = rv - rv_fit
	
	return coef, rv_fit, residuals, omega

def rv_model(t, per, tp, ecc, om, k):
    """RV Drive
    Args:
        t (array of floats): times of observations (JD)
        per (float): orbital period (days)
        tp (float): time of periastron (JD)
        ecc (float): eccentricity
        om (float): argument of periatron (degree)s
        k (float): radial velocity semi-amplitude (m/s)
        
    Returns:
        rv: (array of floats): radial velocity model
    """

    omega = np.pi * om / 180.
    # Performance boost for circular orbits
    if ecc == 0.0:
        m = 2 * np.pi * (((t - tp) / per) - np.floor((t - tp) / per))
        return k * np.cos(m + omega)

    if per < 0:
        per = 1e-4
    if ecc < 0:
        ecc = 0
    if ecc > 0.99:
        ecc = 0.99

    # Calculate the approximate eccentric anomaly, E1, via the mean anomaly  M.
    nu = true_anomaly(t, tp, per, ecc)
    rv = k * (np.cos(nu + omega) + ecc * np.cos(omega))

    return rv
	
def kepler(Marr, eccarr):
    """Solve Kepler's Equation
    Args:
        Marr (array): input Mean anomaly
        eccarr (array): eccentricity
    Returns:
        array: eccentric anomaly
    """

    conv = 1.0e-12  # convergence criterion
    k = 0.85

    Earr = Marr + np.sign(np.sin(Marr)) * k * eccarr  # first guess at E
    # fiarr should go to zero when converges
    fiarr = ( Earr - eccarr * np.sin(Earr) - Marr)
    convd = np.where(np.abs(fiarr) > conv)[0]  # which indices have not converged
    nd = len(convd)  # number of unconverged elements
    count = 0

    while nd > 0:  # while unconverged elements exist
        count += 1

        M = Marr[convd]  # just the unconverged elements ...
        ecc = eccarr[convd]
        E = Earr[convd]

        fi = fiarr[convd]  # fi = E - e*np.sin(E)-M    ; should go to 0
        fip = 1 - ecc * np.cos(E)  # d/dE(fi) ;i.e.,  fi^(prime)
        fipp = ecc * np.sin(E)  # d/dE(d/dE(fi)) ;i.e.,  fi^(\prime\prime)
        fippp = 1 - fip  # d/dE(d/dE(d/dE(fi))) ;i.e.,  fi^(\prime\prime\prime)

        # first, second, and third order corrections to E
        d1 = -fi / fip
        d2 = -fi / (fip + d1 * fipp / 2.0)
        d3 = -fi / (fip + d2 * fipp / 2.0 + d2 * d2 * fippp / 6.0)
        E = E + d3
        Earr[convd] = E
        fiarr = ( Earr - eccarr * np.sin( Earr ) - Marr) # how well did we do?
        convd = np.abs(fiarr) > conv  # test for convergence
        nd = np.sum(convd is True)

    if Earr.size > 1:
        return Earr
    else:
        return Earr[0]


def true_anomaly(t, tp, per, ecc):
    """
    Calculate the true anomaly for a given time, period, eccentricity.
    Args:
        t (array): array of times in JD
        tp (float): time of periastron, same units as t
        per (float): orbital period in days
        ecc (float): eccentricity
    Returns:
        array: true anomoly at each time
    """

    # f in Murray and Dermott p. 27
    m = 2 * np.pi * (((t - tp) / per) - np.floor((t - tp) / per))
    eccarr = np.zeros(t.size) + ecc
    e1 = kepler(m, eccarr)
    n1 = 1.0 + ecc
    n2 = 1.0 - ecc
    nu = 2.0 * np.arctan((n1 / n2)**0.5 * np.tan(e1 / 2.0))

    return nu
    
def fit_rv_model(t, rv_data, P_fit):

    import numpy as np
    from scipy.optimize import curve_fit

    p0 = [P_fit, t.max(), 0.0, 0.0, np.std(rv_data)]

    params_fit, _ = curve_fit(rv_model,
    t,
    rv_data,
    p0=p0,
    bounds=([0, -np.inf, 0., -360., 0],
    [90., np.inf, 0.5, 360., np.inf]
    )
    )
    P_fit, tp_fit, ecc_fit, om_fit, k_fit = params_fit
    rv_model_fit =  rv_model(t, *params_fit)

    print(f'Período: {P_fit:.3f} dias')
    print(f'tp: {tp_fit:.3f}')
    print(f'e: {ecc_fit:.3f}')
    print(f'\u03C9: {om_fit:.3f}\u00B0')
    print(f'K: {k_fit:.3f} m/s')

    return rv_model_fit, params_fit
    
def remove_planet(t, rv_data, P_init, planet='b', save_plots=False):
    '''
    Fit a Keplerian model to the RV data using rv_model() and returns the residuals and the best-fit parameters.

    Parameters
    ----------
    t : array_like
        Observation times (e.g., in JD or BJD).
    rv_data : array_like
        Observed radial velocities.
    P_init : float
        Initial guess for the orbital period [days].
    planet : str, optional
        Planet label used for naming output files (default: 'b').
    save_plots : bool or str, optional
        If False, show plot. If string, save plot as 'output/<save_plots>_<planet>_phase.png'.

    Returns
    -------
    rv_residuals : ndarray
        Residual RVs after subtracting the fitted Keplerian model.
    params_fit : ndarray
        Best-fit parameters [P, tp, e, omega, K].
    '''
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    os.makedirs('output', exist_ok=True)	

    rv_centered = rv_data - rv_data.median()
    rv_model_centered, params_fit = fit_rv_model(t, rv_centered, P_init)
    rv_model_full = rv_model_centered + rv_data.mean()

    P_fit = params_fit[0]

    t_fit = np.linspace(t.min(), t.max(), 10000)
    rv_model_interp = rv_model(t_fit, *params_fit) + rv_data.mean()

    phase_fit = (t_fit % P_fit) / P_fit

    phase_data = (t % P_fit) / P_fit

    idx = np.argsort(phase_fit)
    phase_sorted = phase_fit[idx]
    rv_model_interp_sorted = rv_model_interp[idx]

    label_params = (
    f'P={params_fit[0]:.3f} d, '
    f'tp={params_fit[1]:.3f}, '
    f'e={params_fit[2]:.3f}, '
    fr'$\omega$={params_fit[3]:.2f}$^\circ$, '
    f'K={params_fit[4]:.3f} m/s'
    )

    plt.figure(figsize=[10,8])
    plt.scatter(phase_data, rv_data, label='Observed', c='k',s=30)
    plt.plot(phase_sorted, rv_model_interp_sorted, label=f'Model ({label_params})', c='tab:cyan', lw=3.)
    plt.xlabel('Orbital Phase', fontsize=14.)
    plt.ylabel('RV m/s', fontsize=14.)
    plt.legend(fontsize=12)
    if save_plots:
        target_name = str(save_plots)
        plt.savefig(f'output/{target_name}_{planet}_phase.png')
        plt.close()
    else:
        plt.show()
    
    return rv_data - rv_model_full, params_fit
