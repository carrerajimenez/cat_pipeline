"""
File:         line_fitters.py
Description:  A module for fitting spectral line profiles.
              - A global continuum is fit first using pre-defined bands.
              - Fits each of the three CaT lines in separate regions
                on the continuum-normalized spectrum using a
                generic, user-defined profile.
Author(s):    R. Carrera (INAF-OAS)
Contact:      jimenez.carrera@inaf.it
Version:      1.0.0
Date:         18-Nov-2025
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from lmfit import Parameters, minimize
from astropy.stats import sigma_clipped_stats, sigma_clip
from scipy.optimize import curve_fit
import warnings
from uncertainties import ufloat
from numpy import pi, sqrt, exp
import scipy.integrate  # Added import

log = logging.getLogger(__name__)

def cole_profile(x, center, gauss_amp, gauss_sigma, lorentz_amp, lorentz_sigma):
    gauss_sigma = max(gauss_sigma, 1e-10)
    lorentz_sigma = max(lorentz_sigma, 1e-10)
    gauss = gauss_amp * np.exp(-(x - center)**2 / (2 * gauss_sigma**2))
    lorentz = lorentz_amp / (1 + ((x - center) / lorentz_sigma)**2)
    return gauss + lorentz
def gaussian_profile(x, center, gauss_amp, gauss_sigma):
    gauss_sigma = max(gauss_sigma, 1e-10)
    return gauss_amp * exp(-(x - center)**2 / (2 * gauss_sigma**2))
def rutledge_profile(x, center, rutledge_amp, rutledge_sigma):
    rutledge_sigma = max(rutledge_sigma, 1e-10)
    return rutledge_amp / ((1 + ((x - center) / rutledge_sigma)**2))**2.5
def cole_residuals(params, x, data):
    p = params.valuesdict()
    model = cole_profile(x, p['center'], p['gauss_amp'], p['gauss_sigma'],
                         p['lorentz_amp'], p['lorentz_sigma'])
    return data - model
def gaussian_residuals(params, x, data):
    p = params.valuesdict()
    model = gaussian_profile(x, p['center'], p['gauss_amp'], p['gauss_sigma'])
    return data - model
def rutledge_residuals(params, x, data):
    p = params.valuesdict()
    model = rutledge_profile(x, p['center'], p['rutledge_amp'], p['rutledge_sigma'])
    return data - model 
def create_params_cole(center_guess, region, flux_slice):
    params = Parameters()
    params.add('center', value=center_guess, min=center_guess-2, max=center_guess+2)
    params.add('gauss_amp', value=0.5, min=0.001, max=1.2) 
    params.add('gauss_sigma', value=0.5, min=0.001, max=12.0)
    params.add('lorentz_amp', value=0.5, min=0.005, max=1.2)
    params.add('lorentz_sigma', value=2, min=0.001, max=7.0)
    return params
def create_params_gaussian(center_guess, region, flux_slice):
    params = Parameters()
    params.add('center', value=center_guess, min=region[0], max=region[1])
    params.add('gauss_amp', value=np.nanmax(flux_slice), min=0, max=1.5) 
    params.add('gauss_sigma', value=0.5, min=0.01, max=2.0)
    return params
def create_params_rutledge(center_guess, region, flux_slice):
    params = Parameters()
    params.add('center', value=center_guess, min=region[0], max=region[1])
    params.add('rutledge_amp', value=np.nanmax(flux_slice), min=0, max=1.5) 
    params.add('rutledge_sigma', value=0.5, min=0.01, max=2.0)
    return params   
def idl_tabulate(x, y, step=5):
    """Integrate using Newton-Cotes rule, similar to IDL tabulate."""
    def newton_cotes_integrate(x_segment, y_segment):
        if x_segment.shape[0] < 2:
            return 0
        rn = (x_segment.shape[0] - 1) * (x_segment - x_segment[0]) / (x_segment[-1] - x_segment[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x_segment[-1] - x_segment[0]) / (x_segment.shape[0] - 1) * np.dot(weights, y_segment)
    result = 0
    for idx in range(0, x.shape[0], step - 1):
        result += newton_cotes_integrate(x[idx:idx + step], y[idx:idx + step])
    return result
def make_ufloat(param):
    """Helper to safely create ufloat, handling None/0.0 stderr."""
    if param.stderr is not None and param.stderr > 0:
        return ufloat(param.value, param.stderr)
    return ufloat(param.value, np.nan)
def calculate_ew_cole(result_params, x, snr):
    """Calculates EW and error from 'cole' fit parameters."""
    try:
        center_val = result_params['center'].value
        g_amp_val = result_params['gauss_amp'].value
        g_sig_val = result_params['gauss_sigma'].value
        l_amp_val = result_params['lorentz_amp'].value
        l_sig_val = result_params['lorentz_sigma'].value
        eqw = idl_tabulate(x, cole_profile(x, center_val, g_amp_val, g_sig_val, l_amp_val, l_sig_val))
        center_err = result_params['center'].stderr
        if center_err is None: center_err = np.nan
        error = np.sqrt(3.55 * np.abs(center_err)) / snr
        return eqw, error
    except Exception as e:
        log.warning(f"Could not calculate EW for cole profile: {e}")
        return np.nan, np.nan
def calculate_ew_gaussian(result_params, x, snr):
    """Calculates EW and error from 'gaussian' fit parameters."""
    try:
        center_val = result_params['center'].value
        g_amp_val = result_params['gauss_amp'].value
        g_sig_val = result_params['gauss_sigma'].value
        eqw = idl_tabulate(x, gaussian_profile(x, center_val, g_amp_val, g_sig_val))
        center_err = result_params['center'].stderr
        if center_err is None: center_err = np.nan
        error = np.sqrt(3.55 * np.abs(center_err)) / snr
        return eqw, error
    except Exception as e:
        log.warning(f"Could not calculate EW for gaussian profile: {e}")
        return np.nan, np.nan
def calculate_ew_rutledge(result_params, x, snr):
    """Calculates EW and error from 'rutledge' fit parameters."""
    try:
        center_val = result_params['center'].value
        r_amp_val = result_params['rutledge_amp'].value
        r_sig_val = result_params['rutledge_sigma'].value
        eqw = idl_tabulate(x, rutledge_profile(x, center_val, r_amp_val, r_sig_val))
        center_err = result_params['center'].stderr
        if center_err is None: center_err = np.nan
        error = np.sqrt(3.55 * np.abs(center_err)) / snr
        return eqw, error
    except Exception as e:
        log.warning(f"Could not calculate EW for rutledge profile: {e}")
        return np.nan, np.nan
PROFILE_REGISTRY = {
    'cole': {
        'profile_func': cole_profile,
        'residual_func': cole_residuals,
        'param_creator': create_params_cole,
        'ew_calculator': calculate_ew_cole,
    },
    'gaussian': {
        'profile_func': gaussian_profile,
        'residual_func': gaussian_residuals,
        'param_creator': create_params_gaussian,
        'ew_calculator': calculate_ew_gaussian,
    },
    'rutledge': { 
        'profile_func': rutledge_profile,
        'residual_func': rutledge_residuals,
        'param_creator': create_params_rutledge,
        'ew_calculator': calculate_ew_rutledge,
    }
}
def _poly_fit_func(x, a, b, c):
    """A 2nd-degree polynomial for continuum fitting."""
    #return a * (x**2) + b * x + c
    return a * x + b
def fit_continuum(wavelength, flux, blue_edges, red_edges, force_fit=False):
    """
    Fits a global continuum to pre-defined bands.
    If flux is already ~1.0 AND force_fit=False, assumes pre-normalized
    and returns 1.0.
    """
    flux_median = np.nanmedian(flux)
    
    if np.isfinite(flux_median) and abs(flux_median - 1.0) < 0.1 and not force_fit:
        log.warning("Flux appears pre-normalized (median ~ 1.0) and force_fit=False. "
                 "Skipping global continuum fit.")
        return np.ones_like(wavelength)
    elif np.isfinite(flux_median) and abs(flux_median - 1.0) < 0.1 and force_fit:
         log.warning("Flux appears pre-normalized, but force_local_continuum=True. "
                  "Proceeding with fit.")
    
    try:
        n_points = len(blue_edges) #if red_edges[0] < np.min(wavelength) or blue_edges[-1] > np.max(wavelength) else len(blue_edges) - 1
        x_center = np.full(n_points, np.nan, dtype='f4')
        y_percentile = np.full(n_points, np.nan, dtype='f4')
        y_median = np.full(len(blue_edges), np.nan, dtype='f4')
        y_error = np.full(len(blue_edges), np.nan, dtype='f4')
        for i in range(n_points):
            if red_edges[i] < np.min(wavelength) or blue_edges[i] > np.max(wavelength):
                continue
            idx = np.where((wavelength >= blue_edges[i]) & (wavelength <= red_edges[i]))[0]
            if len(idx) > 0:
                x_center[i] = (red_edges[i] - blue_edges[i]) / 2 + blue_edges[i]
                median, mean, std = sigma_clipped_stats(flux[idx], sigma_lower=3, sigma_upper=5)
                y_median[i] = median
                y_error[i] = std
                clipped_flux = sigma_clip(flux[idx], sigma_lower=3, sigma_upper=5, masked=False)
                y_percentile[i] = np.percentile(clipped_flux, 95,method='median_unbiased') if len(clipped_flux) > 1 else median
            log.debug(f"Continuum fitting {blue_edges[i]:.2f} - {red_edges[i]:.2f} nm: Median {y_median[i]:.2f}, 95 percentile {np.percentile(clipped_flux, 95,method='median_unbiased'):.2f} 5 percentile {np.percentile(clipped_flux, 5,method='median_unbiased'):.2f}")
        valid_mask = np.isfinite(x_center) & np.isfinite(y_percentile)
        if np.sum(valid_mask) < 4:
             log.warning("Not enough valid continuum points (< 4) to fit parabola. "
                         "Falling back to median.")
             return np.full_like(wavelength, np.nanmedian(flux))
        x_data = x_center[valid_mask]
        y_data = y_median[valid_mask]
        x_mean = np.mean(x_data)
        x_scaled = x_data - x_mean
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=np.RankWarning)
            warnings.simplefilter("ignore", category=Warning)
            try:
                #coeffs, _ = curve_fit(_poly_fit_func, x_scaled, y_data)
                coeffs, _ = curve_fit(_poly_fit_func, x_center[:n_points], y_data[:n_points])
            except Exception as e:
                log.error(f"curve_fit for continuum failed: {e}. Falling back to median.")
                return np.full_like(wavelength, np.nanmedian(flux))
        continuum_fit = _poly_fit_func(wavelength, *coeffs)
        return continuum_fit
    except Exception as e:
        log.error(f"Global continuum fit failed: {e}. Returning median.")
        return np.full_like(wavelength, np.nanmedian(flux))


class CaTripletFitter: 
    """
    Fits the Ca II Triplet (3 lines) using a *generic* profile
    on a pre-normalized spectrum. Uses lmfit.minimize and emcee.
    """
    LINE_CENTERS = {'v1': 8498.02, 'v2': 8542.09, 'v3': 8662.14} 
    EW_KEYS = {
        'v1': ('EW_8498', 'Err_8498', 'redchi_8498', 'success_8498'),
        'v2': ('EW_8542', 'Err_8542', 'redchi_8542', 'success_8542'),
        'v3': ('EW_8662', 'Err_8662', 'redchi_8662', 'success_8662')
    }
    
    # --- MODIFICADO: __init__ acepta nuevos argumentos ---
    def __init__(self, wave: np.ndarray, flux: np.ndarray, snr: float,
                 line_regions: dict, cont_bands: dict, 
                 profile_name: str = 'cole',
                 fitter_method: str = 'emcee',
                 fitter_prefit_method: str = 'leastsq', # <-- NUEVO ARGUMENTO
                 force_local_continuum: bool = False):
        self.wave = wave
        self.flux = flux
        self.snr = snr
        self.results = {}
        self.global_continuum = None
        self.flux_normalized = None
        self.REGIONS = line_regions
        self.CONT_BLUE_EDGES = cont_bands['blue']
        self.CONT_RED_EDGES = cont_bands['red']
        self.wave_slices = {} 
        
        if profile_name not in PROFILE_REGISTRY:
            log.error(f"Profile '{profile_name}' not found in PROFILE_REGISTRY. "
                      f"Available: {list(PROFILE_REGISTRY.keys())}")
            raise ValueError(f"Profile '{profile_name}' not recognized.")
        
        log.info(f"Initializing fitter with '{profile_name}' profile.")
        self.profile_name = profile_name
        self.fitter_method = fitter_method
        self.fitter_prefit_method = fitter_prefit_method # <-- ALMACENA EL MÉTODO DE PRE-AJUSTE
        self.force_local_continuum = force_local_continuum
        
        funcs = PROFILE_REGISTRY[profile_name]
        self.profile_func = funcs['profile_func']
        self.residual_func = funcs['residual_func']
        self.param_creator_func = funcs['param_creator']
        self.ew_calculator_func = funcs['ew_calculator']

    def fit(self):
        """ 
        First, fits a global continuum.
        Then, fits each line on the normalized spectrum (1-flux)
        using minimize() and the *chosen fitter method*.
        """
        self.results = {}
        self.wave_slices = {}
        self.global_continuum = fit_continuum(
            self.wave, self.flux, self.CONT_BLUE_EDGES, self.CONT_RED_EDGES,
            force_fit=self.force_local_continuum
        )
        if self.global_continuum is None:
            return

        self.flux_normalized = self.flux / self.global_continuum
        log.debug(f"Normalized flux min {self.flux_normalized.min():.2f}, max {self.flux_normalized.max():.2f}")
        
        for prefix, region in self.REGIONS.items():
            center_guess = self.LINE_CENTERS[prefix]
            
            mask = (self.wave >= region[0]) & (self.wave <= region[1])
            if not np.any(mask) or len(self.wave[mask]) < 5:
                self.results[prefix] = None
                continue
                
            wave_slice = self.wave[mask]
            flux_norm_slice = 1.0 - self.flux_normalized[mask]
            self.wave_slices[prefix] = wave_slice 
            
            params = self.param_creator_func(center_guess, region, flux_norm_slice)
            
            try:
                if self.fitter_method == 'emcee':
                    # 1. Run the PRE-FIT first (leastsq or nelder)
                    log.info(f"Running pre-fit method '{self.fitter_prefit_method}' for {prefix}...")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        results_prefit = minimize(self.residual_func, params, 
                                                   args=(wave_slice, flux_norm_slice),
                                                   method=self.fitter_prefit_method, 
                                                   nan_policy="omit")
                    
                    # 2. Run emcee using pre-fit results
                    emcee_params = results_prefit.params.copy()
                    emcee_params.add('__lnsigma', value=np.log(0.01), min=np.log(0.001), max=np.log(0.1))

                    log.info(f"Running emcee for {prefix}...")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore") 
                        self.results[prefix] = minimize(self.residual_func, params=emcee_params, 
                                                        method='emcee', args=(wave_slice, flux_norm_slice),
                                                        steps=1000, burn=500, thin=10, 
                                                        is_weighted=False, progress=False)
                    log.info(f"emcee for {prefix} finished.")
                
                else:
                    # Use any other method (e.g., 'leastsq', 'nelder')
                    log.info(f"Running {self.fitter_method} for {prefix}...")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        self.results[prefix] = minimize(self.residual_func, params, 
                                                       args=(wave_slice, flux_norm_slice),
                                                       method=self.fitter_method,
                                                       nan_policy="omit")
                    log.info(f"{self.fitter_method} for {prefix} finished.")
            
            except Exception as e:
                log.error(f"CaT ({self.profile_name}) fitting with method "
                          f"'{self.fitter_method}' failed for {prefix}: {e}")
                self.results[prefix] = None

    def get_results(self) -> dict:
        """ 
        Extracts EWs, errors, and fit statistics.
        - If method is 'emcee', use the 'ew_calculator_func'.
        - If method is 'leastsq'/'nelder', use continuum perturbation.
        """
        ew_data = {}
        ews = []
        errs = []

        for prefix in self.REGIONS.keys():
            result = self.results.get(prefix)
            ew_key, err_key, redchi_key, success_key = self.EW_KEYS[prefix]
            
            ew_val, ew_err = np.nan, np.nan
            redchi_val, success_val = np.nan, False
            
            wave_slice = self.wave_slices.get(prefix, None)
            
            is_successful = (result and (result.success or self.fitter_method == 'emcee') 
                             and wave_slice is not None)

            if is_successful:
                success_val = result.success
                redchi_val = result.redchi if np.isfinite(result.redchi) else np.nan
                
                ew_val, ew_err_default = self.ew_calculator_func(result.params, wave_slice, self.snr)

                if self.fitter_method == 'emcee':
                    ew_err = ew_err_default
                    log.debug(f"[{prefix}] Using emcee error method: {ew_err:.4f}")
                    
                else:
                    log.debug(f"[{prefix}] Using continuum perturbation for error...")
                    try:
                        mask = (self.wave >= self.REGIONS[prefix][0]) & (self.wave <= self.REGIONS[prefix][1])
                        flux_norm_slice_inverted = 1.0 - self.flux_normalized[mask]
                        params_orig = result.params.copy()
                        
                        flux_high = flux_norm_slice_inverted + (1.0 / self.snr)
                        flux_low  = flux_norm_slice_inverted - (1.0 / self.snr)

                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            res_high = minimize(self.residual_func, params_orig, 
                                                args=(wave_slice, flux_high),
                                                method=self.fitter_method, nan_policy="omit")
                        
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            res_low = minimize(self.residual_func, params_orig, 
                                               args=(wave_slice, flux_low),
                                               method=self.fitter_method, nan_policy="omit")

                        if not res_high.success or not res_low.success:
                            raise Exception("Re-fit de perturbación fallido")

                        ew_high, _ = self.ew_calculator_func(res_high.params, wave_slice, self.snr)
                        ew_low, _  = self.ew_calculator_func(res_low.params, wave_slice, self.snr)
                        
                        ew_err = max(abs(ew_high - ew_val), abs(ew_low - ew_val))
                        log.debug(f"[{prefix}] Perturbation errors: +{ew_high-ew_val:.4f} / {ew_low-ew_val:.4f}. Final error: {ew_err:.4f}")
                        
                    except Exception as e:
                        log.warning(f"[{prefix}] Perturbación del continuo para error falló: {e}. Usando NaN.")
                        ew_err = np.nan
            
            ew_data[ew_key] = ew_val
            ew_data[err_key] = ew_err
            ew_data[redchi_key] = redchi_val
            ew_data[success_key] = success_val
            
            ews.append(ew_val)
            errs.append(ew_err)

        ew_data['Sum_EW'] = np.nansum(ews)
        ew_data['Sum_Err'] = np.sqrt(np.nansum(np.square(errs)))
        if not np.isfinite(ew_data['Sum_EW']): ew_data['Sum_EW'] = np.nan
        if not np.isfinite(ew_data['Sum_Err']): ew_data['Sum_Err'] = np.nan

        return ew_data

    def plot(self, star_id: str, ew_data: dict, rv_kms: float) -> plt.Figure:
        """ 
        Generates a 4-panel diagnostic plot using the NORMALIZED flux.
        """
        fig = plt.figure(figsize=(15, 10))
        gs = GridSpec(2, 3, figure=fig)
        ax_main = fig.add_subplot(gs[0, :])
        ax1 = fig.add_subplot(gs[1, 0])
        ax2 = fig.add_subplot(gs[1, 1])
        ax3 = fig.add_subplot(gs[1, 2])
        
        rv_str = f", RV = {rv_kms:.2f} km/s" if np.isfinite(rv_kms) else ""
        title = f"Ca II Triplet Fit ({self.profile_name}, Global Cont.) for {star_id} (SNR: {self.snr:.1f}){rv_str}"
        fig.suptitle(title, fontsize=16)
        
        wave_nm = self.wave / 10.0
        
        if self.flux_normalized is None:
            log.warning("flux_normalized is None. Cannot plot.")
            ax_main.set_title("ERROR: flux_normalized is None. Fit may have failed.")
            return fig
            
        ax_main.plot(wave_nm, self.flux_normalized, 'bo', markersize=3, label='Data (Normalized)')
        ax_main.axhline(1.0, color='k', ls='--', label='Continuum')

        for i, (blue, red) in enumerate(zip(self.CONT_BLUE_EDGES, self.CONT_RED_EDGES)):
            ax_main.axvspan(blue / 10.0, red / 10.0, color='gray', alpha=0.2, 
                            label='Cont. Band' if i == 0 else None)
        for i, (prefix, region) in enumerate(self.REGIONS.items()):
            ax_main.axvspan(region[0] / 10.0, region[1] / 10.0, color='red', alpha=0.15, 
                            label='Line Region' if i == 0 else None)
        
        for i, (prefix, region) in enumerate(self.REGIONS.items()):
            result = self.results.get(prefix)
            mask = (self.wave >= region[0]) & (self.wave <= region[1])
            
            if result and (result.success or self.fitter_method == 'emcee') and np.any(mask):
                wave_slice_nm = self.wave[mask] / 10.0
                wave_slice = self.wave[mask]
                flux_norm_slice = 1.0 - self.flux_normalized[mask]
                
                residuals = self.residual_func(result.params, wave_slice, flux_norm_slice)
                best_fit_model_norm = flux_norm_slice - residuals
                fit_normalized = 1.0 - best_fit_model_norm 
                
                ax_main.plot(wave_slice_nm, fit_normalized, 'r-', 
                             label="Total Fit" if i == 0 else None)
        
        ax_main.set_ylabel("Normalized Flux")
        ax_main.legend()
        ax_main.set_xlim(846.5, 879.9)
        plot_min = np.nanmin(self.flux_normalized) - 0.1
        ax_main.set_ylim(min(0.1, plot_min), 1.2) 

        # --- Individual Plots ---
        line_axes = [ax1, ax2, ax3]
        for ax, prefix in zip(line_axes, self.REGIONS.keys()):
            region = self.REGIONS[prefix]
            result = self.results.get(prefix)
            ew_key, err_key, redchi_key, success_key = self.EW_KEYS[prefix]
            mask = (self.wave >= region[0]) & (self.wave <= region[1])
            
            if not np.any(mask): continue
            
            wave_slice_nm = self.wave[mask] / 10.0
            flux_norm_slice_region = self.flux_normalized[mask]
            
            ax.plot(wave_slice_nm, flux_norm_slice_region, 'bo', markersize=4, label='Data')
            ax.axhline(1.0, color='k', ls='--', label='Continuum')

            if result and (result.success or self.fitter_method == 'emcee'):
                flux_norm_slice_inverted = 1.0 - flux_norm_slice_region
                residuals = self.residual_func(result.params, self.wave[mask], flux_norm_slice_inverted)
                best_fit_model_norm = flux_norm_slice_inverted - residuals
                fit_normalized = 1.0 - best_fit_model_norm
                
                ax.plot(wave_slice_nm, fit_normalized, 'r-', label='Total Fit')

            ax.set_title(f"Line {self.LINE_CENTERS[prefix]/10.0:.1f} nm")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_xlim(region[0] / 10.0-0.5, region[1] / 10.0+0.5)
            plot_min_region = np.nanmin(flux_norm_slice_region) - 0.1
            ax.set_ylim(min(0.1, plot_min_region), 1.2)

            for (blue, red) in zip(self.CONT_BLUE_EDGES, self.CONT_RED_EDGES):
                ax.axvspan(blue / 10.0, red / 10.0, color='gray', alpha=0.2)

            ew_val = ew_data.get(ew_key, np.nan)
            ew_err = ew_data.get(err_key, np.nan)
            text_str = f"EW = {ew_val:.2f} $\pm$ {ew_err:.2f} $\AA$"
            ax.text(0.05, 0.05, text_str, transform=ax.transAxes, ha='left', 
                    va='bottom', fontsize=10, bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.2'))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        return fig