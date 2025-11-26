"""
File:         rv_tools.py
Description:  A module for loading spectra (various formats), performing 
              radial velocity corrections, and other spectrum utilities.
Author(s):    R. Carrera (INAF-OAS)
Contact:      jimenez.carrera@inaf.it
Version:      1.0.0
Date:         18-Nov-2025
"""

import logging
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from pathlib import Path
from numpy.typing import ArrayLike

# --- Constants & Setup ---
CLIGHT = 299792.458  # Speed of light in km/s
log = logging.getLogger(__name__)

# --- SNR Estimation Function ---
def der_snr(flux: np.ndarray) -> float:
    """
    Computes the Signal-to-Noise Ratio (SNR) using the
    DER_SNR algorithm. Robust against continuum trends.
    """
    if flux is None or len(flux) < 5:
        log.warning("Cannot compute der_snr: flux array is too small.")
        return np.nan
    signal = np.nanmedian(flux)
    if len(flux) > 4:
        n = len(flux)
        noise = 1.4826 * np.nanmedian(np.abs(
            2.0 * flux[1:n-1] - flux[0:n-2] - flux[2:n]
        )) / 6.0
    else:
        noise = np.nanstd(flux[1:] - flux[:-1]) / np.sqrt(2.0)
    if noise == 0 or not np.isfinite(noise):
        log.warning("Failed to compute noise for der_snr (noise=0).")
        return np.nan
    return signal / noise

# --- Data Loading Function ---
def load_template_fits(template_path: str | Path) -> (np.ndarray, np.ndarray):
    """
    Loads a FITS template file (e.g., from git_correlation.py's load_template).
    """
    try:
        filepath = str(template_path) + '.fits'
        with fits.open(filepath) as hdu:
            data = hdu[1].data
        log.debug(f"Successfully loaded RV template: {filepath}")
        return data['wavelength'], data['flux']
    except Exception as e:
        log.error(f"Failed to read template FITS {template_path}: {e}")
        return None, None

# --- Class 1: Spectrum Data Holder ---
class Spectrum:
    """
    Handles loading (from various formats), calibration, and 
    slicing of a single spectrum.
    """
    def __init__(self, filepath: Path, spec_config: dict):
        self.filepath = filepath
        self.filename = filepath.name
        self.config = spec_config
        self.wave = None
        self.flux = None
        self.variance = None
        self.header = None
        self.snr = None
        log.debug(f"Spectrum object created for {self.filename}")

    def load(self, name_file_type: str = 'default') -> bool:
        """
        Dispatcher method to load spectrum based on format in config.
        """
        file_format = self.config.get('format', 'fits_primary')
        log.debug(f"Loading spectrum {self.filename} using format: {file_format}")
        try:
            if file_format == 'fits_primary':
                return self._load_fits_primary(name_file_type)
            elif file_format in ['csv', 'fits_bintable', 'ascii']:
                return self._load_table_format()
            else:
                log.error(f"Unknown spectrum format in config: {file_format}")
                return False
        except Exception as e:
            log.error(f"Failed to load file {self.filename}: {e}")
            return False

    def _load_fits_primary(self, name_file_type: str) -> bool:
        """
        Loads data from the primary HDU of a FITS file.
        """
        with fits.open(self.filepath) as hdul:
            self.flux = hdul[0].data
            self.header = hdul[0].header
            
        # Ensure flux is *always* a 1D numpy array of a standard float type
        # .flatten() handles cubes (e.g., [1,1,N]) or 2D images
        # .astype(np.float64) handles non-standard dtypes that SciPy/NumPy might not like
        try:
            self.flux = np.array(self.flux).astype(np.float64).flatten()
            log.debug(f"Flux array successfully flattened to 1D shape: {self.flux.shape}")
        except Exception as e:
            log.error(f"Failed to flatten or cast flux data for {self.filename}: {e}")
            return False
        # --- END OF MODIFICATION ---
            
        self._calibrate_wavelength_from_header()
                    
        log.debug(f"Loaded spectrum data for: {self.filename}")
        return True

    def _load_table_format(self) -> bool:
        col_wave = self.config.get('col_wave')
        col_flux = self.config.get('col_flux')
        col_var = self.config.get('col_var', None)
        if not col_wave or not col_flux:
            log.error(f"Config for table format is missing 'col_wave' or 'col_flux'.")
            return False
        tbl = Table.read(self.filepath)
        if col_wave not in tbl.colnames:
            log.error(f"Column '{col_wave}' (wavelength) not found in {self.filename}")
            return False
        if col_flux not in tbl.colnames:
            log.error(f"Column '{col_flux}' (flux) not found in {self.filename}")
            return False
        self.wave = np.asarray(tbl[col_wave])
        self.flux = np.asarray(tbl[col_flux])
        self.header = tbl.meta
        if col_var and col_var in tbl.colnames:
            self.variance = np.asarray(tbl[col_var])
            log.debug(f"Loaded variance from column '{col_var}'")
        log.debug(f"Loaded spectrum data for: {self.filename}")
        return True

    def _calibrate_wavelength_from_header(self):
        log.debug(f"Calibrating wavelength for {self.filename}...")
        
        # --- MODIFICATION: Ensure NAXIS1 is read from the *correct* header ---
        # If flux was flattened, NAXIS1 might be wrong. Get it from the *data's* shape.
        npix = np.arange(self.flux.shape[-1]) # Use the last dimension
        
        crpix = self.header.get('CRPIX1', 1.0)
        crval = self.header['CRVAL1']
        
        cdelt = self.header.get('CD1_1')
        if cdelt is None:
            cdelt = self.header.get('CDELT1') 
            if cdelt is None:
                log.error(f"Cannot calibrate wavelength for {self.filename}: "
                          f"Header is missing both 'CD1_1' and 'CDELT1'.")
                raise KeyError(f"Header missing 'CD1_1' and 'CDELT1'")
            else:
                log.debug(f"Found 'CDELT1' (value: {cdelt}) for wavelength step.")
        else:
            log.debug(f"Found 'CD1_1' (value: {cdelt}) for wavelength step.")
        
        if crval is None:
            log.error(f"Cannot calibrate wavelength for {self.filename}: "
                      f"Header is missing 'CRVAL1'.")
            raise KeyError(f"Header missing 'CRVAL1'")

        self.wave = (npix - crpix + 1) * cdelt + crval
        if 'DC-FLAG' in self.header and self.header['DC-FLAG'] == 1:
            self.wave = 10**self.wave
            log.debug(f"Applied 10**wave for log-scale: {self.filename}")
        log.debug(f"Wavelength calibration complete for {self.filename} min {self.wave.min():.2f} max {self.wave.max():.2f}")


    def get_region(self, w_min: float, w_max: float) -> (np.ndarray, np.ndarray):
        if self.wave is None:
            log.warning("Spectrum not loaded, cannot get region.")
            return None, None
        mask = (self.wave >= w_min) & (self.wave <= w_max)
        if not np.any(mask):
            log.warning(f"No data found in region {w_min}-{w_max} for {self.filename}")
            return None, None
        log.debug(f"Values min {self.wave[mask].min():.2f} max {self.wave[mask].max():.2f} and flux min {self.flux[mask].min():.2f} max {self.flux[mask].max():.2f} for {self.filename}")
        return self.wave[mask], self.flux[mask]

# --- Class 2: Radial Velocity Corrector ---

class RVCorrector:
    """
    Performs cross-correlation and applies Doppler shift to spectra.
    """
    def __init__(self, template_wave: np.ndarray, template_flux: np.ndarray):
        self.wtemp = template_wave
        self.ftemp = template_flux
        log.debug("RVCorrector initialized with template.")

    @staticmethod
    def _gaussian(x, mu, sigma, amp):
        """Gaussian function for peak fitting."""
        return amp * np.exp(-0.5 * ((x - mu) / sigma)**2)

    def correct(self, wave: np.ndarray, flux: np.ndarray) -> (np.ndarray, np.ndarray, float):
        """
        Corrects the observed flux for radial velocity.
        
        Returns:
            (wave, shifted_flux, rv_kms) tuple
        """
        if wave is None or len(wave) < 50:
            log.warning("RV correction skipped: invalid input wavelength array.")
            return wave, flux, np.nan

        # 1. Create linear-in-log-space grid
        w_min_log = np.log10(wave.min())
        w_max_log = np.log10(wave.max())
        n_points = len(wave)
        log_wave_grid = np.linspace(w_min_log, w_max_log, num=n_points)
        d_log_w = log_wave_grid[1] - log_wave_grid[0]
        linear_wave_grid = 10**log_wave_grid

        # 2. Resample
        obs_interp = interp1d(wave, flux, bounds_error=False, fill_value=np.nanmedian(flux))
        temp_interp = interp1d(self.wtemp, self.ftemp, bounds_error=False, fill_value=np.nanmedian(self.ftemp))
        obs_resampled = obs_interp(linear_wave_grid)
        temp_resampled = temp_interp(linear_wave_grid)

        # 3. Cross-correlate
        corr = np.correlate(1 - obs_resampled, 1 - temp_resampled, mode='full')
        lags = (np.arange(len(corr)) - len(corr) // 2) * d_log_w

        # 4. Find peak
        peak_idx = corr.argmax()
        try:
            fit_window = 20
            start = max(0, peak_idx - fit_window)
            end = min(len(corr), peak_idx + fit_window + 1)
            p0 = [lags[peak_idx], d_log_w * 5, corr[peak_idx]]
            params, _ = curve_fit(self._gaussian, lags[start:end], corr[start:end], p0=p0)
            log_w_shift = params[0]
        except Exception as e:
            log.warning(f"Gaussian peak fit failed in RV correction: {e}. Using peak index.")
            log_w_shift = lags[peak_idx]

        # v/c = 10^(delta_log_w) - 1
        rv_kms = CLIGHT * (10**log_w_shift - 1)
        log.info(f"Calculated RV shift: {rv_kms:.2f} km/s")

        # Apply shift: wave_rest = wave_obs / (1 + v/c) = wave_obs / (10^log_w_shift)
        shifted_wave_grid = wave / (10**log_w_shift)

        # 6. Interpolate back to original grid
        final_flux = np.interp(wave, shifted_wave_grid, flux)
        
        return wave, final_flux, rv_kms