"""
File:         pre_normalizer.py
Description:  A module for pre-normalizing spectra using various methods
              based on the IDL continuum.pro logic.
Author(s):    R. Carrera (INAF-OAS)
Contact:      jimenez.carrera@inaf.it
Version:      1.0.0
Date:         05-Nov-2025
"""

import logging
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import ndimage
from astropy.stats import sigma_clipped_stats, sigma_clip

log = logging.getLogger(__name__)

class SpectrumNormalizer:
    """
    Applies spectral normalization using various methods controlled
    by the 'order' parameter, mimicking the IDL continuum.pro script.
    
    order >= 0:   Iterative polynomial fit [cite: 3]
    order == -1:  Global median of the entire spectrum 
    -1 < order < 0: Running percentile filter [cite: 4, 29]
    order == -2:  Piecewise fit (PEM), robust high-point 
    order == -3:  Running average (smooth) 
    -200 < order < -100: Piecewise fit (PEM) with percentile [cite: 6, 29]
    """
    
    def __init__(self, order: float, niter: int, lowrej: float, highrej: float):
        """
        Args:
            order (float): Controls the normalization method.
            niter (int): Controls iterations, filter width, or # of pieces.
            lowrej (float): Sigma rejection threshold for low points (order >= 0).
            highrej (float): Sigma rejection threshold for high points (order >= 0).
        """
        self.order = order
        self.niter = int(niter)
        self.lowrej = lowrej
        self.highrej = highrej
        log.info(f"SpectrumNormalizer initialized: order={order}, "
                 f"niter={niter}, rej={lowrej}/{highrej}")

    def normalize(self, spectrum_obj):
        """
        Fits a continuum and normalizes the spectrum flux in-place.
        This is the main dispatcher method.
        """
        wave = spectrum_obj.wave
        flux = spectrum_obj.flux
        
        continuum = None
        
        if self.order >= 0:
            # Iterative polynomial fit [cite: 3]
            log.debug("Using iterative polynomial normalizer (order >= 0)")
            continuum = self._fit_polynomial(wave, flux)
            
        elif self.order == -1:
            # Global median 
            log.debug("Using global median normalizer (order == -1)")
            continuum = self._fit_global_median(wave, flux)
            
        elif -1 < self.order < 0:
            # Percentile filter [cite: 4, 29]
            log.debug("Using percentile filter normalizer (-1 < order < 0)")
            continuum = self._fit_percentile_filter(wave, flux)
            
        elif self.order == -2:
            # PEM (Piecewise fit), default robust 
            log.debug("Using PEM normalizer (order == -2)")
            continuum = self._fit_pem(wave, flux)
            
        elif self.order == -3:
            # Running average (smooth) 
            log.debug("Using running average normalizer (order == -3)")
            continuum = self._fit_running_average(wave, flux)
            
        elif -200 < self.order < -100:
            # PEM (Piecewise fit) with percentile [cite: 6, 29]
            log.debug("Using PEM (percentile) normalizer (-200 < order < -100)")
            continuum = self._fit_pem(wave, flux)
            
        else:
            log.warning(f"Invalid order '{self.order}'. Skipping normalization.")
            return

        if continuum is not None and np.any(continuum > 0):
            spectrum_obj.flux = spectrum_obj.flux / continuum
            log.info("Pre-normalization complete. Spectrum flux is now normalized.")
        else:
            log.error("Continuum fit returned None or zero. Flux was not normalized.")

    def _fit_polynomial(self, wave, flux):
        """ Iterative polynomial fit (order >= 0) [cite: 3] """
        yy = flux.copy()
        xx = wave.copy()
        mask = np.isfinite(yy)
        
        if np.sum(mask) < self.order + 1:
            log.error("Not enough valid data points to fit continuum.")
            return None

        for i in range(self.niter):
            try:
                coeffs = P.polyfit(xx[mask], yy[mask], int(self.order))
            except np.linalg.LinAlgError as e:
                log.warning(f"Polyfit failed on iteration {i+1}: {e}. Using last valid fit.")
                break 

            yfit = P.polyval(xx, coeffs)
            residuals = yy[mask] - yfit[mask]
            sigma = np.std(residuals)
            
            # Find points to reject [cite: 9, 10]
            low_mask = residuals < -self.lowrej * sigma
            high_mask = residuals > self.highrej * sigma
            
            reject_indices = np.where(mask)[0][low_mask | high_mask]
            
            if reject_indices.size == 0:
                log.debug(f"Continuum converged on iteration {i+1}.")
                break 
            
            yy[reject_indices] = yfit[reject_indices]
            
        try:
            final_coeffs = P.polyfit(xx[mask], yy[mask], int(self.order))
        except np.linalg.LinAlgError as e:
            log.error(f"Final continuum polyfit failed: {e}")
            return None
            
        continuum = P.polyval(wave, final_coeffs)
        return continuum

    def _fit_global_median(self, wave, flux):
        """ Global median (order == -1)  """
        median_val = np.nanmedian(flux)
        if np.isfinite(median_val):
            return np.full_like(wave, median_val)
        else:
            log.error("Could not calculate median of flux.")
            return None

    def _fit_percentile_filter(self, wave, flux):
        """ Running percentile filter (-1 < order < 0) [cite: 4, 29] """
        percentile = -self.order * 100.0  # order=-0.5 -> percentile=50
        width = self.niter # niter is filter width [cite: 7]
        log.debug(f"Running percentile_filter: width={width}, percentile={percentile}")
        log.debug(f"(type of input array {type(flux)} with dimensions {flux.shape})")
        return ndimage.percentile_filter(flux, percentile, size=width, mode='mirror')

    def _fit_running_average(self, wave, flux):
        """ Running average (order == -3)  """
        width = self.niter # niter is filter width [cite: 8]
        log.debug(f"Running uniform_filter1d (smooth): width={width}")
        return ndimage.uniform_filter1d(flux, size=width, mode='mirror')

    def _fit_pem(self, wave, flux):
        """ 
        Piecewise (PEM) fit (order == -2 or -200 < order < -100) [cite: 5, 6, 29]
        This is a Python equivalent, as pem.pro is not a standard IDL func.
        """
        n_pieces = self.niter # niter is number of pieces [cite: 7]
        
        if self.order == -2:
            percentile = 95.0 # Default: fit to high points
            log.debug(f"Running PEM fit: n_pieces={n_pieces}, percentile={percentile}")
        else:
            # order is -100 to -200
            percentile = -(100.0 + self.order) # e.g. order=-150 -> percentile=50 [cite: 6]
            log.debug(f"Running PEM fit: n_pieces={n_pieces}, percentile={percentile}")

        if n_pieces < 2:
            log.error(f"PEM fit requires niter >= 2 (for n_pieces). Got {n_pieces}.")
            return None
            
        # Split arrays into pieces
        wave_chunks = np.array_split(wave, n_pieces)
        flux_chunks = np.array_split(flux, n_pieces)
        
        x_points = []
        y_points = []
        
        for w, f in zip(wave_chunks, flux_chunks):
            if len(f) > 0:
                try:
                    # Find robust high point in the chunk
                    clipped_f = sigma_clip(f, sigma_lower=3, sigma_upper=3, masked=False)
                    y_points.append(np.nanpercentile(clipped_f, percentile))
                    x_points.append(np.nanmedian(w))
                except Exception as e:
                    log.warning(f"Could not find percentile in chunk: {e}")
        
        x_points = np.array(x_points)
        y_points = np.array(y_points)
        
        valid = np.isfinite(x_points) & np.isfinite(y_points)
        if np.sum(valid) < 2:
            log.error("Not enough valid chunks to perform PEM interpolation.")
            return None
        
        # Interpolate between the anchor points
        try:
            # Use linear interpolation, fill edges with nearest value
            interpolator = interp1d(x_points[valid], y_points[valid], 
                                    kind='linear', bounds_error=False, 
                                    fill_value=(y_points[valid][0], y_points[valid][-1]))
            
            return interpolator(wave)
        except Exception as e:
            log.error(f"PEM interpolation failed: {e}")
            return None