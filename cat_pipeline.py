"""
File:         cat_pipeline.py
Description:  The main orchestration script for the Ca II Triplet EW pipeline.
              Reads a single input table, runs analysis, and saves results.
Author(s):    R. Carrera (INAF-OAS)
Contact:      jimenez.carrera@inaf.it
Version:      1.0.0
Date:         18-Nov-2025
"""

import os
import re
import csv
import logging
from pathlib import Path
import time
import argparse
import yaml

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

# --- Local Module Imports ---
from rv_tools import Spectrum, RVCorrector, load_template_fits, der_snr
from line_fitters import CaTripletFitter 
from pre_normalizer import SpectrumNormalizer 

# --- Setup root logger for the script ---
log = logging.getLogger(__name__)


class CaTPipeline:
    """
    Orchestrates the entire process from file discovery to saving results.
    """
    
    REGION_PRESCRIPTIONS = {
        'cenarro2001': {
            'line_regions': {
                'v1': [8484.0, 8513.0],
                'v2': [8522.0, 8562.0],
                'v3': [8642.0, 8682.0]
            },
            'cont_bands': {
                'blue': [8474.0, 8563.0, 8619.0, 8700.0, 8776.0],
                'red':  [8484.0, 8577.0, 8642.0, 8725.0, 8792.0]
            }
        },
        'battaglia2008': {
            'line_regions': {
                'v1': [8498.0-15.0, 8498.0+15.0],
                'v2': [8542.0-15.0, 8542.0+15.0],
                'v3': [8662.0-15.0, 8662.0+15.0] 
            },
            'cont_bands': { 
                'blue': [8474.0, 8563.0, 8619.0, 8700.0, 8776.0],
                'red':  [8484.0, 8577.0, 8642.0, 8725.0, 8792.0]
            }
        },
        'armandroffzinn1988': {
            'line_regions': {
                'v1': [8490.0, 8506.0],
                'v2': [8532.0, 8552.0],
                'v3': [8653.0, 8671.0]
            },
            'cont_bands': {
                'blue': [8474.0, 8521.0, 8555.0, 8626.0, 8695.0],
                'red':  [8489.0, 8531.0, 8595.0, 8650.0, 8725.0]
            }
        },
    }
    
    def __init__(self, config: dict):
        self.config = config
        self.workdir = Path(config['workdir']).resolve()
        log.info(f"Working directory: {self.workdir}")
        
        # Config for single table
        self.input_table_path = self.workdir / config['input_table']
        self.col_id = config['col_id']
        self.col_filename = config['col_filename']
        
        # Spectrum config
        self.spec_config = config.get('spectrum_config', 
                                    {'format': 'fits_primary'})
        self.file_ext = config.get('spectra_file_ext', '.fits')
        if self.file_ext and not self.file_ext.startswith('.'):
            self.file_ext = '.' + self.file_ext

        # SNR config
        self.col_snr = config.get('col_snr', None) 
        self.compute_snr = config.get('compute_snr_if_missing', False)
        
        # RV Config
        self.skip_rv = config.get('skip_rv_correction', False)

        # Region Config
        prescription = config.get('region_prescription', 'cenarro2001')
        if isinstance(prescription, str):
            if prescription not in self.REGION_PRESCRIPTIONS:
                log.warning(f"Unknown region prescription '{prescription}'. "
                            f"Using 'cenarro2001'.")
                self.active_regions = self.REGION_PRESCRIPTIONS['cenarro2001']
            else:
                log.info(f"Using region prescription: '{prescription}'")
                self.active_regions = self.REGION_PRESCRIPTIONS[prescription]
        elif isinstance(prescription, dict):
            log.info("Using custom region prescription from config.")
            self.active_regions = prescription
        else:
            log.error("Invalid 'region_prescription' format. Using 'cenarro2001'.")
            self.active_regions = self.REGION_PRESCRIPTIONS['cenarro2001']

        self.fitter_type = config.get('fitter_type', 'cole').lower() 
        
        self.fitter_method = config.get('fitter_method', 'emcee').lower()
        self.fitter_prefit_method = config.get('fitter_prefit_method', 'leastsq').lower()
        log.info(f"Using fitting engine: {self.fitter_method}")
        if self.fitter_method == 'emcee':
            log.info(f"Using pre-fit '{self.fitter_prefit_method}' para emcee.")
        
        self.force_local_continuum = config.get('force_local_continuum', False)
        if self.force_local_continuum:
            log.info("Forcing local continuum fit (force_local_continuum=True).")

        self.spectra_path = self.workdir / config['spectra_file']
        self.template_path = self.workdir / config['template_file']
        self.name_file = config['name_file'] # For output
        
        # Pre-Normalization
        self.pre_norm_config = config.get("pre_normalization", {"apply": False})
        self.normalizer = None
        if self.pre_norm_config.get("apply", False):
            self.normalizer = SpectrumNormalizer(
                order=self.pre_norm_config.get("order", 2),
                niter=self.pre_norm_config.get("niter", 5),
                lowrej=self.pre_norm_config.get("lowrej", 1.0),
                highrej=self.pre_norm_config.get("highrej", 3.0)
            )
            log.info("Pre-normalization step ENABLED.")
        else:
            log.info("Pre-normalization step DISABLED.")
        
        self.figdir = self.workdir / 'eqw_figures'
        self.outdir = self.workdir / 'eqw_result'
        self.figdir.mkdir(exist_ok=True)
        self.outdir.mkdir(exist_ok=True)
        
        self.jobs = []
        self.rv_corrector = None

    def _load_input_table(self):
        try:
            tbl = Table.read(self.input_table_path)
        except Exception as e:
            log.critical(f"Failed to read input table: {self.input_table_path} - {e}")
            return
        required_cols = [self.col_id, self.col_filename]
        for col in required_cols:
            if col not in tbl.colnames:
                log.critical(f"Input table does not have the required column: '{col}'")
                return
        for row in tbl:
            snr_val = None
            if self.col_snr and self.col_snr in tbl.colnames:
                snr_val = row[self.col_snr]
                if not np.isfinite(snr_val):
                    snr_val = None
            self.jobs.append({
                'id': str(row[self.col_id]),
                'filename': str(row[self.col_filename]),
                'snr': snr_val
            })
        log.info(f"Loaded {len(self.jobs)} jobs from {self.input_table_path}")

    def load_dependencies(self):
        self._load_input_table()
        if not self.skip_rv:
            wtemp, ftemp = load_template_fits(self.template_path)
            if wtemp is not None:
                self.rv_corrector = RVCorrector(wtemp, ftemp)
                log.info("RV Corrector ready.")
            else:
                log.error("Failed to load RV template. RV correction will fail.")
        else:
            log.info("Skipping RV template loading (configuration).")

    def process_spectrum(self, job: dict, i: int, total_jobs: int) -> (dict, plt.Figure):
        """
        Full processing pipeline for a single spectrum job.
        """
        star_id = job['id']
        filename = job['filename']
        has_extension = os.path.splitext(filename)[1] != ""
        if not has_extension and self.file_ext:
            filename = filename + self.file_ext
            
        filepath = self.spectra_path / filename
        
        log.info(f"--- Processing {i} of {total_jobs}: {star_id} ({filename}) ---")

        if not filepath.exists():
            log.warning(f"File not found: {filepath}. Skipping.")
            return None, None
            
        # 1. Load Spectrum
        spectrum = Spectrum(filepath, self.spec_config)
        if not spectrum.load(name_file_type=self.name_file):
            log.warning(f"Failed to load spectrum {filepath}. Skipping.")
            return None, None
            
        # 1b. Pre-Normalization
        if self.normalizer:
            log.info("Applying iterative pre-normalization...")
            try:
                self.normalizer.normalize(spectrum) 
            except Exception as e:
                log.error(f"Pre-normalization failed for {filename}: {e}")
                return None, None
        
        # 2. Handle SNR
        snr = job['snr']
        computed_snr_val = np.nan
        
        if (snr is None) and self.compute_snr:
            log.info(f"SNR for {star_id} not found. Calculating with der_snr...")
            _, flux_reg_snr = spectrum.get_region(8450., 8725.)
            if flux_reg_snr is not None and len(flux_reg_snr) > 10:
                snr = der_snr(flux_reg_snr)
                computed_snr_val = snr
                log.info(f"Calculated SNR: {snr:.2f}")
            else:
                log.warning("Cannot calculate SNR, insufficient data in region.")
                snr = np.nan
        elif snr is None:
            log.warning(f"No SNR for {star_id} and compute_snr=False. Skipping.")
            return None, None
        spectrum.snr = snr
            
        # 3. Get CaT Region
        all_regions = self.active_regions['line_regions'].values()
        all_bands = self.active_regions['cont_bands'].values()
        all_values = [item for sublist in all_regions for item in sublist] + \
                     [item for sublist in all_bands for item in sublist]
        min_wave = min(all_values)
        max_wave = max(all_values)
        log.info(f"Extracting CaT region: {min_wave} - {max_wave} Å")
        
        wave_reg, flux_reg = spectrum.get_region(min_wave, max_wave)
        if wave_reg is None:
            log.warning(f"Insufficient data in CaT region for {filename}.")
            return None, None
            
        # 4. Handle RV Correction
        rv_shift = np.nan
        if not self.skip_rv and self.rv_corrector:
            log.info("Performing RV correction...")
            wave_rv, flux_rv, rv_shift = self.rv_corrector.correct(wave_reg, flux_reg)
        elif self.skip_rv:
            log.info("Skipping RV correction (configuration). Fitting in 'observed frame'.")
            wave_rv, flux_rv = wave_reg, flux_reg
        else:
            log.warning("RV correction requested but corrector failed. "
                        "Fitting in 'observed frame'.")
            wave_rv, flux_rv = wave_reg, flux_reg
        
        # 5. Fit Lines
        try:
            fitter = CaTripletFitter(
                wave=wave_rv, 
                flux=flux_rv, 
                snr=spectrum.snr,
                line_regions=self.active_regions['line_regions'],
                cont_bands=self.active_regions['cont_bands'],
                profile_name=self.fitter_type,
                fitter_method=self.fitter_method,
                fitter_prefit_method=self.fitter_prefit_method,
                force_local_continuum=self.force_local_continuum
            )
        except ValueError as e:
            log.error(f"Failed to initialize fitter: {e}. Skipping star.")
            return None, None
            
        fitter.fit()
        
        # 6. Get Results
        ew_data = fitter.get_results()
        
        # 7. Create Plot
        fig = fitter.plot(star_id, ew_data=ew_data, rv_kms=rv_shift)
        
        # 8. Combine all info
        result_data = {
            'ID': star_id, 
            'RV_kms': rv_shift,
            'Input_SNR': snr,
            'Computed_SNR': computed_snr_val,
            **ew_data 
        }
        
        return result_data, fig

    def save_results(self, all_results: list):
        if not all_results:
            log.warning("No results generated. Output file will be empty.")
            return
        output_path = self.outdir / f"{self.name_file}.csv"
        header_list = [
            'ID', 'RV_kms', 'Input_SNR',
            'EW_8498', 'Err_8498', 'redchi_8498', 'success_8498',
            'EW_8542', 'Err_8542', 'redchi_8542', 'success_8542',
            'EW_8662', 'Err_8662', 'redchi_8662', 'success_8662',
            'Sum_EW', 'Sum_Err', 'Computed_SNR'
        ]
        final_table_data = []
        for row in all_results:
            final_table_data.append({k: row.get(k, np.nan) for k in header_list})
        try:
            result_table = Table(final_table_data)[header_list]
            for col in result_table.colnames:
                if col not in ['ID', 'success_8498', 'success_8542', 'success_8662']:
                    result_table[col].format = '%.4f'
            for col in ['success_8498', 'success_8542', 'success_8662']:
                 if col in result_table.colnames:
                    result_table[col].format = '%s'

            result_table.write(output_path, format='csv', overwrite=True)
            log.info(f"Results saved successfully to {output_path} ({len(all_results)} objects)")
        except Exception as e:
            log.error(f"Failed to write output CSV: {e}")

    def run(self):
        start_time = time.time() 
        
        self.load_dependencies()
        if not self.jobs:
            log.error("Pipeline cannot start: No jobs loaded from input table.")
            return
            
        pdf_path = self.figdir / f'{self.name_file}.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path)
        all_results = []
        
        total_jobs = len(self.jobs)
        log.info(f"Starting processing of {total_jobs} objects...")
        
        for i, job in enumerate(self.jobs, start=1):
            result_data, fig = self.process_spectrum(job, i, total_jobs)
            
            if result_data and fig:
                all_results.append(result_data)
                pdf.savefig(fig)
                plt.close(fig)
        
        pdf.close()
        log.info(f"Processing finished. PDF saved to {pdf_path}")
        
        self.save_results(all_results)
        
        end_time = time.time()
        total_duration_sec = end_time - start_time
        total_minutes = int(total_duration_sec // 60)
        remaining_seconds = total_duration_sec % 60
        
        log.info("--- Pipeline Finished ---")
        log.info(f"Total execution time: {total_minutes} minutes and {remaining_seconds:.2f} seconds.")

# --- Main execution block ---

if __name__ == "__main__":
    """
    This block runs if the script is executed directly.
    Loads configuration from a YAML file.
    """
    
    # --- 0. Configure Logging (Default INFO) ---
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # --- 1. Parse Command-Line Arguments ---
    parser = argparse.ArgumentParser(
        description="Pipeline for measuring EW of the Ca II Triplet."
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        required=True,
        help="Path to the YAML configuration file (e.g., config.yaml)"
    )
    args = parser.parse_args()

    # --- 2. Load Configuration from YAML ---
    config_path = Path(args.config)
    if not config_path.exists():
        log.critical(f"ERROR: Configuration file not found at: {config_path}")
        exit()

    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        log.info(f"Configuration loaded from {config_path}")
        
        # --- NEW: Set logging level from config ---
        log_level_str = config.get('logging_level', 'INFO').upper()
        
        # Validate level
        valid_levels = {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'}
        if log_level_str not in valid_levels:
            log.warning(f"Invalid logging_level '{log_level_str}' in config. Defaulting to INFO.")
            log_level_str = 'INFO'
            
        # Apply new level
        new_level = getattr(logging, log_level_str)
        logging.getLogger().setLevel(new_level)
        log.info(f"Logging level updated to {log_level_str}")

    except Exception as e:
        log.critical(f"ERROR: Could not read or parse YAML file: {e}")
        exit()
    
    # --- 3. Start the Pipeline ---
    log.info(f"Starting pipeline with configuration...")

    try:
        pipeline = CaTPipeline(config)
        pipeline.run()
    except KeyError as e:
        log.critical(f"ERROR: Missing key in configuration file: {e}")
        log.critical("Please check your config.yaml file against the README.")
    except Exception as e:
        log.critical(f"ERROR: Pipeline failed with an unexpected error: {e}")