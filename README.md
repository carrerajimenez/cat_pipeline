### 🛠️ Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your_username/cat-pipeline.git](https://github.com/your_username/cat-pipeline.git)
    cd cat-pipeline
    ```

2.  **(Optional) Create a virtual environment:**
    It is recommended to use a virtual environment to avoid conflicts with other system packages.
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Linux/macOS
    # venv\Scripts\activate   # On Windows
    ```

3.  **Install dependencies:**
    Install all required Python packages listed in `requirements.txt`.
    ```bash
    pip install -r requirements.txt
    ```

### ▶️ Usage

The pipeline is designed to be run from the command line, controlled entirely by the configuration file.

1.  **Prepare your data:** Ensure your spectra files are in the `spectra/` directory and your master list (CSV) is available.
2.  **Configure:** Edit the `config.yaml` file to match your directory structure and scientific requirements.
3.  **Run:** Execute the main script, passing the configuration file as an argument.

```bash
python cat_pipeline.py -c config.yaml
```

### ⚙️ Configuration (config.yaml)
The config.yaml file controls every aspect of the pipeline. Below is a description of the available parameters:

```yaml
# --- Core Paths ---
workdir: "."                  # Base path for the project. "." uses the current directory.
spectra_file: "spectra"       # Folder containing spectrum files relative to workdir.
template_file: "templates/rv" # Path to the RV template (without .fits extension).
name_file: "run_results"      # Base name for output files (PDF and CSV).

# --- Input Table ---
input_table: "master.csv"     # The master list of targets.
col_id: "Star_ID"             # Column name for the star ID.
col_filename: "file"          # Column name for the spectrum filename.
spectra_file_ext: ".fits"     # Extension to append to filenames (optional).

# --- Processing Options ---
col_snr: "snr"                # Column for SNR. Set to null if missing.
compute_snr_if_missing: true  # Calculate SNR using der_snr if not provided.
skip_rv_correction: false     # Set true to skip cross-correlation.

# --- Pre-Normalization ---
pre_normalization:
  apply: true                 # Enable iterative continuum fitting.
  order: 2                    # Polynomial order (>=0) or method flag (<0).
  niter: 5                    # Number of iterations or filter width.
  lowrej: 1.0                 # Lower sigma rejection threshold.
  highrej: 3.0                # Upper sigma rejection threshold.

# --- Fitter Settings ---
region_prescription: "cenarro2001" # Passbands definition (continuum/line regions).
fitter_type: "cole"           # Profile model: 'cole', 'gaussian', 'rutledge'.
fitter_method: "leastsq"      # Minimizer: 'emcee', 'leastsq', 'nelder'.
fitter_prefit_method: "leastsq" # Pre-minimizer used before 'emcee'.
force_local_continuum: false  # Force local fit even if pre-normalized.

# --- File Format ---
spectrum_config:
  format: "fits_primary"      # 'fits_primary', 'csv', 'fits_bintable', 'ascii'.
  col_wave: "WAVELENGTH"      # Column name for wavelength (if table format).
  col_flux: "FLUX"            # Column name for flux (if table format).
```

### 📊 Output Description
The pipeline generates a CSV file containing the following columns for each star processed:

| Column | Description |
| :--- | :--- |
| `ID` | Unique identifier of the star from the input table. |
| `RV_kms` | Radial Velocity shift applied (km/s). `NaN` if skipped. |
| `Input_SNR` | SNR value provided in the input table. |
| `EW_8498` | Equivalent Width for the 8498 Å line. |
| `Err_8498` | Estimated error for the 8498 Å line EW. |
| `redchi_8498` | Reduced Chi-Square goodness-of-fit statistic for the 8498 Å line. |
| `success_8498` | Boolean flag indicating if the fit converged for the 8498 Å line. |
| `EW_8542` | Equivalent Width for the 8542 Å line. |
| `Err_8542` | Estimated error for the 8542 Å line EW. |
| `redchi_8542` | Reduced Chi-Square goodness-of-fit statistic for the 8542 Å line. |
| `success_8542` | Boolean flag indicating if the fit converged for the 8542 Å line. |
| `EW_8662` | Equivalent Width for the 8662 Å line. |
| `Err_8662` | Estimated error for the 8662 Å line EW. |
| `redchi_8662` | Reduced Chi-Square goodness-of-fit statistic for the 8662 Å line. |
| `success_8662` | Boolean flag indicating if the fit converged for the 8662 Å line. |
| `Sum_EW` | Sum of the EWs of all three lines. |
| `Sum_Err` | Propagated error for the sum of EWs. |
| `Computed_SNR` | SNR calculated by the pipeline using `der_snr` (if enabled). |
```