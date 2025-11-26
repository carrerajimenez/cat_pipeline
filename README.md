# Ca II Triplet EW Pipeline

**A flexible Python tool for measuring the Equivalent Width (EW) of the Calcium II Triplet (CaT) lines (~8500 Å) in stellar spectra.**

This pipeline automates the reduction and analysis process, including radial velocity correction, robust continuum normalization, and line profile fitting, all controlled via a single YAML configuration file.

## 🚀 Key Features

* **Centralized Management:** Reads a single master list (CSV or FITS table) to manage batch analysis.
* **Flexible Input Formats:** Supports spectra in standard FITS format (WCS header) or table-based formats (CSV, binary FITS tables, ASCII).
* **YAML Configuration:** All pipeline parameters are controlled via a `config.yaml` file; no code modification is required.
* **Robust Pre-Normalization:** Includes an optional module based on the iterative IDL algorithm `continuum.pro` to normalize "raw" spectra before analysis.
* **RV Correction:** Performs cross-correlation with a template to automatically correct Radial Velocity (optional).
* **SNR Handling:** Uses pre-calculated SNR from the input table or estimates it in-situ using the robust `der_snr` algorithm.
* **Modular Fitting:**
    * Fits each CaT line in isolated regions for maximum precision.
    * Supports multiple profiles: `cole` (Gaussian+Lorentzian sum), `gaussian`, `rutledge`.
    * Allows selection of minimization engine: `emcee` (Robust MCMC) or `leastsq` (Fast Levenberg-Marquardt).
* **Detailed Output:** Generates multi-panel diagnostic plots (PDF) and a final results table (CSV).

## 📂 Project Structure

The pipeline expects the following directory structure to run:

```text
/your_project_directory/
├── spectra/                 # Your spectrum files (.fits, .csv, etc.)
├── rv_template/             # Template for RV correction
│   └── cat_template_rv.fits
├── cat_pipeline.py          # Main script
├── line_fitters.py          # Fitting logic and profile definitions
├── rv_tools.py              # Loading and RV tools
├── pre_normalizer.py        # Pre-normalization module
├── config.yaml              # Configuration file
├── master_list.csv          # Table with the list of stars
├── requirements.txt         # Python dependencies
└── README.md

## 🛠️ Installation

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

## ▶️ Usage

The pipeline is designed to be run from the command line, controlled entirely by the configuration file.

1.  **Prepare your data:** Ensure your spectra files are in the `spectra/` directory and your master list (CSV) is available.
2.  **Configure:** Edit the `config.yaml` file to match your directory structure and scientific requirements.
3.  **Run:** Execute the main script, passing the configuration file as an argument.

```bash
python cat_pipeline.py -c config.yaml
```

## ⚙️ Configuration (config.yaml)
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

## 📊 Output Description
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

## 🧩 How to Add a New Profile

This pipeline is designed to be easily extended. To add a new analytical line profile (e.g., "MyVoigt"), follow these steps in `line_fitters.py`:

1. Define the Profile Function: Write a function that calculates the line shape based on input parameters.

```python
def my_voigt_profile(x, center, amplitude, sigma, gamma):
    # ... math here ...
    return calculated_flux
```

2. Define the Residual Function: Write a function that returns `data - model`.

```python
def my_voigt_residuals(params, x, data):
    p = params.valuesdict()
    model = my_voigt_profile(x, p['center'], p['amplitude'], p['sigma'], p['gamma'])
    return data - model
```

3. Define Parameter Creator: Write a function to initialize `lmfit`. Parameters with reasonable bounds.

```python
def create_params_my_voigt(center_guess, region, flux_slice):
    params = Parameters()
    params.add('center', value=center_guess, min=region[0], max=region[1])
    params.add('amplitude', value=0.5, min=0)
    # ... add other parameters ...
    return params
```
4. Define EW Calculator: Write a function to calculate the EW and its error. You can use `idl_tabulate` for numerical integration.
```python
def calculate_ew_my_voigt(result_params, x, snr):
    p = result_params.valuesdict()
    # Integrate model
    eqw = idl_tabulate(x, my_voigt_profile(x, p['center'], p['amplitude'], ...))
    # Estimate error
    center_err = result_params['center'].stderr
    error = np.sqrt(3.55 * np.abs(center_err)) / snr
    return eqw, error
```
## 📏 How to Add a New Region Prescription

The pipeline comes with standard region definitions (e.g., `cenarro2001`), but you can easily add your own prescriptions from other papers or for specific science cases.

You have two options:

### Option 1: Ad-hoc (Directly in Config)
For a quick test, you can define the dictionary directly in your `config.yaml` file without touching the code:

```yaml
# In config.yaml
region_prescription:
  line_regions:
    v1: [8490.0, 8510.0]
    v2: [8530.0, 8550.0]
    v3: [8650.0, 8670.0]
  cont_bands:
    blue: [8470.0, 8560.0] # Flat list of edges is NOT supported here, use the Python method below for complex bands
    red:  [8490.0, 8580.0]

### Option 2: Permanent (Add to Code)

To add a named prescription that you can reuse (like `mypaper2024`), edit `cat_pipeline.py`.

1. Open `cat_pipeline.py`.

2. Locate the `REGION_PRESCRIPTIONS` dictionary inside the `CaTPipeline` class.

3. Add your new definition using the same structure:

```Python
REGION_PRESCRIPTIONS = {
        'cenarro2001': { ... },
        
        # --- YOUR NEW PRESCRIPTION ---
        'mypaper2024': {
            'line_regions': {
                'v1': [8490.0, 8510.0], # [Start, End] for 8498 line
                'v2': [8530.0, 8550.0], # [Start, End] for 8542 line
                'v3': [8650.0, 8670.0]  # [Start, End] for 8662 line
            },
            'cont_bands': {
                # List of start points for continuum bands
                'blue': [8470.0, 8560.0, 8600.0], 
                # List of end points for continuum bands (must match 'blue' length)
                'red':  [8480.0, 8570.0, 8610.0]
            }
        },
    }

```

4. Use it: In `config.yaml`, set region_prescription: `"mypaper2024"`.

## 📝 Credits
Author: R. Carrera (INAF-OAS); 

Contributions: M. Navabi (University of Surrey)

## 📄 License
This project is licensed under the MIT License. See the LICENSE file for details.