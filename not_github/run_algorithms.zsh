#!/usr/bin/zsh

# === CONFIGURATION ===
# --- Adjust these variables ---

# Path to your python script
PYTHON_SCRIPT="/media/rcarrera/scratch1/CaTriplet/2025/code_to_publish/cat_pipeline.py"

# --- NEW: Define your list of YAML files here ---
# Add or remove as many files as you want.
YAML_FILES=(
    "/media/rcarrera/scratch1/CaTriplet/2025/fixedcont_spectra/R_8500/config_85k_snr100_nelder.yaml"
    "/media/rcarrera/scratch1/CaTriplet/2025/fixedcont_spectra/R_8500/config_85k_snr100_lsq.yaml"
    "/media/rcarrera/scratch1/CaTriplet/2025/fixedcont_spectra/R_8500/config_85k_snr100_genetic.yaml"
    # You could add more here:
    # "/path/to/another/config.yaml"
)

# Python executable (e.g., python3, or /path/to/your/venv/bin/python)
PYTHON_EXECUTABLE="python"

# Log file for nohup output
LOG_FILE="run_on_algoriths.log"

# === SCRIPT START ===

echo "--- Generic sequential script launcher started ---"

# --- 1. Check if all YAML files exist before starting ---
echo "Checking for input YAML files..."
if [ ${#YAML_FILES[@]} -eq 0 ]; then
    echo "ERROR: The YAML_FILES list is empty. Nothing to run."
    echo "Aborting script."
    exit 1
fi

for yaml_file in "${YAML_FILES[@]}"; do
    if [ ! -f "$yaml_file" ]; then
        echo "ERROR: Required YAML file not found: $yaml_file"
        echo "Aborting script."
        exit 1
    fi
done
echo "All ${#YAML_FILES[@]} YAML files found. Proceeding..."


# --- 2. Run Python scripts sequentially with nohup ---
echo "Starting Python scripts in the background..."

# We pass the variables and the array of files to the 'zsh -c' shell
# $0 = "zsh" (script name for the shell)
# $1 = PYTHON_SCRIPT
# $2 = PYTHON_EXECUTABLE
# $3... = The list of YAML files
nohup zsh -c '
    # Read the arguments passed from the main script
    PYTHON_SCRIPT_NOHUP="$1"
    PYTHON_EXECUTABLE_NOHUP="$2"
    
    # Shift away the first 2 args, so "$@" is just the YAML list
    shift 2 
    
    TOTAL_FILES=$#
    echo "--- Batch process started at $(date) ---"
    echo "Total files to process: $TOTAL_FILES"

    COUNT=0
    # Loop over all remaining arguments (the YAML files)
    for yaml_file in "$@"; do
        COUNT=$((COUNT + 1))
        echo ""
        echo "--- Starting Run $COUNT/$TOTAL_FILES: $yaml_file at $(date) ---"
        
        # 1. Try to run the Python command
        if "$PYTHON_EXECUTABLE_NOHUP" "$PYTHON_SCRIPT_NOHUP" -c "$yaml_file"; then
            echo "--- Run $COUNT/$TOTAL_FILES Succeeded ---"
        else
            # 2. If it failed, log it and break the loop
            echo "--- !!! ERROR: Run $COUNT/$TOTAL_FILES FAILED ($yaml_file) at $(date) !!! ---"
            echo "Aborting remaining runs."
            break # Stop the loop
        fi
    done
    
    echo ""
    echo "--- All runs finished at $(date) ---"
    
' "zsh" "$PYTHON_SCRIPT" "$PYTHON_EXECUTABLE" "${YAML_FILES[@]}" > "$LOG_FILE" 2>&1 &

# Get the Process ID (PID) of the backgrounded nohup job
PID=$!

echo "--- Process successfully launched! ---"
echo "PID: $PID"
echo "Output is being logged to: $LOG_FILE"
echo "You can now safely close this terminal."
