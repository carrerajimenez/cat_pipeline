#!/usr/bin/zsh

# === CONFIGURATION ===
# --- Adjust these variables ---

# Path to your python script
PYTHON_SCRIPT="/media/rcarrera/scratch1/mahdieh/revision_catpaper/code_to_publish/cat_pipeline.py"

# --- Specify your three YAML files ---
YAML_FILE_1="/media/rcarrera/scratch1/mahdieh/revision_catpaper/referenced_spectra/config_85K.yaml"
YAML_FILE_2="/media/rcarrera/scratch1/mahdieh/revision_catpaper/referenced_spectra/config_10K.yaml"
YAML_FILE_3="/media/rcarrera/scratch1/mahdieh/revision_catpaper/referenced_spectra/config_5K.yaml"

# Python executable (e.g., python3, or /path/to/your/venv/bin/python)
PYTHON_EXECUTABLE="python"

# Log file for nohup output
LOG_FILE="run_on_references.log"

# === SCRIPT START ===

echo "--- Sequential script launcher started ---"

# --- 1. Check if all YAML files exist before starting ---
echo "Checking for input YAML files..."
for yaml_file in "$YAML_FILE_1" "$YAML_FILE_2" "$YAML_FILE_3"; do
    if [ ! -f "$yaml_file" ]; then
        echo "ERROR: Required YAML file not found: $yaml_file"
        echo "Aborting script."
        exit 1
    fi
done
echo "All 3 YAML files found. Proceeding..."


# --- 2. Run Python scripts sequentially with nohup ---
echo "Starting Python scripts in the background..."

nohup zsh -c "
    echo \"--- Starting Run 1: $YAML_FILE_1 at \$(date) ---\"

    # 1. Try to run the FIRST Python command
    if \"$PYTHON_EXECUTABLE\" \"$PYTHON_SCRIPT\" -c \"$YAML_FILE_1\"; then
        
        echo \"--- Finished Run 1. Starting Run 2: $YAML_FILE_2 at \$(date) ---\"
        
        # 2. If first succeeded, run the SECOND Python command
        if \"$PYTHON_EXECUTABLE\" \"$PYTHON_SCRIPT\" -c \"$YAML_FILE_2\"; then

            echo \"--- Finished Run 2. Starting Run 3: $YAML_FILE_3 at \$(date) ---\"
            
            # 3. If second succeeded, run the THIRD Python command
            if \"$PYTHON_EXECUTABLE\" \"$PYTHON_SCRIPT\" -c \"$YAML_FILE_3\"; then
                echo \"--- Run 3 Succeeded. All runs complete. ---\"
            else
                echo \"--- ERROR: Run 3 FAILED (\$(basename $YAML_FILE_3)). ---\"
            fi
        else
            echo \"--- ERROR: Run 2 FAILED (\$(basename $YAML_FILE_2)). Run 3 will not start. ---\"
        fi
    else
        echo \"--- ERROR: Run 1 FAILED (\$(basename $YAML_FILE_1)). Runs 2 and 3 will not start. ---\"
    fi

    echo \"--- All runs finished at \$(date) ---\"
" > "$LOG_FILE" 2>&1 &

# Get the Process ID (PID) of the backgrounded nohup job
PID=$!

echo "--- Process successfully launched! ---"
echo "PID: $PID"
echo "Output is being logged to: $LOG_FILE"
echo "You can now safely close this terminal."