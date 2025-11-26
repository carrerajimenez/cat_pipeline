#!/usr/bin/zsh

# === CONFIGURATION ===
# --- Adjust these variables ---

# Path to your python script
PYTHON_SCRIPT="/media/rcarrera/scratch1/mahdieh/revision_catpaper/code_to_publish/cat_pipeline.py"

# --- Specify your YAML file ---
YAML_FILE="/media/rcarrera/scratch1/mahdieh/revision_catpaper/code_to_publish/config_5K.yaml"

# Python executable (e.g., python3, or /path/to/your/venv/bin/python)
PYTHON_EXECUTABLE="python"

# Log file for nohup output
LOG_FILE="run_on_references_1yaml.log"

# === SCRIPT START ===

echo "--- Single-run script launcher started ---"

# --- 1. Check if the YAML file exists before starting ---
echo "Checking for input YAML file: $YAML_FILE"
if [ ! -f "$YAML_FILE" ]; then
    echo "ERROR: Required YAML file not found: $YAML_FILE"
    echo "Aborting script."
    exit 1
fi
echo "YAML file found. Proceeding..."


# --- 2. Run Python script with nohup ---
echo "Starting Python script in the background..."

nohup zsh -c "
    echo \"--- Starting Run: $YAML_FILE at \$(date) ---\"

    # 1. Try to run the Python command
    if \"$PYTHON_EXECUTABLE\" \"$PYTHON_SCRIPT\" -c \"$YAML_FILE\"; then
        
        echo \"--- Run Succeeded: $YAML_FILE at \$(date) ---\"
        
    else
        echo \"--- ERROR: Run FAILED (\$(basename $YAML_FILE)) at \$(date) ---\"
    fi

    echo \"--- Script finished at \$(date) ---\"
" > "$LOG_FILE" 2>&1 &

# Get the Process ID (PID) of the backgrounded nohup job
PID=$!

echo "--- Process successfully launched! ---"
echo "PID: $PID"
echo "Output is being logged to: $LOG_FILE"
echo "You can now safely close this terminal."