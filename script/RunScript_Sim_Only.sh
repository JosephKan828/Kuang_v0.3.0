#!/bin/sh

# Metadata
PROJECT="Kuang2008_v0.3.0"
EXP_NAME="Rad(0.0,0.0,0.1)"

MODEL_TYPE="full"
RAD_TYPE="all"

TimeStamp=$(date +%Y-%m-%d-%H-%M-%S)

# determine label for each run
RUN_ID=${1:-$TimeStamp}

# Setup directory for saving results
HOME_DIR="${HOME}/${PROJECT}/"
OUTPUT_DIR="${WORK}/${PROJECT}/${MODEL_TYPE}/${EXP_NAME}/${RUN_ID}/"
Fig_DIR="${HOME}/${PROJECT}_Analysis/${MODEL_TYPE}/${EXP_NAME}/${RUN_ID}/"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${Fig_DIR}"

# ====================================================
# Major simulation execution
# ====================================================

julia --quiet \
    --project="$HOME_DIR"\
    --threads 2\
    "$HOME_DIR/script/Simulation.jl" "$MODEL_TYPE" "$RAD_TYPE" "$OUTPUT_DIR" & 
    PID=$!

echo "Started Julia simulation with PID $PID. Monitoring memory and CPU usage..."
echo "Output directory: $OUTPUT_DIR"
echo "Timestamp - Memory (%) - CPU (%)"
echo "---------------------------------"

# Monitor in a loop until the process ends
while kill -0 $PID 2>/dev/null; do
    # Get memory % and CPU % using ps
    INFO=$(ps -p $PID -o %mem,%cpu --no-headers 2>/dev/null)
    if [ -n "$INFO" ]; then
        MEM=$(echo $INFO | awk '{print $1}')
        CPU=$(echo $INFO | awk '{print $2}')
        echo "$(date +%H:%M:%S) - $MEM% - $CPU%"
    fi
    sleep 0.5  # Check every second; adjust for less/more frequency
done

echo "Simulation finished."

# ====================================================
# Calculate eigenvalues and eigenvectors
# ====================================================
echo "Starting eigenvalue and eigenvector analysis..."
julia --quiet \
    --project="$HOME_DIR"\
    --threads 2\
    "$HOME_DIR/script/EigenAnalysis.jl" "$MODEL_TYPE" "$RAD_TYPE" "$OUTPUT_DIR" &
    PID=$!

echo "Started Julia eigen analysis with PID $PID. Monitoring memory and CPU usage..."
echo "Timestamp - Memory (%) - CPU (%)"
echo "---------------------------------"
# Monitor in a loop until the process ends
while kill -0 $PID 2>/dev/null; do
    # Get memory % and CPU % using ps
    INFO=$(ps -p $PID -o %mem,%cpu --no-headers 2>/dev/null)
    if [ -n "$INFO" ]; then
        MEM=$(echo $INFO | awk '{print $1}')
        CPU=$(echo $INFO | awk '{print $2}')
        echo "$(date +%H:%M:%S) - $MEM% - $CPU%"
    fi
    sleep 0.5  # Check every second; adjust for less/more frequency
done

ln -sfn "$OUTPUT_DIR" "${WORK}/${PROJECT}/${MODEL_TYPE}/${EXP_NAME}/latest"


# ====================================================
# Post-processing and visualization
# ====================================================
# activate existed conda environment
. /work/b11209013/external/miniforge3/etc/profile.d/conda.sh
conda activate /work/b11209013/external/miniforge3/envs/dataana

# Plot growth rate and frequency

echo "Starting visualizing growth rate and phase speed..."

python "$HOME_DIR/postprocess/Plot_Eigen.py" "$OUTPUT_DIR" "$Fig_DIR"

echo "Visualization completed. Figures saved to $Fig_DIR"
