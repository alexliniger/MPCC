# Lap Time Optimization (LTO) Pre-processing

This directory contains a Python-based lap time optimization tool used to generate optimal racing lines for the MPCC controller.

## Setup

1.  Create a virtual environment:
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Structure

- `casadi_lto.py`: Main optimization script. The model formulation (Magic Formula and dynamics) is aligned with the `ADCodeGen` implementation for consistency.
- `track_utils.py`: Utilities for loading and handling tracks using periodic splines and CasADi interpolants.
- `track_pre_processing.py`: Script to process raw track data into the standard JSON format used by the controller and LTO.

## Usage

1.  (Optional) Process track data:
    ```bash
    python3 track_pre_processing.py
    ```
2.  Run the LTO script:
    ```bash
    python3 casadi_lto.py
    ```

The script will generate an optimal racing line and velocity profile based on the `ADCodeGen` physics model.
