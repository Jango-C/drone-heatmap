# FPV Heatmap Generator

This application generates a signal strength (RSSI) heatmap for FPV drones based on terrain data. It helps visualize potential signal quality for a given flight path and receiver location.

## How to Run
1.  **Clone the repository.**
2.  **Create a virtual environment (recommended):**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```
3.  **Install the required libraries:**
    ```bash
    pip install -r requirements.txt
    ```
4.  **Run the application:**
    Make sure the `elevation_map.tif` file is in the same directory.
    ```bash
    python fpv_map_app.py
    ```