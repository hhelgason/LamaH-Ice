Code used to process the Copernicus Arctic Regional Reanalysis (CARRA) for version 1.5 of the LamaH-Ice dataset.

The processing relies on the xesmf Python package, which uses ESMF under the hood. Due to its dependency on compiled libraries, the environment can be sensitive to package versions and platform compatibility. Environment files are provided to ensure reproducibility. Setting up the environment on Linux is strongly recommended, as the code may not work on Windows.

Environment setup
To recreate the environment:
conda env create -f environment.yml
To recreate the environment with exact package versions and builds:
conda create --name myenv --file env.txt

Run order
To process the CARRA data (available at http://ftp.betravedur.is/LV/icebox/carra/), run the scripts in the following order:
process_carra_parallel.py – regrids and processes raw CARRA files in parallel
combined_watershed_series.py – combines sub-watershed time series into catchment-level series
create_daily_series.py – converts sub-daily data to daily resolution
