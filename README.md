# Snowfall-Point-Processing
Data and scripts for processing lidar point clouds with instrument data from Sesame utilized during my thesis research at the University of Houston.

# Data

Sesame2020.csv
- Meteorological instrument data collected at Sesame

20200526-ParsivelData.csv
- Parsivel disdrometer data collected from Sesame

InstrumentData.csv
- Disdrometer and Meteorological instrument data records from Sesame

Plookup.mat
- Date of disdrometer data

Pdata.mat
- Date and data. Column 7 contains a cell with a 32 x 32 double of all size and speed estimates

# Scripts

# Python

SnowPoint_Processing.py
- Batch transform of point clouds
- Batch export of point data to CSV

SurfaceStatistics.py
- Extract surface statistics (elevation and intensity) from an area within a point cloud

InstrumentData.py
- Export CSV files of instrument data (adjusted to UTC) around scan time of point cloud

ExtractTarget_withStatistics.py
- Extract target coordinates from a point cloud along with intensity statistics



# Matlab Scripts

ParsivelSpectrums.m
- Matlab function to generate speed and size bar graph counts and median estimates from disdrometer data. 
- Uses Plookup.mat and Pdata.mat 

LoadSesameData.m
- matlab generated function to read CSV of instrument data generated from InstrumentData.py

StreakVelocities.m
- matlab function to automatically estimate streak velocities from a point cloud.
- point cloud is filtered within 5 meters of scanner and has a scalar field of linearity calculated and saved from CloudCompare

SpectralonStats.m
- Read CSV of statistics generated from 'ExtractTarget_withStatistics.py' to a table

LineScanProcessing.m
- Use directory of point cloud data (CSV) and instrument data (CSV) generated from 'SnowPoint_Processing.py' to analyze line scans
- Uses 'ParsivelSpectrums.m'
- Outputs various data to Results table, saves point data and snow mask for every line scan

SnowfallStatistics.m
- Compute snowfall statistics based on snow mask and point data generated from LineScanProcessing.m
