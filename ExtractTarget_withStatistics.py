
# import math
import numpy as np
import pdal
import os
import json
import datetime
import csv

# Extract point cloud of interest

# Set directory to .laz files
path =r'D:\Research\OperationSnowflake\Data\Frames\FrameScans/'

# Set directory to save point cloud extract
savepath = r'D:\Research\OperationSnowflake\Data\Frames\FrameScans\SepctralonFrameFilter_2ndReturns/'

# Create list of .laz files
Scans = [os.path.splitext(file)[0] for file in os.listdir(path) if os.path.splitext(file)[1] == '.laz']

# Create empty FrameStats list
FrameStats = []

# Extract point clouds from .laz list
for i in range(len(Scans)):
    infilename = Scans[i] + '.laz'


    # Import .laz file and extract point returns of interest with statistics
    json_string = [
                {
                    "type": "readers.las",
                    "filename": path + infilename,
                    "compression": "laszip"
                },
                {
                    "type": "filters.range",
                    # Set coordinate range and return number of interest
                    "limits": "X[849250.125:849250.175],Y[4174337.055:4174337.5],Z[2727.745:2728.15],ReturnNumber[2:2]"
                },
                {
                    "type": "filters.stats",
                    "dimensions": "Intensity"
                },
                {
                    "compression": "laszip",
                    "minor_version": 2,
                    "scale_x": 0.001,
                    "scale_y": 0.001,
                    "scale_z": 0.001,
                    "offset_x": "auto",
                    "offset_y": "auto",
                    "offset_z": "auto",
                    "type": "writers.las",
                    "filename": savepath +Scans[i]+'_FrameSpectralon2ndReturns.las'  # Set output filename
                }
            ]
    pipeline = pdal.Pipeline(json.dumps(json_string))
    pipeline.validate()
    pipeline.execute()
    m = json.loads(pipeline.metadata)

    Date = datetime.datetime(year=int(Scans[i][0:4]),month=int(Scans[i][4:6]),day=int(Scans[i][6:8]),hour=int(Scans[i][9:11]),minute=int(Scans[i][11:13]))
   # # Z Values
   #  Avg = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('average')
    Cnt = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('count')
   #  Max = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('maximum')
   #  Min = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('minimum')
   #  Std = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('stddev')
   #  Var = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('variance')
   #  # Intensity Values
    AvgInt = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('average')
   #  MaxInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('maximum')
   #  MinInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('minimum')
    StdInt = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('stddev')
   #  VarInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('variance')
    FrameStats.append([Date,Cnt,AvgInt,StdInt])


# Generate CSV of statistics from point cloud extracts
file = open(savepath+'SpectralonFrameStats_2ndReturns.csv', 'w+', newline='')
with file:
    write = csv.writer(file)
    write.writerows(FrameStats)