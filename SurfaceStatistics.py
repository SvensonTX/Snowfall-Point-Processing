
# import math
import numpy as np
import pdal
import os
import json
import datetime
import csv

# Generate surface statistics from a frame scan





# Set directory to .laz files
path =r'D:\Research\OperationSnowflake\Data\Frames\FrameScans/'

# Set directory to save ground extracts
savepath = r'D:\Research\OperationSnowflake\Data\Frames\FrameArea1/'

# Create list of .laz files
Scans = [os.path.splitext(file)[0] for file in os.listdir(path) if os.path.splitext(file)[1] == '.laz']

# Create empty FrameStats list
FrameStats = []

# Batch extract of ground surface AOI to generate statistics
for i in range(len(Scans)):
    infilename = Scans[i] +'.laz'

    # Set crop bounds of interest and filter ground points with elevation and intensity statistics
    json_string = [
                {
                    "type": "readers.las",
                    "filename": path +infilename,
                    "compression": "laszip"
                },
                {
                    "type": "filters.crop",
                    # Bounding coordinates ([xmin, xmax], [ymin, ymax])
                    "bounds": "([849241,849242],[4174330.5,4174331.5])"
                },
                {
                    "type":"filters.smrf"  # Classify ground returns
                },
                {
                    "type":"filters.range",
                    "limits":"Classification[2:2]"
                },
                {
                    "type": "filters.stats",
                    "dimensions": "Z,Intensity"
                },
                {
                    "type": "writers.las",
                    "filename": savepath +infilename+'_Area1.las'
                }
            ]
    pipeline = pdal.Pipeline(json.dumps(json_string))
    pipeline.validate()
    pipeline.execute()
    m = json.loads(pipeline.metadata)

    # Set date time of scan
    Date = datetime.datetime(year=int(Scans[i][0:4]),month=int(Scans[i][4:6]),day=int(Scans[i][6:8]),hour=int(Scans[i][9:11]),minute=int(Scans[i][11:13]))

    # Statistics
    # Z Values
    Avg = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('average')
    Cnt = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('count')
    Max = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('maximum')
    Min = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('minimum')
    Std = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('stddev')
    Var = m.get('metadata').get('filters.stats')[0].get('statistic')[0].get('variance')
    # Intensity Values
    AvgInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('average')
    MaxInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('maximum')
    MinInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('minimum')
    StdInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('stddev')
    VarInt = m.get('metadata').get('filters.stats')[0].get('statistic')[1].get('variance')

    # Append statistics to FrameStats list
    FrameStats.append([Date,Cnt,Min,Max,Avg,Std,Var,MinInt,MaxInt,AvgInt,StdInt,VarInt])

# Create CSV of surface statistics
file = open(savepath+'FrameStatsArea1.csv', 'w+', newline='')
with file:
    write = csv.writer(file)
    write.writerows(FrameStats)