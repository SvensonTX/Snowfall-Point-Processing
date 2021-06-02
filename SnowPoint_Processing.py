import pdal
import json
import os
import numpy as np


# Batch transform of point clouds

# Set directory path to .laz files
path = r'C:\Users\sasorhus\PycharmProjects\Snowflake\LineData\9to12_2019/'

# Set directory path to save transformed files
save = r'D:\Research\OperationSnowflake\Data\LineScansSOCS/'

# Create list of .laz files within a directory
Scans = [os.path.splitext(file)[0] for file in os.listdir(path) if os.path.splitext(file)[1] == '.laz']

# Transform list of .laz files to SOCS at Sesame
for i in range(len(Scans)):
    infilename = path + Scans[i] +'.laz'  # Input point cloud
    outfilename = save + Scans[i] +'_SOCS.laz'  # Set output point cloud name

    json_array = [
        {
            "type": "readers.las",
            "compression": "laszip",
            "filename": infilename  # Read in point cloud
        },
        {
            "type": "filters.reprojection",
            "in_srs": "EPSG:32610",  # Input coordinate system - UTM zone 10N
            "out_srs": "EPSG:4978"   # Output coordinate system - WGS84
        },
        {

            "type": "filters.transformation",  # POP transformation
            "matrix": "0.874268495900 -0.485442681556 0.000000000000 0.000000000000 0.296522721271 0.534029007700 0.791762145284 20675.785437908489 -0.384355138961 -0.692212699869 0.610829522285 -6370166.497704507783 0.000000000000 0.000000000000 0.000000000000 1.000000000000"
        },
        {
            "type": "filters.transformation",  # SOP transformation
            "matrix": "0.243466990948 0.834583002822 0.494160940895 -1371.356756454145 -0.958570643796 0.284725256939 -0.008593455998 -7.178285891526 -0.147872052968 -0.471595948900 0.869327968815 -2362.346996182266 0.000000000000 0.000000000000 0.000000000000 1.000000000000"
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
            "filename": outfilename  # Write output point cloud
        }

    ]

    pipeline = pdal.Pipeline(json.dumps(json_array))
    pipeline.validate()
    pipeline.execute()





# Batch export of .laz point data to CSV

# Set directory path for desired point cloud
path =r'D:\Research\OperationSnowflake\Data\LineScansSOCS/'

# Set directory path to save CSV
savepath = r'D:\Research\OperationSnowflake\Data\LineDataCSV\IntensityReturns\Sep_Nov/'

# Create list of .laz files within a directory
Scans = [os.path.splitext(file)[0] for file in os.listdir(path) if os.path.splitext(file)[1] == '.laz']


# Export point data from list of .laz files to CSV
for i in range(len(Scans)):

    infilename = Scans[i] +'.laz'  # Input point cloud

    json_string = [
                {
                    "type": "readers.las",
                    "filename": path +infilename,
                    "compression": "laszip"
                }
            ]
    pipeline = pdal.Pipeline(json.dumps(json_string))
    pipeline.validate()
    pipeline.execute()
    arrays = pipeline.arrays
    view = arrays[0]

    # Create array of point cloud data
    x = view['X']
    y = view['Y']
    z = view['Z']
    Intensity = view['Intensity']
    time = view['GpsTime']
    ReturnNum = view['ReturnNumber']
    Returns = view['NumberOfReturns']

    # Point data array
    Points = np.array([x,y,z,Intensity,time,ReturnNum,Returns])
    Points = Points.T

    # Save CSV file
    np.savetxt(savepath + Scans[i] +'.csv',Points,delimiter=',')

