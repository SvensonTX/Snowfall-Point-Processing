
# import math
import numpy as np
import pdal
import json
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime


# Export instrument data for each line scan time (+/- 5 min) as CSV

# Read and filter Sesame instrument data (received data from Ned Bair)
SesameData = pd.read_csv(r'D:\Research\OperationSnowflake\Data/Sesame2020.csv')
SesameData = SesameData.drop([0,1])  # Drop
SesameData['TIMESTAMP'] = pd.to_datetime(SesameData['TIMESTAMP'])
SesameData = SesameData.set_index("TIMESTAMP")
SesameData = SesameData.drop(SesameData.index[:21000])
SesameData = SesameData.drop(SesameData.index[370000:])

# Read disdrometer data (received data from Ned Bair)
Disdrometer = pd.read_csv(r'D:\Research\OperationSnowflake\Data\20200526-ParsivelData.csv')

# Set directory of .laz files of interest
path= r'D:\Research\OperationSnowflake\Data\LineScansSOCS/'

# Set directory to save disdrometer CSV files
output = r'D:\Research\OperationSnowflake\Data\LineDataCSV\IntensityReturns\Sep_Nov\DisdroCSV/'

# Create list of .laz files
Scans = [os.path.splitext(file)[0] for file in os.listdir(path) if os.path.splitext(file)[1] == '.laz']


# Convert DataFrame column to date time
Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] = pd.to_datetime(Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"])

# Create new DataFrame and set date time index
DisTimeIdx = Disdrometer.set_index("Date Time (mm/dd/yy hh:mm:ss)")

# Convert DataFrame columns to date time
DisTimeIdx["Date"] = pd.to_datetime(DisTimeIdx['Date'])
DisTimeIdx["Time"] = pd.to_datetime(DisTimeIdx["Time"])

# Merge disdrometer data with instrument data collected at Sesame
Merge = DisTimeIdx.merge(SesameData, left_index=True, right_index=True, how='inner')
Merge['DateTime'] = Merge.index

# Set daylight savings time stamp
DaylightSavings = pd.Timestamp(year=2020, month=3, day=8,hour=2)

# Open list for empty data
EmptySets = []

# Batch process for list of .laz files
for i in range(len(Scans)):

    # Extract date and time from file name
    day = Scans[i][6:8]
    month = Scans[i][4:6]
    year = Scans[i][2:4]
    hour = Scans[i][9:11]
    minute = Scans[i][11:13]
    second = Scans[i][14:16]

    # Round up minute if second is larger than 30 and minute is not equal to 59
    if int(second) > 30 and int(minute) is not 59:
        minute = str(int(minute)+1).zfill(2)

    # Create time stamp of scan time
    scanTime = pd.Timestamp(year=int('20'+year),month=int(month),day=int(day),hour=int(hour),minute=int(minute))

    # Create array of times 5 minutes before and after scanTime
    LineTime = pd.Series(Disdrometer.loc[
        ((Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] >= scanTime + pd.Timedelta('-5 minutes')))
        & ((Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] <= scanTime + pd.Timedelta('5 minutes'))),
        'Date Time (mm/dd/yy hh:mm:ss)']).values

    # If LineTime has no values, add to EmptySets list
    if not LineTime.size:

        EmptySets.append(Scans[i]+'_SOCS.laz')
        continue

    # Check if scan happened before daylights savings, adjust disdrometer date time accordingly
    else:
        if LineTime[5] < DaylightSavings:
            print('PST')  # Printout that data processed is standard time
            LineTime = Disdrometer.loc[
                ((Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] >= scanTime + pd.Timedelta('-5 minutes')))
                & ((Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] <= scanTime + pd.Timedelta('5 minutes'))),
                'Date Time (mm/dd/yy hh:mm:ss)']

            LineData = Merge[
                ((Merge["DateTime"] >= scanTime + pd.Timedelta('-8 hours 5 minutes')))
                & ((Merge["DateTime"] <= scanTime + pd.Timedelta('-7 hours 55 minutes')))]

        # Check if scan happened after daylights savings, adjust disdrometer date time accordingly
        else:
            print('PDT')
            LineTime = Disdrometer.loc[
                ((Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] >= scanTime + pd.Timedelta('-5 minutes')))
                & ((Disdrometer["Date Time (mm/dd/yy hh:mm:ss)"] <= scanTime + pd.Timedelta('5 minutes'))),
                'Date Time (mm/dd/yy hh:mm:ss)']
            LineData = Merge[
                ((Merge["DateTime"] >= scanTime + pd.Timedelta('-7 hours 5 minutes')))
                & ((Merge["DateTime"] <= scanTime + pd.Timedelta('-6 hours 55 minutes')))]



    LineData = LineData.drop("Date Time (mm/dd/yy hh:mm:ss)", axis=1)
    LineTime = LineTime.reset_index(drop=True)
    LineData = LineData.reset_index(drop=True)

    # Combine LineTime (time of scan) and LineData (corresponding instrument data)
    Results = pd.concat([LineTime, LineData], axis=1)

    # # Drop columns of no interest
    # Results = Results.drop(
    # ["Weather code SYNOP WaWa", "Weather code METAR/SPECI", "Weather code NWS", "Temperature in sensor (∞C)",
    #  "Heating current (A)", "Sensor voltage (V)", "Kinetic Energy"], axis=1)

    # Set filename for instrument data CSV
    savecsv = year+month+day+hour+minute+'00.csv'

    # Export CSV
    Results.to_csv(output + savecsv)



# Remove files in EmptySets from directory
# for f in EmptySets:
#     f = f.replace('_SOCS','',1)
#     os.remove(r'D:\Research\OperationSnowflake\Data\LineScansSOCS/' + f)



