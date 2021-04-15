import pandas as pd
import numpy as np

import credentials
import hurdat as h
import tc_functions as fun

# Import storms since 2015 from Kim's dataset that includes shortest distance to land, and do some maneuvering to make
# things closer to the HURDAT data straight from the website for east integration with the Hurdat() class (like column names, 
# naming convention for unnamed storms, and some additional columns)
nal_15 = pd.read_table("data/DTL_NAL_2015.txt", delimiter = ' ')
nal_16 = pd.read_table("data/DTL_NAL_2016.txt", delimiter = ' ')
nal_17 = pd.read_table("data/DTL_NAL_2017.txt", delimiter = ' ')
nal_18 = pd.read_table("data/DTL_NAL_2018.txt", delimiter = ' ')
nal_19 = pd.read_table("data/DTL_NAL_2019.txt", delimiter = ' ')
enp_15 = pd.read_table("data/DTL_ENP_2015.txt", delimiter = ' ')
enp_16 = pd.read_table("data/DTL_ENP_2016.txt", delimiter = ' ')
enp_17 = pd.read_table("data/DTL_ENP_2017.txt", delimiter = ' ')
enp_18 = pd.read_table("data/DTL_ENP_2018.txt", delimiter = ' ')
enp_19 = pd.read_table("data/DTL_ENP_2019.txt", delimiter = ' ')
total = pd.concat([nal_15, nal_16, nal_17, nal_18, nal_19, enp_15, enp_16, enp_17, enp_18, enp_19], ignore_index = True)

total.columns = ['YEAR', "PART_ID", "NAME", "DATETIME", "LAT", "LON", "WIND", "PRESSURE", "CATEGORY", "DISTANCE"]

total['DATETIME'] = pd.to_datetime(total['DATETIME'], format = "%Y%m%d%H")
total['NAME'] = total['NAME'].str.upper()
name_convert = {'UNNAMED04': 'FOUR', 
                'UNNAMED11': 'ELEVEN',
                'UNNAMED08': 'EIGHT'}
total['NAME'].loc[total['NAME'].str.contains('UNNAMED')] = [name_convert[x] for x in total['NAME'].loc[total['NAME'].str.contains('UNNAMED')]]
total['ID'] = total['PART_ID'] + total['YEAR'].astype(str)
hurdat_all = h.Hurdat(data = total)
# Filter storms by minimum intensity of 50
hurdat = hurdat_all.genesis_to_lysis_filter(minimum_wind = 50)
# Add RI / RW labels to all storm observations
hurdat.identify_events(threshold = 25)
# Add 250 km from land restriction as a column
hurdat.distance_to_land_label(min_distance = 250)

row = hurdat.storms.iloc[626]
row

year = row["DATETIME"].year
month = row["DATETIME"].month
day = row["DATETIME"].day
hour = row["DATETIME"].hour

print(row["ID"])

gfs_data = fun.gfs_access(year, month, day, hour, 
                            credentials.RDA_USER, credentials.RDA_PASSWORD)

# Find center in GFS data using best track as a seed:
center_lat, center_lon = fun.vorticity_centroid(row['LAT'], row['LON'], 
                                                pressure = 850, search_radius = 200, 
                                                calc_radius = 150, dataset = gfs_data)

p_lat, p_lon = fun.pressure_min(row['LAT'], row['LON'], search_radius = 200, dataset = gfs_data)

#fun.storm_radius(center_lat, center_lon, gfs_data, max_radius = 1200,
#                 pressure = 850, stride = 10, h = 25, plot = True)
fun.storm_radius(row['LAT'], row['LON'], gfs_data, max_radius = 1200,
                 pressure = 850, stride = 10, h = 25, plot = True)