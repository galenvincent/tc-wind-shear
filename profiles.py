import pandas as pd
import numpy as np
import xarray as xr
from pandarallel import pandarallel
import time

import credentials
import hurdat as h
import tc_functions as fun
import plotting_functions as tcplt

storm_data = pd.read_csv('data/filtered_storm_list.csv')
storm_data["DATETIME"] = pd.to_datetime(storm_data["DATETIME"])

def profile_storm(id, storm_data, shear_plt_folder, profiles_folder, plot = False):
    
    storm = storm_data[storm_data['ID'].str.match(id)]
    storm = storm.reset_index(drop = True)
    vws = []
    centers = []
    profiles = []
    badIndex = -1
    for index, datapoint in storm.iterrows():
        year = datapoint["DATETIME"].year
        month = datapoint["DATETIME"].month
        day = datapoint["DATETIME"].day
        hour = datapoint["DATETIME"].hour
        
        if year == 2018 and month == 5 and day == 26 and hour == 18:
            centers.append([datapoint['LAT'], datapoint['LON']])
            vws.append(None)
            badIndex = index
            continue

        gfs_data = fun.gfs_access(year, month, day, hour, 
                                    credentials.RDA_USER, credentials.RDA_PASSWORD)
        
        ### TODO: Figure out this whole storm centering thing. 
        # Find center in GFS data using best track as a seed:
        #centers.append(fun.vorticity_centroid(datapoint['LAT'], datapoint['LON'], 
        #                                      pressure = 850, search_radius = 200, calc_radius = 150,
        #                                      dataset = gfs_data))
        centers.append([datapoint['LAT'], datapoint['LON']])

        print("Doing #" + str(index+1) + "/" + str(storm.shape[0]))

        # Get wind shear with vortex removed 
        vws.append(fun.shear_stamp(centers[index][0], centers[index][1], 800, gfs_data,
                vortex_rm = True, vortex_rm_rad = 650))
        
    if badIndex != -1:
        vws[badIndex] = vws[badIndex+1]
        badIndex = -1

    # Calculate storm direction
    for ii in range(1, storm.shape[0]):
        direction = (round(centers[ii][1] - centers[ii - 1][1], 5), round(centers[ii][0] - centers[ii - 1][0], 5))
        
        if direction == (0, 0):
            vws[ii].attrs['storm_direction'] = (0.001, 0)
        else:
            vws[ii].attrs['storm_direction'] = direction

    vws[0].attrs['storm_direction'] = vws[1].attrs['storm_direction']

    # Sectorize according to both storm direction and shear direction, and then calculate radial profiles based on this
    profiles = []
    for ii in range(len(vws)):
        vws[ii] = fun.sectorize(vws[ii], vws[ii].attrs['avg_shear'][0], vws[ii].attrs['avg_shear'][1], "shear")
        vws[ii] = fun.sectorize(vws[ii], vws[ii].attrs['storm_direction'][0], vws[ii].attrs['storm_direction'][1], "velocity")
        profiles.append(fun.radial_profile(vws[ii], stride = 10, h = 35, sector_labels = ["shear", "velocity"]))
        if plot:
            tcplt.shear_map(vws[ii], shear_plt_folder + id + "_" + str(ii) + ".png")

    # Collect all storm info into one xarray DataSet
    stamp_radius = profiles[0].attrs['stamp_radius']
    radius = np.array(profiles[0].radius)
    component = np.array(profiles[0].component)
    sector = np.array(profiles[0].sector)
    time = np.array(storm["DATETIME"])
    wind = np.array(storm["WIND"])
    ri = np.array(storm["RI"])
    rw = np.array(storm["RW"])
    leading = np.array(storm["LEADING"]) # Might need to remove this if there is no leading column
    near_land = np.array(storm["NEAR_LAND"])
    dist_to_land = np.array(storm["DISTANCE"])
    center_lat = [profile.attrs['center_lat'] for profile in profiles]
    center_lon = [profile.attrs['center_lon'] for profile in profiles]
    velocity_direction_u = [profile.attrs['sector_velocity_direction'][0] for profile in profiles]
    velocity_direction_v = [profile.attrs['sector_velocity_direction'][1] for profile in profiles]
    shear_direction_u = [profile.attrs['sector_shear_direction'][0] for profile in profiles]
    shear_direction_v = [profile.attrs['sector_shear_direction'][1] for profile in profiles]
    shear_magnitude = np.sqrt(np.power(shear_direction_u, 2) + np.power(shear_direction_v, 2))

    shear_oriented = np.array([profile.values[:,:,0,:] for profile in profiles])
    velocity_oriented = np.array([profile.values[:,:,1,:] for profile in profiles])

    storm_profiles = xr.Dataset(
                    {
                        'shear_oriented': (['time', 'component', 'radius', 'sector'], shear_oriented),
                        'velocity_oriented': (['time', 'component', 'radius', 'sector'], velocity_oriented)
                    },
                    coords = {
                        'time': time,
                        'component': component,
                        'radius': radius,
                        'sector': sector,
                        'ri': ('time', ri),
                        'rw': ('time', rw),
                        'near_land': ('time', near_land),
                        'wind': ('time', wind),
                        'dist_to_land': ('time', dist_to_land),
                        'leading': ('time', leading),
                        'center_lat': ('time', center_lat),
                        'center_lon': ('time', center_lon),
                        'shear_u': ('time', shear_direction_u),
                        'shear_v': ('time', shear_direction_v),
                        'shear_mag': ('time', shear_magnitude),
                        'velocity_u': ('time', velocity_direction_u),
                        'velocity_v': ('time', velocity_direction_v)  
                    },
                    attrs = {
                        'long_name': "Set of all radial profiles for an entire storm.",
                        'storm_id': id,
                        'stamp_radius': stamp_radius
                    }
    )

    storm_profiles.to_netcdf(profiles_folder+"prof_"+id+".nc")

shear_plt_folder = "/glade/work/galenv/shear_figs/"
profiles_folder = "/glade/work/galenv/shear_profiles/"

#shear_plt_folder = "figures/"
#profiles_folder = "data/profiles/"

unique_storms = pd.Series(np.unique(storm_data['ID']))

print("Getting GFS data warmup...")
gfs_data = fun.gfs_access(2016, 12, 12, 0, credentials.RDA_USER, credentials.RDA_PASSWORD)
print("GFS data has been gotten! On to the parallel")

time.sleep(3)

print("Setting up parallel env.")
pandarallel.initialize()
print("Parallel env set up... starting parallel computations.")
unique_storms.parallel_apply(profile_storm, args = (storm_data, shear_plt_folder, profiles_folder, False))

print("All done!")

#unique_storms.iloc[2:5].apply(profile_storm, args = (storm_data, shear_plt_folder, profiles_folder, False))