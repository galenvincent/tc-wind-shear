from numpy.core.numeric import normalize_axis_tuple
import pandas as pd
import numpy as np
import xarray as xr
from pandarallel import pandarallel
import time

import credentials
import tc_functions as fun
import plotting_functions as tcplt

storm_data = pd.read_csv('data/filtered_storm_list_keep-leading-5.csv')
storm_data["DATETIME"] = pd.to_datetime(storm_data["DATETIME"])

def int_circulation_storm(id, storm_data, r, normalize, plt_folder, data_folder, upper = False, plot = False):
    
    storm = storm_data[storm_data['ID'].str.match(id)]
    storm = storm.reset_index(drop = True)
    int_circ = []
    for index, datapoint in storm.iterrows():

        year = datapoint["DATETIME"].year
        month = datapoint["DATETIME"].month
        day = datapoint["DATETIME"].day
        hour = datapoint["DATETIME"].hour

        gfs_data = fun.gfs_access(year, month, day, hour, 
                                    credentials.RDA_USER, credentials.RDA_PASSWORD)

        print("Doing #" + str(index) + "/" + str(storm.shape[0]-1))

        # Use upper level winds or shear?
        if upper:
            vws = fun.wind_stamp(datapoint['LAT'], datapoint['LON'], 800, 200, gfs_data,        
                              vortex_rm = False, vortex_rm_rad = 650)
        else: 
            vws = fun.shear_stamp(datapoint['LAT'], datapoint['LON'], 800, gfs_data,
                              vortex_rm = False, vortex_rm_rad = 650)

        ic = fun.integrated_circulation(vws, r, normalize)
        int_circ.append(ic) # Use this later if you want

        if plot:
            tcplt.two_shade_map(vws, ic, 
                                shading = np.arange(-2.,2.,.05), 
                                ticks = np.arange(-2.,2.,0.5), 
                                savefile = plt_folder + id + "_" + str(index) + ".png",
                                legend_title = "Integrated Circulation")
        
        np.save(data_folder + id + "_" + str(index) + ".npy", ic)

#plt_folder = "/glade/work/galenv/int_circ_figs/"
#data_folder = "/glade/work/galenv/int_circ_data/"

radius = 150
normalize_option = "log"

plt_folder = "data/test/"
data_folder = "data/test/"

unique_storms = pd.Series(np.unique(storm_data['ID']))

print("Getting GFS data warmup...")
gfs_data = fun.gfs_access(2016, 12, 12, 0, credentials.RDA_USER, credentials.RDA_PASSWORD)
print("GFS data has been gotten! On to the parallel")

time.sleep(3)

print("Setting up parallel env.")
pandarallel.initialize()
print("Parallel env set up... starting parallel computations.")
unique_storms.parallel_apply(int_circulation_storm, 
                            args = (storm_data, radius, normalize_option, plt_folder, data_folder, False, False))
print("All done!")

#unique_storms.iloc[3:7].parallel_apply(int_circulation_storm, 
#                                        args = (storm_data, radius, normalize_option, plt_folder, data_folder, False, False))