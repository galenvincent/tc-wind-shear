import credentials
import bounding_box as bb

from siphon.catalog import TDSCatalog
from siphon.http_util import session_manager
from datetime import datetime
import numpy as np
import xarray as xr 
import math

# Set RDA login from credentials.py
session_manager.set_session_options(auth = (credentials.RDA_USER, credentials.RDA_PASSWORD))

def get_shear_stamp(year, month, day, hour, lat, lon, radius, username, password):
    """
    Retrieves a circular stamp of wind shear from GFS analysis (not forecast) 
    data associated with input parameter specifications. 

    Parameters:
        - year (int): after 2015, which is historical GFS data avalability for 
        the best resolution.
        - month (int)
        - day (int)
        - hour (int): 0, 6, 12, or 18. These are the hours that GFS is run.
        - lat (float): latitude for the center of the stamp (-90 - 90)
        - lon (float):  longitude for the center of the stamp (-180 - 180).
        - radius (float): radius of the circular stamp (km).
        - username (str): RDA username for accessing GFS data. 
        - password (str): RDA password for accessing GFS data.

    Returns: 

    """

    # Set RDA account login
    session_manager.set_session_options(auth = (username, password))

    # Convert time info into datetime object for easy maneuvering
    assert year >= 2015, "Year must be >= 2015"
    assert month >= 1 and month <= 12, "Month must be in [1, 12]"
    assert hour in [0, 6, 12, 18], "Hour must be one of 0, 6, 12, 18"
    ymdh = datetime(year, month, day, hour)

    # Set catalog URL and dataset name from time info
    cat_url = ymdh.strftime("https://rda.ucar.edu/thredds/catalog/files/g/ds084.1/%Y/%Y%m%d/catalog.xml")
    dataset_name = ymdh.strftime("gfs.0p25.%Y%m%d%H.f000.grib2")

    # Get remote access to data
    catalog = TDSCatalog(cat_url)
    ds = catalog.datasets[dataset_name]
    dataset = ds.remote_access()

    # Pull zonal and meridional wind data at isobaric surfaces
    zonal_isobaric = dataset.variables["u-component_of_wind_isobaric"]
    merid_isobaric = dataset.variables["v-component_of_wind_isobaric"]

    iso_number = zonal_isobaric.dimensions[1]
    assert iso_number == merid_isobaric.dimensions[1], "Zonal and meridional wind vectors have different pressure sets in data."

    pressure_levels = dataset.variables[iso_number][:]
    data_lat = dataset.variables["lat"][:]
    data_lon = dataset.variables["lon"][:]

    # Get 200 & 850 hPa indices
    pressure_200_ind = np.where(pressure_levels == 20000)[0][0]
    pressure_850_ind = np.where(pressure_levels == 85000)[0][0]

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth, hence the fancy code in bounding_box)
    center = bb.GeoLocation.from_degrees(lat, lon)
    min_lat, max_lat, min_lon, max_lon = center.bounding_locations(radius)

    # Get indices for the bounding lat/lon
    lat_rows = np.where((data_lat >= min_lat) & (data_lat <= max_lat))[0]
    lon_rows = np.where((data_lon >= min_lon) & (data_lon <= max_lon))[0]
    min_lat_ind, max_lat_ind = [min(lat_rows), max(lat_rows) + 1]
    min_lon_ind, max_lon_ind = [min(lon_rows), max(lon_rows) + 1]

    # Get the wind data
    u_200 = zonal_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
    u_850 = zonal_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
    v_200 = merid_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
    v_850 = merid_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]

    # Get lat and lon corresponding to just our data of interest and turn it into 
    # a grid to match the wind data
    lat_subset = data_lat[min_lat_ind:max_lat_ind]
    lon_subset = data_lon[min_lon_ind:max_lon_ind]
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat_subset, lon_subset)]

    # Subset out a circular stamp of the data, putting nan everywhere else
    
    # Distance from each point to the center
    dist_mat = great_circ_dist(lat, lon, lat_grid, lon_grid)
    dist_tf = dist_mat > radius
    # Set all other points to nan
    u_200_stamp = np.where(dist_tf, np.nan, u_200)
    u_850_stamp = np.where(dist_tf, np.nan, u_850)
    v_200_stamp = np.where(dist_tf, np.nan, v_200)
    v_850_stamp = np.where(dist_tf, np.nan, v_850)

    da = xr.DataArray(np.array([[u_200, v_200], [u_850, v_850]]), 
                      coords = [[200, 850], ["u", "v"], lat_subset, lon_subset],
                      dims = ["Pressure", "Direction", "lat", "lon"],
                      name = "wind",
                      attrs = {
                          "long_name": "Wind values for u (zonal) and v (meridional) directions.",
                          "Units": "meters/second"
                      })

    #da = xr.DataArray(np.array([[u_200_stamp, v_200_stamp], [u_850_stamp, v_850_stamp]]), 
    #                  coords = [[200, 850], ["u", "v"], lat_subset, lon_subset],
    #                  dims = ["Pressure", "Direction", "lat", "lon"],
    #                  name = "wind",
    #                  attrs = {
    #                      "long_name": "Wind values for u (zonal) and v (meridional) directions.",
    #                      "Units": "meters/second"
    #                  })                  

     # Create masked arrays (for computations)
    u_200_stamp = np.ma.array(u_200_stamp, mask = np.isnan(u_200_stamp))
    u_850_stamp = np.ma.array(u_200_stamp, mask = np.isnan(u_850_stamp))
    v_200_stamp = np.ma.array(u_200_stamp, mask = np.isnan(v_200_stamp))
    v_850_stamp = np.ma.array(u_200_stamp, mask = np.isnan(v_850_stamp))
    
    return da

def great_circ_dist(lat1, lon1, lat2, lon2):
    """
    Vectorized function calculates the great circle distance between two points 
    of (lat, lon).

    Parameters:
        - lat1 (float): Latitude of first point (degrees).
        - lon1 (float): Longitude of first point (degrees).
        - lat2 (float): Latitude of second point (degrees).
        - lon2 (float): Longitude of second point (degrees).
        - radius (float): The radius of the sphere for which you want to compute.
        Default is the columetric mean radius of the earth in kilometers.

    Returns:
        - (float) a scalar (or matrix, depending on inputs) of the distance 
        between point 1 and point 2.

    """
    earth_radius = 6371.0

    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    dist = earth_radius * np.arccos(
                            np.sin(lat1_rad) * np.sin(lat2_rad) + 
                            np.cos(lat1_rad) * np.cos(lat2_rad) * np.cos(lon1_rad - lon2_rad)
                          )
    return dist

def circular_cut(center_lat, center_lon, lat_grid, lon_grid, stamp_radius):
    """
    Helper function to get a circular stamp of wind data from a grid of data.
    """
    
    dist_mat = great_circ_dist(center_lat, center_lon, lat_grid, lon_grid)

    return dist_mat <= stamp_radius



get_shear_stamp(2019, 1, 10, 18, 
                39.656456, -105.012831, 800, 
                credentials.RDA_USER, credentials.RDA_PASSWORD)





