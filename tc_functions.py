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

def get_shear_stamp(year, month, day, hour, 
                    lat, lon, radius, 
                    username, password):
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
        (DataArray) A xarray DataArray object containing the zonal (u)/ 
        meridional (v) directions and magnitude for a circular stamp of vertical
        wind shear calculated between the 200 and 850 hPa isobars. 

    Additional Details:
        Use the quadrantize() function to split up the wind shear data into 
        quadrants from a specified direction.
    """

    # Convert time info into datetime object for easy maneuvering
    assert year >= 2015, "Year must be >= 2015"
    assert month >= 1 and month <= 12, "Month must be in [1, 12]"
    assert hour in [0, 6, 12, 18], "Hour must be one of 0, 6, 12, 18"
    ymdh = datetime(year, month, day, hour)

    # Set catalog URL and dataset name from time info
    cat_url = ymdh.strftime("https://rda.ucar.edu/thredds/catalog/files/g/ds084.1/%Y/%Y%m%d/catalog.xml")
    dataset_name = ymdh.strftime("gfs.0p25.%Y%m%d%H.f000.grib2")

    # Get remote access to data
    dataset = get_gfs_data(cat_url, dataset_name, username, password)

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
    min_lat_ind, max_lat_ind = [min(lat_rows) - 1, max(lat_rows) + 1]
    min_lon_ind, max_lon_ind = [min(lon_rows) - 1, max(lon_rows) + 1]

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

    # Calculate wind shear components and magnitude
    u_wind_shear = u_200_stamp - u_850_stamp
    v_wind_shear = v_200_stamp - v_850_stamp
    magnitude_wind_shear = np.sqrt(np.power(u_wind_shear, 2) + np.power(v_wind_shear, 2))

    shear = xr.DataArray(
        data = np.array([magnitude_wind_shear, u_wind_shear, v_wind_shear]),
        coords = {
            "component": ["magnitude", "u", "v"],
            "lat": lat_subset,
            "lon": lon_subset,
        },
        dims = ["component", "lat", "lon"],
        name = "200-850hPa_wind_shear",
        attrs = {
            "long_name": "200_850 hPa vertical wind shear",
            "units": "meters/second",
            "center_lat": lat,
            "center_lon": lon,
            "stamp_radius": radius,
            "time": ymdh

        }
    )

    return shear

def get_wind_stamp(year, month, day, hour, 
                   lat, lon, radius, pressure, 
                   username, password):
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
        - pressure (float): pressure level for wind data (hPa). Closest 
        available pressure level will be returned. 
        - username (str): RDA username for accessing GFS data. 
        - password (str): RDA password for accessing GFS data.

    Returns:
        (DataArray) A xarray DataArray object containing the zonal (u)/ 
        meridional (v) directions and magnitude for a circular stamp of GFS
        isobaric wind data at the specified time, location, and radius. 

    Additional Details:
        Use the quadrantize() function to split up the wind shear data into 
        quadrants from a specified direction.
    """

    # Convert time info into datetime object for easy maneuvering
    assert year >= 2015, "Year must be >= 2015"
    assert month >= 1 and month <= 12, "Month must be in [1, 12]"
    assert hour in [0, 6, 12, 18], "Hour must be one of 0, 6, 12, 18"
    ymdh = datetime(year, month, day, hour)

    # Set catalog URL and dataset name from time info
    cat_url = ymdh.strftime("https://rda.ucar.edu/thredds/catalog/files/g/ds084.1/%Y/%Y%m%d/catalog.xml")
    dataset_name = ymdh.strftime("gfs.0p25.%Y%m%d%H.f000.grib2")

    # Get remote access to data
    dataset = get_gfs_data(cat_url, dataset_name, username, password)

    # Pull zonal and meridional wind data at isobaric surfaces
    zonal_isobaric = dataset.variables["u-component_of_wind_isobaric"]
    merid_isobaric = dataset.variables["v-component_of_wind_isobaric"]

    iso_number = zonal_isobaric.dimensions[1]
    assert iso_number == merid_isobaric.dimensions[1], "Zonal and meridional wind vectors have different pressure sets in data."

    pressure_levels = dataset.variables[iso_number][:]
    data_lat = dataset.variables["lat"][:]
    data_lon = dataset.variables["lon"][:]

    # Get closest pressure index
    pressure = pressure * 100 # convert from hPa to Pa
    pressure_ind = np.argmin(np.abs(pressure_levels - pressure))
    pressure_used = pressure_levels[pressure_ind]/100 # get closest pressure in hPa

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth, hence the fancy code in bounding_box)
    center = bb.GeoLocation.from_degrees(lat, lon)
    min_lat, max_lat, min_lon, max_lon = center.bounding_locations(radius)

    # Get indices for the bounding lat/lon
    lat_rows = np.where((data_lat >= min_lat) & (data_lat <= max_lat))[0]
    lon_rows = np.where((data_lon >= min_lon) & (data_lon <= max_lon))[0]
    min_lat_ind, max_lat_ind = [min(lat_rows) - 1, max(lat_rows) + 1]
    min_lon_ind, max_lon_ind = [min(lon_rows) - 1, max(lon_rows) + 1]

    # Get the wind data
    u = zonal_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
    v = merid_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]

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
    u_stamp = np.where(dist_tf, np.nan, u)
    v_stamp = np.where(dist_tf, np.nan, v)

    # Calculate wind magnitude
    magnitude = np.sqrt(np.power(u_stamp, 2) + np.power(v_stamp, 2))

    wind = xr.DataArray(
        data = np.array([magnitude, u_stamp, v_stamp]),
        coords = {
            "component": ["magnitude", "u", "v"],
            "lat": lat_subset,
            "lon": lon_subset,
        },
        dims = ["component", "lat", "lon"],
        name = "wind_stamp",
        attrs = {
            "long_name": "Circular stamp of wind data from GFS.",
            "units": "meters/second",
            "center_lat": lat,
            "center_lon": lon,
            "pressure_level": pressure_used,
            "stamp_radius": radius,
            "time": ymdh

        }
    )

    return wind

def get_gfs_data(catalog_url, dataset_name, username, password):
    """
    Helper function to return a remote access object with the desired dataset.

    Parameters:
        - catalog_url (str): URL of the RDA catalog the dataset is in.
        - dataset_name (str): Name of the dataset within the catalog. 
        - username (str): RDA account username. 
        - password (str): RDA accound password. 

    Returns:
        An object returned from remote_access function from the siphon module. 
    """
    # Set RDA username and password
    session_manager.set_session_options(auth = (username, password))

    # Get remote access to data
    catalog = TDSCatalog(catalog_url)
    ds = catalog.datasets[dataset_name]
    dataset = ds.remote_access()

    return dataset

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
        (float): a scalar (or matrix, depending on inputs) of the distance 
        between point 1 and point 2, in kilometers.

    """
    # Set earth's radius in kilometers
    earth_radius = 6371.0

    # Convert degrees to radians
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    # Calculate geodesic distance between to lat/lon points on a sphere.
    # Equation given here: http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates
    dist = earth_radius * np.arccos(
            np.sin(lat1_rad) * np.sin(lat2_rad) + 
            np.cos(lat1_rad) * np.cos(lat2_rad) * np.cos(lon1_rad - lon2_rad)
        )

    return dist

def quadrantize(x, dir_u, dir_v):
    """
    Determines quadrants based on an initial direction vector.

    Takes in an xarray DataArray object (x) with lat/lon coordinates and
    an attribute indicating the center latitude and longitude. Then computes
    quadrants based on some direction vector given by (dir_u, dir_v), where
    the quadrants 1, 2, 3, 4 are labeled in a counter-clockwise direction
    starting at (dir_u, dir_v). These quadrant labels are added to the DataArray
    and returned.

    Parameters:
        - x (DataArray): Wind shear or wind data output from the get_shear_stamp()
        or get_wind_stamp() functions. 
        - dir_u (float): zonal (east) component of the direction vector.
        - dir_v (float): meridional (north) component of the direction vector.
    
    Returns:
        (DataArray): An updated DataArray object with a new coordinate to 
        indicate which quadrant each data point is in.
    """
    
    # Get latitude and longitude arrays from the DataArray object
    lat = x.lat.values
    lon = x.lon.values

    # Create lat/lon grids out of the arrays, so we have every individual 
    # (lon, lat) coordinate combination
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat, lon)]    

    # Get the center of the stamp from the DataArray
    center_lat = x.attrs["center_lat"]
    center_lon = x.attrs["center_lon"]

    # Convert to 0-360 scaled longitude to match data from DataArray
    if center_lon < 0:
        center_lon = center_lon + 360

    # Create centered lat/lon data and put into one array
    uv_grid = np.array([lon_grid - center_lon, lat_grid - center_lat])

    # Compute angle between each point in grid and the direction given in 
    # function input.
    # See https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors/16544330#16544330
    # for explanation of how I am calculating the angle between two vectors here
    # to preserve the full 360 degrees of identifiability needed to divide into
    # four quadrants
    direction = np.array([dir_u, dir_v])
    direction_for_wacky_product = np.array([-1*dir_v, dir_u])

    dot = np.einsum("i,ijk", direction, uv_grid)
    det = np.einsum("i,ijk", direction_for_wacky_product, uv_grid)    
    angles = np.arctan2(det, dot)

    # Create array containing the quadrant number for each (lon, lat) pair. 
    # Numbering is 1, 2, 3, 4 counterclockwise starting from the given direction
    conditions = [angles >= np.pi/2, 
                  (angles < np.pi/2) & (angles >= 0), 
                  (angles < 0) & (angles >= -np.pi/2),
                  angles < -np.pi/2]
    choices = [2, 1, 4, 3]

    quadrants = np.select(conditions, choices, default = np.nan)

    # Update the DataArray and return 
    x = x.assign_coords({"quadrant": (("lat", "lon"), quadrants)})
    x.attrs["quadrant_direction"] = (dir_u, dir_v)

    return x


y = get_shear_stamp(2019, 1, 10, 18, 
                42.756, -105.012831, 800, 
                credentials.RDA_USER, credentials.RDA_PASSWORD)

quadrantize(y, 1, -1)
