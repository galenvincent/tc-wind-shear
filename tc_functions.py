import credentials
import GeoLocation as geo
import kernels as k
import RemoveVortexFunction as rm_vort

from siphon.catalog import TDSCatalog
from siphon.http_util import session_manager
from datetime import datetime
import numpy as np
import xarray as xr 
import math
import matplotlib.pyplot as plt

def gfs_access(year, month, day, hour, username, password):
    """
    Return a remote access object for a GFS dataset.

    Parameters:
        - year (int): Year of requested data access (2015 - present).
        - month (int): Month of requested data access. 
        - day (int): Day of requested data access.
        - hour (int): Hour of analysis for dataset. One of [0, 6, 12, 18].
        - username (str): RDA account username. 
        - password (str): RDA account password. 

    Returns:
        An object returned from remote_access function from the siphon module. 
    """

    # Convert time info into datetime object for easy maneuvering
    assert year >= 2015, "Year must be >= 2015"
    assert month >= 1 and month <= 12, "Month must be in [1, 12]"
    assert hour in [0, 6, 12, 18], "Hour must be one of 0, 6, 12, 18"
    ymdh = datetime(year, month, day, hour)

    # Set catalog URL and dataset name from time info
    catalog_url = ymdh.strftime("https://rda.ucar.edu/thredds/catalog/files/g/ds084.1/%Y/%Y%m%d/catalog.xml")
    dataset_name = ymdh.strftime("gfs.0p25.%Y%m%d%H.f000.grib2")

    # Set RDA username and password
    session_manager.set_session_options(auth = (username, password))

    # Get remote access to data
    catalog = TDSCatalog(catalog_url)
    ds = catalog.datasets[dataset_name]
    dataset = ds.remote_access()

    return dataset

def shear_stamp(lat, lon, radius, dataset, vortex_rm = False, vortex_rm_rad = None):
    """
    Retrieves a circular stamp of wind shear from GFS analysis (not forecast) 
    data associated with input parameter specifications. 

    Parameters:
        - lat (float): latitude for the center of the stamp (-90, 90)
        - lon (float):  longitude for the center of the stamp (-180, 180).
        - radius (float): radius of the circular stamp (km).
        - dataset: A GFS remote dataset as returned by the gfs_access()
        function.
        - vortex_rm (bool): Do you want to remove the hurricane vortex from 
        the wind fields?
        - vortex_rm_rad (float): If vortex_rm = True, what radius do you want
        to remove the vortex over? Given in km.

    Returns:
        (xarray.DataArray) A xarray DataArray object containing the zonal (u)/ 
        meridional (v) directions and magnitude for a circular stamp of vertical
        wind shear calculated between the 200 and 850 hPa isobars. 

    Additional Details:
        Use the sectorize() function to split up the wind shear data into 
        sectors from a specified direction.
    """
    
    assert radius > 0, "radius must be positive"
    assert lat >= -90 and lat <= 90, "Latitude must be in range (-90, 90)."
    assert lon >= -180 and lon <= 180, "Longitude must be in range (-180, 180)."

    # Pull zonal and meridional wind data at isobaric surfaces
    zonal_isobaric = dataset.variables["u-component_of_wind_isobaric"]
    merid_isobaric = dataset.variables["v-component_of_wind_isobaric"]

    iso_number = zonal_isobaric.dimensions[1]
    assert iso_number == merid_isobaric.dimensions[1], "Zonal and meridional wind vectors have different pressure sets in data."

    pressure_levels = dataset.variables[iso_number][:]
    data_lat = dataset.variables["lat"][:]
    data_lon = dataset.variables["lon"][:]

    # Get 200 & 850 hPa indices
    pressure_850_ind = np.argmin(np.abs(pressure_levels - 850 * 100))
    pressure_200_ind = np.argmin(np.abs(pressure_levels - 200 * 100))

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth)
    min_lat, max_lat, min_lon, max_lon = bounding_box(lat, lon, radius)

    # Get indices for the bounding lat
    lat_rows = np.where((data_lat >= min_lat) & (data_lat <= max_lat))[0]
    min_lat_ind = min(lat_rows) if min(lat_rows) == 0 else min(lat_rows) - 1
    max_lat_ind = max(lat_rows) if max(lat_rows) > (data_lat.size - 2) else max(lat_rows) + 2 
    lat_subset = data_lat[min_lat_ind:max_lat_ind]
    
    # Handle the case where longitude spans the transition from 0 to 360.
    if min_lon < 90 and max_lon > 270:
        lon_rows_1 = np.where((data_lon >= max_lon))[0]
        lon_rows_2 = np.where((data_lon <= min_lon))[0]
        
        min_lon_ind_1 = min(lon_rows_1) - 1 
        max_lon_ind_1 = max(lon_rows_1)
        min_lon_ind_2 = min(lon_rows_2)
        max_lon_ind_2 = max(lon_rows_2) + 2

        u_200_1 = zonal_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind_1:max_lon_ind_1]
        u_200_2 = zonal_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind_2:max_lon_ind_2]
        u_200 = np.concatenate((u_200_1, u_200_2), axis = 1)

        u_850_1 = zonal_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind_1:max_lon_ind_1]
        u_850_2 = zonal_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind_2:max_lon_ind_2]
        u_850 = np.concatenate((u_850_1, u_850_2), axis = 1)

        v_200_1 = merid_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind_1:max_lon_ind_1]
        v_200_2 = merid_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind_2:max_lon_ind_2]
        v_200 = np.concatenate((v_200_1, v_200_2), axis = 1)

        v_850_1 = merid_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind_1:max_lon_ind_1]
        v_850_2 = merid_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind_2:max_lon_ind_2]
        v_850 = np.concatenate((v_850_1, v_850_2), axis = 1)
        
        lon_subset_1 = data_lon[min_lon_ind_1:max_lon_ind_1]
        lon_subset_2 = data_lon[min_lon_ind_2:max_lon_ind_2]
        lon_subset = np.concatenate((lon_subset_1, lon_subset_2), axis = 0)
    else:
        lon_rows = np.where((data_lon >= min_lon) & (data_lon <= max_lon))[0]
    
        min_lon_ind = min(lon_rows) if min(lon_rows) == 0 else min(lon_rows) - 1
        max_lon_ind = max(lon_rows) if max(lon_rows) > (data_lon.size - 2) else max(lon_rows) + 2

        # Get the wind data
        u_200 = zonal_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
        u_850 = zonal_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
        v_200 = merid_isobaric[0, pressure_200_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
        v_850 = merid_isobaric[0, pressure_850_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]

        lon_subset = data_lon[min_lon_ind:max_lon_ind]

    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat_subset, lon_subset)]

    # Subset out a circular stamp of the data, putting nan everywhere else
    # Distance from each point to the center
    dist_mat = great_circ_dist(lat, lon, lat_grid, lon_grid)
    dist_tf = dist_mat > radius

    ### TODO: Determine the radius over which to remove the vortex in a data-driven way

    # Remove the hurricane vortex if desired
    if vortex_rm:
        u_200, v_200 = rm_vort.removeTCvortex(u_200, v_200, dist_mat, vortex_rm_rad)
        u_850, v_850 = rm_vort.removeTCvortex(u_850, v_850, dist_mat, vortex_rm_rad)

    # Set all other points to nan
    u_200_stamp = np.where(dist_tf, np.nan, u_200)
    u_850_stamp = np.where(dist_tf, np.nan, u_850)
    v_200_stamp = np.where(dist_tf, np.nan, v_200)
    v_850_stamp = np.where(dist_tf, np.nan, v_850)

    # Calculate wind shear components and magnitude
    u_wind_shear = u_200_stamp - u_850_stamp
    v_wind_shear = v_200_stamp - v_850_stamp
    magnitude_wind_shear = np.sqrt(np.power(u_wind_shear, 2) + np.power(v_wind_shear, 2))

    # Calculate the area-averaged wind shear vector
    u_mean = np.nanmean(u_wind_shear)
    v_mean = np.nanmean(v_wind_shear)
    mag_mean = np.sqrt(np.power(u_mean, 2) + np.power(v_mean, 2))

    # Convert center longitude to (0, 360) scale if in (-180, 180) scale:
    lon = lon + 360 if lon < 0 else lon

    shear = xr.DataArray(
        data = np.array([magnitude_wind_shear, u_wind_shear, v_wind_shear]),
        coords = {
            "component": ["magnitude", "u", "v"],
            "lat": lat_subset,
            "lon": lon_subset
        },
        dims = ["component", "lat", "lon"],
        name = "200-850hPa_wind_shear",
        attrs = {
            "long_name": "200_850 hPa vertical wind shear",
            "units": "meters/second",
            "center_lat": lat,
            "center_lon": lon,
            "stamp_radius": radius,
            "avg_shear": (u_mean, v_mean),
            "avg_magnitude": mag_mean
        }
    )

    return shear

def wind_stamp(lat, lon, radius, pressure, dataset, vortex_rm = False, vortex_rm_rad = None):
    """
    Retrieves a circular stamp of wind shear from GFS analysis (not forecast) 
    data associated with input parameter specifications. 

    Parameters:
        - lat (float): latitude for the center of the stamp (-90, 90)
        - lon (float):  longitude for the center of the stamp (-180, 180).
        - radius (float): radius of the circular stamp (km).
        - pressure (float): pressure level for wind data (hPa). Closest 
        available pressure level will be returned.
        - dataset: A GFS remote dataset as returned by the gfs_access()
        function. 
        - vortex_rm (bool): Do you want to remove the hurricane vortex from 
        the wind fields?
        - vortex_rm_rad (float): If vortex_rm = True, what radius do you want
        to remove the vortex over? Given in km.

    Returns:
        (xarray.DataArray) A xarray DataArray object containing the zonal (u)/ 
        meridional (v) directions and magnitude for a circular stamp of GFS
        isobaric wind data at the specified time, location, and radius. 

    Additional Details:
        Use the sectorize() function to split up the wind shear data into 
        sectors from a specified direction.
    """

    assert radius > 0, "radius must be positive"
    assert lat >= -90 and lat <= 90, "Latitude must be in range (-90, 90)."
    assert lon >= -180 and lon <= 180, "Longitude must be in range (-180, 180)."

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
    min_lat, max_lat, min_lon, max_lon = bounding_box(lat, lon, radius)

    # Get indices for the bounding lat/lon
    lat_rows = np.where((data_lat >= min_lat) & (data_lat <= max_lat))[0]

    min_lat_ind = min(lat_rows) if min(lat_rows) == 0 else min(lat_rows) - 1
    max_lat_ind = max(lat_rows) if max(lat_rows) > (data_lat.size - 2) else max(lat_rows) + 2 
    
    lat_subset = data_lat[min_lat_ind:max_lat_ind]

    # Handle the case where longitude spans the transition from 0 to 360.
    if min_lon < 90 and max_lon > 270:
        lon_rows_1 = np.where((data_lon >= max_lon))[0]
        lon_rows_2 = np.where((data_lon <= min_lon))[0]
        
        min_lon_ind_1 = min(lon_rows_1) - 1 
        max_lon_ind_1 = max(lon_rows_1)
        min_lon_ind_2 = min(lon_rows_2)
        max_lon_ind_2 = max(lon_rows_2) + 2

        u1 = zonal_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind_1:max_lon_ind_1]
        u2 = zonal_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind_2:max_lon_ind_2]
        u = np.concatenate((u1, u2), axis = 1)

        v1 = merid_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind_1:max_lon_ind_1]
        v2 = merid_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind_2:max_lon_ind_2]
        v = np.concatenate((v1, v2), axis = 1)
        
        lon_subset_1 = data_lon[min_lon_ind_1:max_lon_ind_1]
        lon_subset_2 = data_lon[min_lon_ind_2:max_lon_ind_2]
        lon_subset = np.concatenate((lon_subset_1, lon_subset_2), axis = 0)
    else:
        lon_rows = np.where((data_lon >= min_lon) & (data_lon <= max_lon))[0]
    
        min_lon_ind = min(lon_rows) if min(lon_rows) == 0 else min(lon_rows) - 1
        max_lon_ind = max(lon_rows) if max(lon_rows) > (data_lon.size - 2) else max(lon_rows) + 2

        # Get the wind data
        u = zonal_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]
        v = merid_isobaric[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]

        lon_subset = data_lon[min_lon_ind:max_lon_ind]

    # Get lat and lon grid corresponding to our data of interest 
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat_subset, lon_subset)]

    # Subset out a circular stamp of the data, putting nan everywhere else
    # Distance from each point to the center
    dist_mat = great_circ_dist(lat, lon, lat_grid, lon_grid)
    dist_tf = dist_mat > radius

    # Remove the hurricane vortex if desired
    if vortex_rm:
        u, v = rm_vort.removeTCvortex(u, v, dist_mat, vortex_rm_rad)

    # Set all other points to nan
    u_stamp = np.where(dist_tf, np.nan, u)
    v_stamp = np.where(dist_tf, np.nan, v)

    # Calculate wind magnitude
    magnitude = np.sqrt(np.power(u_stamp, 2) + np.power(v_stamp, 2))

    # Convert center longitude to (0, 360) scale if in (-180, 180) scale:
    lon = lon + 360 if lon < 0 else lon

    wind = xr.DataArray(
        data = np.array([magnitude, u_stamp, v_stamp]),
        coords = {
            "component": ["magnitude", "u", "v"],
            "lat": lat_subset,
            "lon": lon_subset
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
        }
    )

    return wind

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

def bounding_box(center_lat, center_lon, radius):
    """
    Get bounding box of longitude and latitude that fits a certain radius 
    within it. Mostly based on code in the GeoLocation module, which was 
    taken from:
        https://github.com/jfein/PyGeoTools/blob/master/geolocation.py

    Parameters:
        - center_lat (float): Latitude for the center of the box (-90, 90).
        - center_lon (float): Longitude for the center of the box (-180, 180).
        - radius (float): Radius to be held within the box (km).

    Returns:
        (list): A list of [minimum latitude, maximum latitude, minimum longitude,
        maximum longitude] for the bounding box, where longitude is now on
        (0, 360) to mesh with the GFS data well.
    """
    center = geo.GeoLocation.from_degrees(center_lat, center_lon)
    sw_loc, ne_loc = center.bounding_locations(radius)

    sw_lat = sw_loc.deg_lat
    sw_lon = sw_loc.deg_lon
    ne_lat = ne_loc.deg_lat
    ne_lon = ne_loc.deg_lon

    # Convert longitude to (0, 360) range to mesh with GFS data
    if sw_lon < 0:
        sw_lon = sw_lon + 360
    if ne_lon < 0:
        ne_lon = ne_lon + 360

    min_lat, max_lat = sorted([sw_lat, ne_lat])
    min_lon, max_lon = sorted([sw_lon, ne_lon])
        
    return [min_lat, max_lat, min_lon, max_lon]

def sectorize(x, dir_u, dir_v, name, n_sector = 4):
    """
    Determines sectors based on an initial direction vector.

    Takes in an xarray DataArray object (x) with lat/lon coordinates and
    an attribute indicating the center latitude and longitude. Then computes
    sectors based on some direction vector given by (dir_u, dir_v), where
    the sectors are labeled increasing in a counter-clockwise direction
    starting at (dir_u, dir_v). These sector labels are added to the DataArray
    and returned.

    Parameters:
        - x (xarray.DataArray): Wind shear or wind data output from the 
        shear_stamp() or wind_stamp() functions. 
        - dir_u (float): zonal (east) component of the direction vector.
        - dir_v (float): meridional (north) component of the direction vector.
        - name (str): what to name the layer of data that will contain the 
        sector numbers.
        - n_sector (int): Number of sectors you want to split the stamp into.
        Defalut is 4 (splitting into quadrants). This number must be a power 
        of 2 (i.e. 2, 4, 8, etc.).
    
    Returns:
        (xarray.DataArray): An updated DataArray object with a new coordinate to 
        indicate which sector each data point is in.
    """
    
    assert dir_u != 0 or dir_v != 0, "One of dir_u or dir_v must be non-zero."
    assert isinstance(x, xr.core.dataarray.DataArray), "x must be a DataArray."
    assert np.floor(np.log2(n_sector)) == np.ceil(np.log2(n_sector)), "n_sector must be a power of 2"

    # Get latitude and longitude arrays from the DataArray object
    lat = x.lat.values
    lon = x.lon.values

    # Create lat/lon grids out of the arrays, so we have every individual 
    # (lon, lat) coordinate combination
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat, lon)]    

    # Get the center of the stamp from the DataArray
    center_lat = x.attrs["center_lat"]
    center_lon = x.attrs["center_lon"]

    # Convert to (0, 360) scaled longitude if in (-180, 180) scale.
    center_lon = center_lon + 360 if center_lon < 0 else center_lon

    # Create centered lat/lon data and put into one array
    uv_grid = np.array([lon_grid - center_lon, lat_grid - center_lat])

    # Compute angle between each point in grid and the direction given in 
    # function input. See https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors/16544330#16544330
    # for explanation of how I am calculating the angle between two vectors here
    # to preserve the full 360 degrees of identifiability needed to divide into
    # n_sector sectors
    direction = np.array([dir_u, dir_v])
    direction_for_wacky_product = np.array([-1*dir_v, dir_u])

    dot = np.einsum("i,ijk", direction, uv_grid)
    det = np.einsum("i,ijk", direction_for_wacky_product, uv_grid)    
    angles = np.arctan2(det, dot)

    # Create array containing the sector number for each (lon, lat) pair. 
    # Numbering is increasing counterclockwise starting from the given direction.
    angle_breaks = np.linspace(np.pi, -np.pi, n_sector+1)
    conditions = []
    for ii in range(0, angle_breaks.size - 1):
        conditions.append( (angles <= angle_breaks[ii]) & (angles > angle_breaks[ii + 1]) )

    # Draw a picture to see why the numbering scheme works like this
    choices = np.concatenate((np.arange(n_sector/2, 0, -1), np.arange(n_sector, n_sector/2, -1)))

    sectors = np.select(conditions, choices, default = 2)

    # Update the DataArray and return 
    x = x.assign_coords({"sector_"+name: (("lat", "lon"), sectors)})
    x.attrs["sector_"+name+"_direction"] = (dir_u, dir_v)

    return x

def pressure_min(seed_lat, seed_lon, search_radius, dataset):
    """
    Determines the center of a TC based on surface-level pressure minimum
    within some search radius of the seed latitude and longitude.

    Parameters:
        - seed_lat (float): Latitude of seed searching point.
        - seed_lon (float): Longitude of seed searching point. 
        - search_radius (float): The distance (km) you would like to search in 
        all directions for the pressure minimum. Only pressure within this 
        radius will be evaluated. 
        - dataset: A GFS remote dataset as returned by the gfs_access()
        function.
    
    Returns:
        (numpy.ndarray): An array of [lat, lon] of the pressure minimum 
        latitude: (-90, 90), longitude: (-180, 180).
    """

    assert search_radius > 0, "radius must be positive"
    assert seed_lat >= -90 and seed_lat <= 90, "Latitude must be in range (-90, 90)."
    assert seed_lon >= -180 and seed_lon <= 180, "Longitude must be in range (-180, 180)."

    # Pull pressure at surface level
    if "Pressure_surface" in dataset.variables:
        pressure_ds = dataset.variables["Pressure_surface"]
    else:
        return [seed_lat, seed_lon]

    lat = dataset.variables["lat"][:]
    lon = dataset.variables["lon"][:]

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth, hence the fancy code in bounding_box)
    min_lat, max_lat, min_lon, max_lon = bounding_box(seed_lat, seed_lon, search_radius)

    # Get indices for the bounding lat/lon
    lat_rows = np.where((lat >= min_lat) & (lat <= max_lat))[0]
    lon_rows = np.where((lon >= min_lon) & (lon <= max_lon))[0]
    min_lat_ind, max_lat_ind = [min(lat_rows) - 1, max(lat_rows) + 1]
    min_lon_ind, max_lon_ind = [min(lon_rows) - 1, max(lon_rows) + 1]

    # Get the pressure data from within the bounding box
    pressure = pressure_ds[0, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]

    lat_subset = lat[min_lat_ind:max_lat_ind]
    lon_subset = lon[min_lon_ind:max_lon_ind]
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat_subset, lon_subset)]

    # Subset out a circular stamp of the data
    # Distance from each point to the center
    dist_mat = great_circ_dist(seed_lat, seed_lon, lat_grid, lon_grid)
    dist_tf = dist_mat > search_radius
    # Make a masked array based on this radius.
    pressure_stamp = np.ma.array(pressure, mask = dist_tf)
    min_ind = np.unravel_index(np.argmin(pressure_stamp), pressure_stamp.shape)

    lat_min = lat_grid[min_ind]
    lon_min = lon_grid[min_ind]

    lon_min = lon_min - 360 if lon_min > 180 else lon_min

    return [lat_min, lon_min]

def vorticity_centroid(seed_lat, seed_lon, pressure, search_radius, calc_radius, dataset):
    """
    Determines the center of a TC based on the isobaric vorticity centroid at
    some specific pressure. Calculations done within some calculation radius of 
    the seed latitude and longitude. Center of the TC is first initialized 
    from the seed location using minimum pressure. 

    Parameters:
        - seed_lat (float): Latitude of seed searching point.
        - seed_lon (float): Longitude of seed searching point. 
        - pressure (float): The isobar to evaluate vorticity at (hPa). Closest 
        available pressure level available in the data will be used.
        - search_radius (float): The distance (km) you would like to search in 
        all directions for the pressure minimum used to initialize the 
        vorticity calculation. 
        - calc_radius (float): The distanec (km) you would like to take the
        vorticity centroid over.
        - dataset: A GFS remote dataset as returned by the gfs_access()
        function.
    
    Returns:
        (numpy.ndarray): An array of [lat, lon] of the vorticity centroid.
    """

    assert search_radius > 0, "radius must be positive"
    assert seed_lat >= -90 and seed_lat <= 90, "Latitude must be in range (-90, 90)."
    assert seed_lon >= -180 and seed_lon <= 180, "Longitude must be in range (-180, 180)."

    # Pull the absolute vorticity dataset
    vort_ds = dataset.variables["Absolute_vorticity_isobaric"]

    iso_number = vort_ds.dimensions[1]
    pressure_levels = dataset.variables[iso_number][:]
    lat = dataset.variables["lat"][:]
    lon = dataset.variables["lon"][:]

    # Get closest pressure index
    pressure = pressure * 100 # convert from hPa to Pa
    pressure_ind = np.argmin(np.abs(pressure_levels - pressure))

    lat = dataset.variables["lat"][:]
    lon = dataset.variables["lon"][:]

    pressure_lat, pressure_lon = pressure_min(seed_lat, seed_lon, search_radius, dataset)
    # Convert long to (-180, 180) for the sake of bounding box

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth, hence the fancy code in bounding_box)
    min_lat, max_lat, min_lon, max_lon = bounding_box(pressure_lat, pressure_lon, calc_radius)

    # Get indices for the bounding lat/lon
    lat_rows = np.where((lat >= min_lat) & (lat <= max_lat))[0]
    lon_rows = np.where((lon >= min_lon) & (lon <= max_lon))[0]
    min_lat_ind, max_lat_ind = [min(lat_rows) - 1, max(lat_rows) + 1]
    min_lon_ind, max_lon_ind = [min(lon_rows) - 1, max(lon_rows) + 1]

    # Get the pressure data from within the bounding box
    vort = vort_ds[0, pressure_ind, min_lat_ind:max_lat_ind, min_lon_ind:max_lon_ind]

    lat_subset = lat[min_lat_ind:max_lat_ind]
    lon_subset = lon[min_lon_ind:max_lon_ind]
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat_subset, lon_subset)]

    # Subset out a circular stamp of the data
    # Distance from each point to the center
    dist_mat = great_circ_dist(pressure_lat, pressure_lon, lat_grid, lon_grid)
    dist_tf = dist_mat > calc_radius
    # Make a masked array based on this radius.
    vort_stamp = np.ma.array(vort, mask = dist_tf)
    
    # Calculate centroid
    lat_min = np.sum(vort_stamp * lat_grid)/np.sum(vort_stamp)
    lon_min = np.sum(vort_stamp * lon_grid)/np.sum(vort_stamp)

    # Convert lon back to (-180 to 180)
    lon_min = lon_min - 360 if lon_min > 180 else lon_min

    return [lat_min, lon_min]

def radial_profile(x, stride = 10, h = 25, sector_labels = [], kernel = k.epanechnikov, normalized = False):
    """
    Compute the radial profiles of a stamp.

    Computes radial profiles for stamps returned from shear_stamp() or 
    wind_stamp() functions. Will compute separate profiles for however many 
    different sectors are encoded in the data (default is 1, but could be
    4 after applying the sectorize() function).
    
    A kernel with bandwidth h around each radius is used to calculate the radial
    profile, so each profile is the average shear/wind around r +- h weighted by 
    the kernel.

    Parameters:
        - x (xarray.DataArray): Wind shear or wind data output from the 
        shear_stamp() or wind_stamp() functions. Or after applying sectorize().
        - stride (float): Distance (km) between radii to calculate radial
        profiles at. 
        - h (float): Bandwidth (km) for the kernel used to calculate each 
        radial profile. 
        - sector_labels (list): List of names of the various sector coordinates 
        that you would like to compute over. If an empty list, no sectors will be 
        profiles will be computed over all data.
        - kernel (function): A function from kernels.py or any vectorized kernel 
        function whose first two arguments are (data, bandwidth).
    
    Returns:
        (xarray.DataArray): A data array containing radial profiles for each 
        sector for both the u (zonal) and v (meridional) directions. 

    """

    assert isinstance(x, xr.core.dataarray.DataArray), "x must be a DataArray."
    assert isinstance(sector_labels, list), "sector_labels must be a list."

    # Get latitude and longitude arrays from the DataArray object
    lat = x.lat.values
    lon = x.lon.values

    # Create lat/lon grids out of the arrays, so we have every individual 
    # (lon, lat) coordinate combination
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat, lon)]

    # Get the center of the stamp from the DataArray
    center_lat = x.attrs["center_lat"]
    center_lon = x.attrs["center_lon"]

    # Create distance matrix of km from the center of the stamp
    dist_mat = great_circ_dist(center_lat, center_lon, lat_grid, lon_grid)

    # Create list of radii that we will calculate the radial profile at
    radii = np.arange(start = math.ceil(h/2), 
                      stop = x.attrs["stamp_radius"] - math.ceil(h/2), 
                      step = stride)

    # Get the sector information from the DataArray
    #sectors = x.sector.values
    #sectors_unique = np.sort(np.unique(sectors))
    
    # Pull the zonal and meridional components
    u_comp = np.ma.masked_invalid(x.sel(component = "u").values)
    v_comp = np.ma.masked_invalid(x.sel(component = "v").values)
    stamp_mask = u_comp.mask

    # If normalized = True, then normalize the components by the average magnitude
    if normalized:
        mag = x.attrs["avg_magnitude"]
        u_comp = u_comp/mag
        v_comp = v_comp/mag

    if not sector_labels:
        sector_labels = ["placeholder"]
        unique_sectors = [1]
        n_unique_sectors = 1
    else:
        n_unique_sectors = []
        for label in sector_labels:
            sectors = x.coords["sector_"+label].values
            unique_sectors = np.sort(np.unique(sectors))
            n_unique_sectors.append(unique_sectors.size)
        
        assert n_unique_sectors.count(n_unique_sectors[0]) == len(n_unique_sectors), "All sectorizations must have the same number of sectors."

        n_unique_sectors = n_unique_sectors[0]

    n_sector_labels = len(sector_labels)

    # Initialize data structure for the profiles
    profiles = np.zeros((2, radii.size, n_sector_labels, n_unique_sectors))

    # Loop through each radius we will calcuate at
    for ii in range(radii.size):
        # Apply the kernel function at the desired radius to get weights
        dist_kerneled = np.ma.array(kernel(dist_mat - radii[ii], h))
        # Multiply weights by the u and v components
        u_weighted = np.multiply(dist_kerneled, u_comp)
        v_weighted = np.multiply(dist_kerneled, v_comp)

        # loop through each sector label
        for jj in range(n_sector_labels):
            if sector_labels[jj] == "placeholder":
                profiles[0, ii, jj, 0] = np.sum(u_weighted)/np.sum(dist_kerneled)
                profiles[1, ii, jj, 0] = np.sum(v_weighted)/np.sum(dist_kerneled)
                break

            sectors = x.coords["sector_"+sector_labels[jj]].values
            
            # loop through each sector
            for kk in range(n_unique_sectors):
                # Create a mask for the specific sector
                sector_mask = np.logical_or(stamp_mask, sectors != unique_sectors[kk])
                u_weighted.mask = sector_mask
                v_weighted.mask = sector_mask
                dist_kerneled.mask = sector_mask

                # Calculate the weighted average for each component in this sector 
                # at this radius value
                profiles[0, ii, jj, kk] = np.sum(u_weighted)/np.sum(dist_kerneled)
                profiles[1, ii, jj, kk] = np.sum(v_weighted)/np.sum(dist_kerneled)

    # Create the DataArray to return
    radial_prof = xr.DataArray(
        data = profiles,
        coords = {
            "component": ["u", "v"],
            "radius": radii,
            "sector_label": sector_labels,
            "sector": unique_sectors
        },
        dims = ["component", "radius", "sector_label", "sector"],
        name = "radial_profile",
        attrs = {
            "long_name": "Radial ORB function for wind shear or wind.",
            "center_lat": x.attrs["center_lat"],
            "center_lon": x.attrs["center_lon"],
            "stamp_radius": x.attrs["stamp_radius"]
        }
    )

    for label in sector_labels:
        if label == "placeholder": 
            break
        
        radial_prof.attrs["sector_"+label+"_direction"] = x.attrs["sector_"+label+"_direction"]
    
    return radial_prof

def storm_radius(center_lat, center_lon, dataset, max_radius = 800, pressure = 850, 
                 stride = 10, h = 25, kernel = k.epanechnikov, plot = False):
    """
    Calculate the radius of a storm to remove a vortex over using the method 
    detailed in this paper: https://doi.org/10.1175/1520-0493(1995)123<2791:IITGHP>2.0.CO;2.

    Parameters:
        - center_lat (float) Latitude of storm center on (-90, 90)
        - center_lon (float) Longitude of storm center on (-180, 180)
        - dataset: A GFS remote dataset as returned by the gfs_access()
        function.
        - max_radius (float): Maximum radius to look out to (km)
        - pressure (float): Pressure isobar to pull wind field from (hPa)
        - stride (float): Distance between radial steps for the tangental wind 
        profile (km).
        - h (float): Bandwidth for the kernel used to calculate the tangential 
        wind profile (km).
        - kernel (function): A function from kernels.py or any vectorized kernel 
        function whose first two arguments are (data, bandwidth).
        - plot (boolean): Whether or not to produce a plot of the tangential 
        widn profile with indicator marking the radius selected.

    Returns:
        (float) The radius selected for vortex removal (km).
    """

    assert max_radius > 0, "radius must be positive"
    assert center_lat >= -90 and center_lat <= 90, "Latitude must be in range (-90, 90)."
    assert center_lon >= -180 and center_lon <= 180, "Longitude must be in range (-180, 180)."

    # Get wind field at specified pressure
    wind = wind_stamp(center_lat, center_lon, max_radius, pressure, dataset)

    # Extract lat/lon grid
    lat = wind.lat.values
    lon = wind.lon.values
    lat_grid, lon_grid = [x.T for x in np.meshgrid(lat, lon)]
    # Important to use the center_lat and center_lon from the wind dataset, as it
    # will have the correct scaling (0, 360) to match the data from GFS
    center_lat = wind.attrs["center_lat"]
    center_lon = wind.attrs["center_lon"]
    dist_mat = great_circ_dist(center_lat, center_lon, lat_grid, lon_grid)

    # Get grid data poisitions relative to storm center
    centered_position_vectors = np.array([lon_grid - center_lon, lat_grid - center_lat])

    # Get wind vector data into a numpy array
    wind_vectors = np.ma.masked_invalid([wind.sel(component = 'u'), wind.sel(component = 'v')])

    def vector_rejection(a, b):
        """
        Calculate the vector rejection (i.e. perpendicular projection) of vector b onto vector a.
        
        Parameters:
            - a (3d array): an array of shape (2, #lon, #lat) containing vector components for different positions on a grid.
            - b (3d array): an array of the same shape as a, containing the vector components for projecting onto a.

        Returns:
            (2d array) an array of shape (2, #lon, #lat) containing the vector rejections.
        """
        dot_ab = np.einsum("ijk,ijk->jk", a, b)
        dot_aa = np.einsum("ijk,ijk->jk", a, a)
        to_subtract = np.multiply(np.divide(dot_ab, dot_aa), a)
        return np.subtract(b, to_subtract)

    # Get tangential wind component
    tangential_component = vector_rejection(centered_position_vectors, wind_vectors)

    # Calculate cross produce between grid position and tangential wind in order
    # to determine CCW (+) or CW (-) direction
    wind_vectors_for_cross = wind_vectors[[1,0],:,:]
    wind_vectors_for_cross[1,:,:] = -1 * wind_vectors_for_cross[1,:,:]
    cross_prod = np.einsum("ijk,ijk->jk", centered_position_vectors, wind_vectors_for_cross)
    cross_prod = np.ma.masked_invalid(cross_prod)

    # Convert tangential vectors to magnitude, with sign indicating CW (-) or CCW (+).
    tangential_magnitude = np.multiply(np.sqrt(np.power(tangential_component[0,:,:], 2) + np.power(tangential_component[1,:,:], 2)), np.sign(cross_prod))

    # Compute tangential wind profile
    radii = np.arange(start = 0 + math.ceil(h/2), stop = max_radius - math.ceil(h/2), step = stride)
    profile = np.zeros(radii.size)
    for ii in range(radii.size):
        # Apply the kernel function at the desired radius to get weights
        dist_kerneled = np.ma.array(kernel(dist_mat - radii[ii], h))
        # Multiply weights by the u and v components
        weighted = np.multiply(dist_kerneled, tangential_magnitude)
        profile[ii] = np.nansum(weighted)/np.sum(dist_kerneled)
        if np.isnan(profile[ii]):
            profile[ii] = 0


    # Numerical derivative of profile
    deriv = np.gradient(profile, stride)

    # Find peak of the profile (which should be just at the eyewall, or nearby)
    max_wind_ind = np.argmax(profile)

    # Handle the radius selection based on conditions outlined in the paper 
    # given in the function description.
    condition1 = np.where(((profile < 6) & (-1*deriv < 4e-3)), 1, 0)
    condition1[0:max_wind_ind] = 0

    # Search for first place where 3 hits in a row: 
    b = np.array([1, 1, 1])
    three_conseq = np.array([x for x in range(condition1.size) if np.array_equal(condition1[x:x+b.size], b)])
    if three_conseq.size != 0:
        check1 = three_conseq[0] + 2
    else:
        check1 = np.inf

    # Search for where wind speed drops below 3 m/s
    condition2 = np.where(profile < 3, 1, 0)
    condition2[0:max_wind_ind] = 0

    low_wind = np.nonzero(condition2)[0]
    if low_wind.size != 0:
        check2 = low_wind[0]
    else: 
        check2 = np.inf

    rf_ind = np.minimum(check1, check2)
    if np.isinf(rf_ind):
        rf_ind = radii.size - 1
    
    rf_ind = int(rf_ind)
    rf = radii[rf_ind]

    # Make a pretty plot 
    if plot:
        plt.plot(radii, profile)
        plt.plot(rf, profile[rf_ind], "ob")
        plt.plot(radii, np.zeros(radii.size), linewidth = 0.5, alpha = 0.5, color = "black")
        plt.plot(radii, 3*np.ones(radii.size), linewidth = 0.5, alpha = 0.5, color = "black")
        plt.plot(radii, 6*np.ones(radii.size), linewidth = 0.5, alpha = 0.5, color = "black")
        plt.ylabel("tangential velocity (m/s)")
        plt.xlabel("radius (km)")
    
    return rf

    