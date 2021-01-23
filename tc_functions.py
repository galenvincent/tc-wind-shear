import credentials
import GeoLocation as geo

from siphon.catalog import TDSCatalog
from siphon.http_util import session_manager
from datetime import datetime
import numpy as np
import xarray as xr 
import math

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

def shear_stamp(lat, lon, radius, dataset):
    """
    Retrieves a circular stamp of wind shear from GFS analysis (not forecast) 
    data associated with input parameter specifications. 

    Parameters:
        - lat (float): latitude for the center of the stamp (-90, 90)
        - lon (float):  longitude for the center of the stamp (-180, 180).
        - radius (float): radius of the circular stamp (km).
        - dataset: A GFS remote dataset as returned by the gfs_access()
        function.

    Returns:
        (DataArray) A xarray DataArray object containing the zonal (u)/ 
        meridional (v) directions and magnitude for a circular stamp of vertical
        wind shear calculated between the 200 and 850 hPa isobars. 

    Additional Details:
        Use the quadrantize() function to split up the wind shear data into 
        quadrants from a specified direction.
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
    pressure_200_ind = np.where(pressure_levels == 20000)[0][0]
    pressure_850_ind = np.where(pressure_levels == 85000)[0][0]

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth)
    min_lat, max_lat, min_lon, max_lon = bounding_box(lat, lon, radius)

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

    # Convert center longitude to (0, 360) scale if in (-180, 180) scale:
    lon += 360 if lon < 0 else lon

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
        }
    )

    return shear

def wind_stamp(lat, lon, radius, pressure, dataset):
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

    Returns:
        (DataArray) A xarray DataArray object containing the zonal (u)/ 
        meridional (v) directions and magnitude for a circular stamp of GFS
        isobaric wind data at the specified time, location, and radius. 

    Additional Details:
        Use the quadrantize() function to split up the wind shear data into 
        quadrants from a specified direction.
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

    # Convert center longitude to (0, 360) scale if in (-180, 180) scale:
    lon += 360 if lon < 0 else lon

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
        sw_lon += 360
    if ne_lon < 0:
        ne_lon += 360

    min_lat, max_lat = sorted([sw_lat, ne_lat])
    min_lon, max_lon = sorted([sw_lon, ne_lon])
        
    return [min_lat, max_lat, min_lon, max_lon]

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
        - x (DataArray): Wind shear or wind data output from the shear_stamp()
        or wind_stamp() functions. 
        - dir_u (float): zonal (east) component of the direction vector.
        - dir_v (float): meridional (north) component of the direction vector.
    
    Returns:
        (DataArray): An updated DataArray object with a new coordinate to 
        indicate which quadrant each data point is in.
    """
    
    assert dir_u != 0 or dir_v != 0, "One of dir_u or dir_v must be non-zero."
    assert isinstance(x, xr.core.dataarray.DataArray), "x must be a DataArray."

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
    center_lon += 360 if center_lon < 0 else center_lon

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
    pressure_ds = dataset.variables["Pressure_surface"]

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

    lon_min -= 360 if lon_min > 180 else lon_min

    return [lat_min, lon_min]

def vorticity_centroid(seed_lat, seed_lon, pressure, search_radius, calc_radius, dataset):
    """
    Determines the center of a TC based on the isobaric vorticity centroid at
    some specific pressure. Calculations done within some calculation radius of 
    the seed latitude and longitude. Center of the TC is firrst initialized 
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

    return [lat_min, lon_min]

# Ideas for the future:
#   - Figure out how to calculate the radius of maximal wind (RMW) and then use 
#   this papers (https://www.mdpi.com/2073-4433/10/7/376/htm) suggestion of 2.25
#   RMW as the distance to calculate vorticity over. 
#   - You might need a very similar function to calculate RMW that you will need
#   to calculate radial ORB! so that can come in handy. 
    
ds = gfs_access(2020, 11, 3, 0, credentials.RDA_USER, credentials.RDA_PASSWORD)
y = shear_stamp(14.3, 277.5 - 360, 800, ds)
quadrantize(y, 1, -1)
