import credentials
import bounding_box as bb

from siphon.catalog import TDSCatalog
from siphon.http_util import session_manager
from datetime import datetime
import numpy as np
import xarray as xr 

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
        - lat/lon (float): lattitude & longitude for the center of stamp (degrees).
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
    lat = dataset.variables["lat"][:]
    lon = dataset.variables["lon"][:]

    # Get 200 & 850 hPa indices
    pressure_200_ind = np.where(pressure_levels == 20000)[0][0]
    pressure_850_ind = np.where(pressure_levels == 85000)[0][0]

    # Get min and max lat/lon based on a bounding box around the circle of 
    # interest. (Note that the radius of this circle is actually the distance
    # along the great circle of the earth, hence the fancy code in bounding_box)

    center = bb.GeoLocation.from_degrees(lat, lon)
    min_lat, max_lat, min_lon, max_lon = center.bounding_locations(radius)









