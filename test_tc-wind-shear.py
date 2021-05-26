# Unit tests file for the tc-wind-shear project
import pytest
import numpy as np

import tc_functions as fun
# You will need to create your own credentials file containing two variables, 
# RDA_USER = your RDA usernameS
# RDA_PASSWORD = your password
import credentials

def test_shear_stamp():
    # Test that shear_stamp outputs appropriate data based on stable 
    # build on 1/22/2021

    dataset = fun.gfs_access(2020, 11, 3, 0, credentials.RDA_USER, credentials.RDA_PASSWORD)
    data = fun.shear_stamp(14.3, 277.5 - 360, 800, dataset)
    
    correct_lat = np.array(
        [21.5 , 21.25, 21.  , 20.75, 20.5 , 20.25, 20.  , 19.75, 19.5 ,
        19.25, 19.  , 18.75, 18.5 , 18.25, 18.  , 17.75, 17.5 , 17.25,
        17.  , 16.75, 16.5 , 16.25, 16.  , 15.75, 15.5 , 15.25, 15.  ,
        14.75, 14.5 , 14.25, 14.  , 13.75, 13.5 , 13.25, 13.  , 12.75,
        12.5 , 12.25, 12.  , 11.75, 11.5 , 11.25, 11.  , 10.75, 10.5 ,
        10.25, 10.  ,  9.75,  9.5 ,  9.25,  9.  ,  8.75,  8.5 ,  8.25,
        8.  ,  7.75,  7.5 ,  7.25]
    )
    correct_lon = np.array(
        [270.  , 270.25, 270.5 , 270.75, 271.  , 271.25, 271.5 , 271.75,
        272.  , 272.25, 272.5 , 272.75, 273.  , 273.25, 273.5 , 273.75,
        274.  , 274.25, 274.5 , 274.75, 275.  , 275.25, 275.5 , 275.75,
        276.  , 276.25, 276.5 , 276.75, 277.  , 277.25, 277.5 , 277.75,
        278.  , 278.25, 278.5 , 278.75, 279.  , 279.25, 279.5 , 279.75,
        280.  , 280.25, 280.5 , 280.75, 281.  , 281.25, 281.5 , 281.75,
        282.  , 282.25, 282.5 , 282.75, 283.  , 283.25, 283.5 , 283.75,
        284.  , 284.25, 284.5 , 284.75]
    )

    correct_mag = np.load("test_files/get_shear_stamp_test_magnitude.npy")

    assert np.array_equal(data.lat.values, correct_lat), "Latitude points do not match."
    assert np.array_equal(data.lon.values, correct_lon), "Longitude points do not match."
    assert np.array_equal(data[0].values, correct_mag, equal_nan=True), "Shear magnitude does not match."
    assert data.attrs["center_lat"] == 14.3, "Center latitude does not match."
    assert data.attrs["center_lon"] == 277.5, "Center longitude does not match."
    assert data.attrs["stamp_radius"] == 800, "Stamp radius does not match."

    # Check that magnitude is actually given by the correct combination of components. 
    u = data.sel(component = "u").values
    v = data.sel(component = "v").values
    mag = np.sqrt(np.power(u, 2) + np.power(v, 2))
    assert np.array_equal(mag, correct_mag, equal_nan=True), "Magnitude doesn't match components."
    
    # Check that the stamp size is around how big you need it to be
    lat_grid, lon_grid = [x.T for x in np.meshgrid(data.lat.values, data.lon.values)]
    dist_mat = fun.great_circ_dist(data.attrs["center_lat"], data.attrs["center_lon"], lat_grid, lon_grid)
    dist_mat = np.ma.array(dist_mat, mask = np.isnan(mag))
    assert np.max(dist_mat) < 800, "Stamp is too large."
    assert np.max(dist_mat) > 750, "Stamp is too small."

def test_wind_stamp():
    # Test that wind_stamp outputs appropriate data based on stable 
    # build on 1/22/2021

    dataset = fun.gfs_access(2020, 11, 3, 0, credentials.RDA_USER, credentials.RDA_PASSWORD)
    data = fun.wind_stamp(14.3, 277.5 - 360, 800, 450, dataset)
    
    correct_lat = np.array(
        [21.5 , 21.25, 21.  , 20.75, 20.5 , 20.25, 20.  , 19.75, 19.5 ,
        19.25, 19.  , 18.75, 18.5 , 18.25, 18.  , 17.75, 17.5 , 17.25,
        17.  , 16.75, 16.5 , 16.25, 16.  , 15.75, 15.5 , 15.25, 15.  ,
        14.75, 14.5 , 14.25, 14.  , 13.75, 13.5 , 13.25, 13.  , 12.75,
        12.5 , 12.25, 12.  , 11.75, 11.5 , 11.25, 11.  , 10.75, 10.5 ,
        10.25, 10.  ,  9.75,  9.5 ,  9.25,  9.  ,  8.75,  8.5 ,  8.25,
        8.  ,  7.75,  7.5 ,  7.25]
    )
    correct_lon = np.array(
        [270.  , 270.25, 270.5 , 270.75, 271.  , 271.25, 271.5 , 271.75,
        272.  , 272.25, 272.5 , 272.75, 273.  , 273.25, 273.5 , 273.75,
        274.  , 274.25, 274.5 , 274.75, 275.  , 275.25, 275.5 , 275.75,
        276.  , 276.25, 276.5 , 276.75, 277.  , 277.25, 277.5 , 277.75,
        278.  , 278.25, 278.5 , 278.75, 279.  , 279.25, 279.5 , 279.75,
        280.  , 280.25, 280.5 , 280.75, 281.  , 281.25, 281.5 , 281.75,
        282.  , 282.25, 282.5 , 282.75, 283.  , 283.25, 283.5 , 283.75,
        284.  , 284.25, 284.5 , 284.75]
    )

    correct_mag = np.load("test_files/get_wind_stamp_test_magnitude.npy")

    assert np.array_equal(data.lat.values, correct_lat), "Latitude points do not match."
    assert np.array_equal(data.lon.values, correct_lon), "Longitude points do not match."
    assert np.array_equal(data[0].values, correct_mag, equal_nan=True), "Shear magnitude does not match."
    assert data.attrs["center_lat"] == 14.3, "Center latitude does not match."
    assert data.attrs["center_lon"] == 277.5, "Center longitude does not match."
    assert data.attrs["stamp_radius"] == 800, "Stamp radius does not match."
    assert data.attrs["pressure_level"] == 450.0

    # Check that magnitude is actually given by the correct combination of components. 
    u = data.sel(component = "u").values
    v = data.sel(component = "v").values
    mag = np.sqrt(np.power(u, 2) + np.power(v, 2))
    assert np.array_equal(mag, correct_mag, equal_nan=True), "Magnitude doesn't match components."
    
    # Check that the stamp size is around how big you need it to be
    lat_grid, lon_grid = [x.T for x in np.meshgrid(data.lat.values, data.lon.values)]
    dist_mat = fun.great_circ_dist(data.attrs["center_lat"], data.attrs["center_lon"], lat_grid, lon_grid)
    dist_mat = np.ma.array(dist_mat, mask = np.isnan(mag))
    assert np.max(dist_mat) < 800, "Stamp is too large."
    assert np.max(dist_mat) > 750, "Stamp is too small."
    

