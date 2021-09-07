import numpy as np
import pandas as pd
import io
import urllib.request as url
from datetime import datetime, timedelta

class Hurdat:
    """
    Class for interfacing with hurdat2 data.

    Initialization:
        - data (Pandas DataFrame): A dataframe with columns giving information 
        on the storms of interest. Columns should include at least "ID", 
        "DATETIME", "LAT", "LON", "WIND", and "DISTANCE"
    """
    def __init__(self, data):
        """
        Import Hurdat2 dataset.
        """

        # Make sure at least the neccesary columns are present
        assert "ID" in data.columns
        assert "DATETIME" in data.columns
        assert "LAT" in data.columns
        assert "LON" in data.columns
        assert "WIND" in data.columns
        assert 'DISTANCE' in data.columns

        self.storms = data

    def get_tc_by_id(self, id):
        """
        Extract the Hurdat2 entry for a single TC by its unique ID.

        Parameters:
            - id: Hurdat ID for TC, e.g. AL152016 for Hurricane Nicole [2016] 
            (15th TC in the North Atlantic in the 2016 season).
        
        Returns:
            Pandas dataframe containing full Hurdat2 data for the TC of interest.
        """
        return self.storms[self.storms['ID'].str.match(id)]

    def linear_interp(self, id, time):
        """
        Get estimated Hurdat2 data for TC at non-synoptic times.

        The Hurdat2 database records best track data for TCs at 6 hour intervals. This function provides a linear
        interpolation for TC location and intensity.

        :param id: Hurdat ID for TC, e.g. AL152016 for Hurricane Nicole [2016] (15th TC in the North Atlantic in the
                   2016 season.
        :param time: Time at which TC data is desired.
        :return: Data frame row containing interpolated TC data.
        """
        storm = self.get_tc_by_id(id)
        last_ts = storm[storm['DATETIME'] <= time].iloc[[-1]]
        if last_ts.iloc[0]['DATETIME'] == time:
            return(last_ts)
        next_ts = storm[storm['DATETIME'] >= time].iloc[[0]]
        if next_ts.iloc[0]['DATETIME'] == time:
            return(next_ts)
        fraction = (
                np.abs(time - last_ts.iloc[0]['DATETIME'])/
                (
                    np.abs(next_ts.iloc[0]['DATETIME'] - time) +
                    np.abs(time - last_ts.iloc[0]['DATETIME'])
                )
        )
        out = last_ts
        for col in ['LAT', 'LON', 'WIND']:
            dif = fraction*(next_ts.iloc[0][col] - last_ts.iloc[0][col])
            out[col] = last_ts.iloc[0][col] + dif

        out['DATETIME'] = time
        out['DATE'] = time.strftime(format='%Y%m%d')
        out['TIME'] = time.strftime(format='%H%M')

        return(out)

    def get_hourly_position(self, id):
        """
        Interpolates Hurdat2 data for a TC from 6-hourly to hourly data.

        The Hurdat2 database records best track data for TCs at 6 hour intervals. Uses Hurdat.linear_interp() to
        estimate hourly track data.

        :param id: Hurdat ID for TC, e.g. AL152016 for Hurricane Nicole [2016] (15th TC in the North Atlantic in the
                   2016 season.
        :return: Data frame containing hourly interpolation of TC data.
        """
        storm = self.get_tc_by_id(id)
        startDate = np.min(storm['DATETIME'])
        endDate = np.max(storm['DATETIME'])
        hourly_position = self.linear_interp(id, startDate)
        date = startDate + timedelta(hours=1)
        while date <= endDate:
            hourly_position = hourly_position.append(
                self.linear_interp(id, date),
                ignore_index=True
            )
            date = date + timedelta(hours=1)

        return(hourly_position)

    def get_half_hourly_position(self, id):
        """
        Interpolates Hurdat2 data for a TC from 6-hourly to half hourly data.

        The Hurdat2 database records best track data for TCs at 6 hour intervals. Uses Hurdat.linear_interp() to
        estimate half hourly track data.

        :param id: Hurdat ID for TC, e.g. AL152016 for Hurricane Nicole [2016] (15th TC in the North Atlantic in the
                   2016 season.
        :return: Data frame containing hourly interpolation of TC data.
        """
        storm = self.get_tc_by_id(id)
        startDate = np.min(storm['DATETIME'])
        endDate = np.max(storm['DATETIME'])
        hourly_position = self.linear_interp(id, startDate)
        date = startDate + timedelta(hours=.5)
        while date <= endDate:
            hourly_position = hourly_position.append(
                self.linear_interp(id, date),
                ignore_index=True
            )
            date = date + timedelta(hours=0.5)

        return(hourly_position)

    def genesis_to_lysis_filter(self, minimum_wind, keep_leading_n = 1):
        """
        Filter's storms based on needing to have wind speed > minimum_wind for
        at least one time step. Then returns all data for the storm between
        the first time it crossed the threshold (genesis) to the last time it
        falls below the threshold (lysis).

        The observation immediately preceeding the first time it crosses the 
        threshold is also returned for purposes of calculating the velocity of
        the storm at the first time point later in the analysis.

        Parameters:
            - minimum_wind (float): Wind threshold for defining genesis and 
            lysis. Units of kt. 
            - keep_leading_n (int): Must be >= 0. Keep this many observations  
            that proceed the genesis (when storm first goes >= minimum_wind).
            This option is useful for keeping the data you might need for 
            looking back in time from a rapid intensification observation. Storms 
            that never go above minimum_wind are always removed from the dataset.
        
        Returns:
            A new Hurdat instance containing the subsetted data. An additional 
            column is added to the data indicating whether or not the data point
            is a "leading" point (points in the leadup to the true genesis point).
        """

        # Get modifiable copy of dataframe of storms
        hurdat = self.storms.copy()

        # Create a column for the maximum wind seen during each storm
        hurdat['MAX_WIND'] = hurdat.groupby('ID')['WIND'].transform('max')

        # Filter to only have storms that at some point go above min. wind
        hurdat = hurdat.loc[hurdat['MAX_WIND'] >= minimum_wind]

        # Create column containing time of the TC only if it is above the min
        # wind boundary
        hurdat['DATETIME_ABOVE_MIN_WIND'] = hurdat['DATETIME'].where(hurdat['WIND'] >= minimum_wind)
        
        # Calculate the start and stop conditions for each storm based on the 
        # maximum and minimum times where the storm goes above the minimum wind
        hurdat['Start'] = hurdat.groupby('ID')['DATETIME_ABOVE_MIN_WIND'].transform('min')
        hurdat['Stop'] = hurdat.groupby('ID')['DATETIME_ABOVE_MIN_WIND'].transform('max')

        # Filter down dataset to only include storms between genesis and lysis
        hurdat = hurdat.loc[(hurdat['DATETIME'] >= hurdat['Start'] - pd.Timedelta(6*keep_leading_n, 'h')) & (hurdat['DATETIME'] <= hurdat['Stop'])]
        hurdat['LEADING'] = hurdat['DATETIME'] < hurdat['Start']
        hurdat.drop(['MAX_WIND', 'DATETIME_ABOVE_MIN_WIND', 'Start', 'Stop'], axis = 1, inplace = True)
        hurdat = hurdat.reset_index(drop = True)

        return Hurdat(data = hurdat)

    def drop_leading(self):
        """
        Function to drop the rows of the DataFrame with column "LEADING" = True.
        For use once you have used these columns to calculate the velocity that
        you need for them, and you want to discard them.
        """
        self.storms = self.storms.loc[np.logical_not(self.storms["LEADING"])]
        self.storms.drop("LEADING", axis = 1, inplace = True)
        self.storms = self.storms.reset_index(drop = True)

    def identify_events(self, threshold, col_names = None, drop_overlap = False):
        """
        Identifies each point in the dataset as either being in a rapid
        intensification event, or being in a rapid weakening event. Both are 
        defined as a change in storm wind speed of at least :threshold:kt within
        a 24 hour period.

        :param threshold: Wind speed change (kt) needed to define RI or RW event.
        :param col_names: List of 2 column names, where the first is for RI and 
            the second is for RW.
        :param drop_overlap: Boolean for whether or not observations marked as 
            both RI and RW should be removed from BOTH lists.  
        :return: Nothing. Just adds two columns to the storms dataframe, one for
                 indication of an RI event, one for indication of a RW event.
        """

        ## TODO: Add an option to remove the last RI point and first RW point
        ## (where there is commonly overlap).

        def single_storm_identify(intensities, thresh):
            TT = len(intensities)
            # If storm is to short for the algorithm, don't identify anything
            if TT < 5:
                return np.full(TT, False)

            ZZ = np.full((TT, TT), False)

            del_1 = np.diff(intensities, n = 1) # lag-1 differences
            del_4_mat = np.zeros((5, TT-4)) 
            for ii in range(1, 5):
                del_4_mat[ii,:] = [intensities[jj+ii] - intensities[jj] for jj in range(TT-4)]
            del_4 = np.amax(del_4_mat, axis = 0) # max of lag 0, 1, 2, 3, 4 differences

            AA = np.where(del_4 >= thresh)[0] # All times with RI
            BB = np.where(del_1 > 0)[0] # All increasing lag-1 times

            # Follow the algorithm outlined elsewhere for calculating RI
            for tt in range(TT - 4):
                if np.isin(tt, AA):
                    ZZ[tt, tt:(tt + 5)] = True
                    for hh in range(4, 0, -1):
                        if not np.isin(tt+hh-1, BB):
                            ZZ[tt, tt+hh] = False
                        else:
                            break
                    for hh in range(4):
                        if not np.isin(tt+hh, BB):
                            ZZ[tt, tt+hh] = False
                        else:
                            break
            YY = np.any(ZZ, 0)
            return YY

        # Get list of storm IDs
        ids = pd.unique(self.storms['ID'])

        # Loop procedure over all storms for RI and RW, keeping track of things
        # in order
        RI = []
        RW = []
        for id in ids:
            intensity = np.array(self.get_tc_by_id(id)['WIND'])
            RI_temp = single_storm_identify(intensity, threshold)
            RI.extend(RI_temp)

            RW_temp = np.flip(single_storm_identify(np.flip(intensity), threshold))
            RW.extend(RW_temp)
        
        # potentially drop the overlap points, where a point is marked as both
        # RI and RW; typically happens at the turning point if RI is immediately
        # followed by RW.
        if drop_overlap:
            RI = np.array(RI)
            RW = np.array(RW)
            where_overlap = (RI == RW) & (RI == True)
            RI[where_overlap] = False
            RW[where_overlap] = False

        # Add new columns to dataframe
        if col_names is None:
            self.storms['RI'] = RI
            self.storms['RW'] = RW
        else:
            self.storms[col_names[0]] = RI
            self.storms[col_names[1]] = RW

    
    def distance_to_land_label(self, min_distance):
        """
        Creates a column in the storms dataframe as a T/F of whether or not the
        storm's center is wihtin min_distance of land. Used to identify which
        events we want to classify as RI/RW events for study. Column is 1 if 
        within min_distance of land, 0 otherwise. 
        """

        assert 'DISTANCE' in self.storms.columns, "Need a DISTANCE column in order to compute this function."

        within_distance = self.storms['DISTANCE'] < min_distance
        self.storms['NEAR_LAND'] = within_distance


        




