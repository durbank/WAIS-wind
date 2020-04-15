# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:05:26 2020

Downscale lower resolution ERA5 data to 8 m resolution using REMA DEM

@author: Eric Johnson
"""

# Load trace lat/lons

# Import all required modules and variables
from myModule import *
import Load_ERA5 # *******Move to myModule**********
import utm
import numpy as np
from geopy.distance import distance
from numpy.matlib import repmat

# Import PAIPR-generated data
PAIPR_dir = ROOT_DIR.joinpath('data/gamma_20111109')
data_0 = import_PAIPR(PAIPR_dir)

accum_long_df = format_PAIPR(data_0, start_yr=1979, end_yr=2009)

accum_long_df = accum_long_df.drop_duplicates(subset=['trace_ID'], keep='first')

accum_easting  = []
accum_northing = []
for x in range(len(accum_long_df)):
    easting, northing, zone_number, zone_letter = utm.from_latlon(accum_long_df['Lat'].iloc[x], accum_long_df['Lon'].iloc[x])
    accum_easting.append(easting)
    accum_northing.append(northing)

accum_long_df['Easting']  = accum_easting
accum_long_df['Northing'] = accum_northing

wind_u_10m, wind_v_10m, precip, vLat, vLon, vTime, elevations = Load_ERA5.Load_ERA5()

def Downscale_ERA5(wind_u_10m, wind_v_10m, precip, elevations, vLat, vLon, vTime, accum_long_df):
    
    # Alternate idea for supplying fewer variables as inputs
    # class ERA5:
    #     pass

    # ERA5_Data = ERA5()
    # ERA5_Data.wind_u_10m = wind_u_10m
    # ERA5_Data.wind_v_10m = wind_v_10m
    # ERA5_Data.precip = precip
    # ERA5_Data.vLat = vLat_ERA
    # ERA5_Data.vLon = vLon_ERA
    # ERA5_Data.time = vTime

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Load ERA5 data

    mElevations_ERA5 = elevations
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Load trace data

    trace_data_df = accum_long_df.set_index('trace_ID')
    trace_data_df['wind_u'] = ""
    trace_data_df['wind_v'] = ""
    trace_data_df['precip'] = ""

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Convert lat/lon to gridded UTM
    
    mLat_ERA = np.transpose(repmat(vLat,len(vLon),1))
    mLon_ERA = repmat(vLon,len(vLat),1)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Downscale ERA at each trace location

    for x in range(len(trace_data_df)):
        
        try:
            # Find closest ERA5 location to each trace ID
            closest_idx = np.sqrt(
                np.square(np.matrix(mLat_ERA - trace_data_df['Lat'][x])) + 
                np.square((np.matrix(mLon_ERA - trace_data_df['Lon'][x])))
                             ).argmin()
            closest_idx = np.unravel_index(closest_idx, np.shape(mLat_ERA))
            
            # Calculate distances to nearest ERA5 grid points
            m3x3_ERA_Lat = vLat[closest_idx[0]-1:closest_idx[0]+2]
            m3x3_ERA_Lat = np.transpose(repmat(m3x3_ERA_Lat,3,1))
            m3x3_ERA_Lon  = vLon[closest_idx[1]-1:closest_idx[1]+2]
            m3x3_ERA_Lon = repmat(m3x3_ERA_Lon,3,1)
            m3x3_Distances = np.zeros([3,3])
            for m in range(3):
                for n in range(3):
                    era_lat = m3x3_ERA_Lat[m,n]
                    era_lon = m3x3_ERA_Lon[m,n]
                    trace_lat = trace_data_df['Lat'][x]
                    trace_lon = trace_data_df['Lon'][x]
                    m3x3_Distances[m,n] = distance([era_lat,era_lon], [trace_lat, trace_lon]).km
            m3x3_Weights = m3x3_Distances / np.sum(m3x3_Distances[:])
            # Calculate weighted mean (based on distance) of surrounding ERA5 values
            # for each time step
            vWind_u = []
            vWind_v = []
            vPrecip = []
            for y in range(np.shape(wind_u_10m)[0]):
                temp = m3x3_Weights * wind_u_10m[y,closest_idx[0]-1:closest_idx[0]+2,closest_idx[1]-1:closest_idx[1]+2]
                vWind_u.append(np.sum(temp[:]))
                temp = m3x3_Weights * wind_v_10m[y,closest_idx[0]-1:closest_idx[0]+2,closest_idx[1]-1:closest_idx[1]+2]
                vWind_v.append(np.sum(temp[:]))
                temp = m3x3_Weights * precip[y,closest_idx[0]-1:closest_idx[0]+2,closest_idx[1]-1:closest_idx[1]+2]
                vPrecip.append(np.sum(temp[:]))
                
            trace_data_df['wind_u'][x] = vWind_u
            trace_data_df['wind_v'][x] = vWind_v
            trace_data_df['precip'][x] = vPrecip
            
        except:
            pass

        print([x/39024*100,'%'])









































































    return
