# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:05:26 2020

Downscale lower resolution ERA5 data to 8 m resolution using REMA DEM

@author: Eric Johnson
"""

def Downscale_ERA5(wind_u_10m, wind_v_10m, precip, vLat_ERA, vLon_ERA, vTime):
    
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

    import utm
    import numpy as np

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Load ERA5 data
    
    # wind_u, wind_v, precip
    # iLat_min, iLat_max, iLong_min, iLong_max
    # Elevation_ERA5

    mElevations_ERA5 = [[1,2,3], [4,5,6], [7,8,9]] # placeholder
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Load REMA DEM
    # mDEM
    # iLat_min, iLat_max, iLong_min, iLong_max
    # iResolution

    mLat_REMA = 0 # placeholder
    mLon_REMA = 0 # placeholder

    vLat_REMA = mREMA_Lat[:,1] # placeholder
    vLon_REMA = mREMA_Lon[1,:] # placeholder
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Mask out areas with radar trace data
    
    mMask_ERA  = 0 # placeholder
    mMask_REMA = 0 # placeholder

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Convert lat/lon to gridded UTM
    
    # ERA5
    mEasting_ERA = np.zeros([len(vLat_ERA),len(vLon_ERA)])
    mNorthing_ERA = np.zeros([len(vLat_ERA),len(vLon_ERA)])
        
    for lat_idx in range(len(vLat_ERA)):
        lat = vLat_ERA(lat_idx)
        for lon_idx in range(len(vLon_ERA)):
            lon = vLon_ERA(lon_idx)
            easting, northing, zone_number, zone_letter = utm.from_latlon(lat, lon)
            mEasting_ERA[lat_idx,lon_idx]  = easting
            mNorthing_ERA[lat_idx,lon_idx] = northing
            
    # REMA
    mEasting_REMA = np.zeros([len(vLat_REMA),len(vLon_REMA)])
    mNorthing_REMA = np.zeros([len(vLat_REMA),len(vLon_REMA)])
        
    for lat_idx in range(len(vLat_REMA)):
        lat = vLat_REMA(lat_idx)
        for lon_idx in range(len(vLon_REMA)):
            lon = vLon_REMA(lon_idx)
            easting, northing, zone_number, zone_letter = utm.from_latlon(lat, lon)
            mEasting_REMA[lat_idx,lon_idx]  = easting
            mNorthing_REMA[lat_idx,lon_idx] = northing

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Downscale ERA at each point onto high resolution (DEM) grid
            
    mWindU_hd  = np.zeros(np.shape(mEasting_REMA))
    mWindV_hd  = np.zeros(np.shape(mEasting_REMA))
    mPrecip_hd = np.zeros(np.shape(mEasting_REMA))

    for m in range(np.shape(mEasting_REMA)[0]):
        for n in range(np.shape(mEasting_REMA)[1]):
            try: # try/except necessary to avoid edges of ERA5 data matrix
                
                # Find indices of nearest ERA grid point
                iERA_e_idx = (np.abs(mEasting_ERA  -  mEasting_REMA(m,n))).argmin() # does this find index or value? want index
                iERA_n_idx = (np.abs(mNorthing_ERA - mNorthing_REMA(m,n))).argmin() # does this find index or value? want index
                # Calculate distances to nearest ERA5 grid points
                m3x3_ERA_Easting  =  mEasting_ERA[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
                m3x3_ERA_Northing = mNorthing_ERA[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
                m3x3_Distances = np.sqrt( (m3x3_ERA_Easting - mEasting_REMA[m,n])**2 + (m3x3_ERA_Northing - mNorthing_REMA[m,n])**2 )
                m3x3_Weights = m3x3_Distances / sum(m3x3_Distances[:])
                # Calculate weighted mean (based on distance) of surrounding ERA5 values
                temp = m3x3_Weights * wind_u_10m[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
                mWindU_hd[m,n] = sum(temp[:])
                temp = m3x3_Weights * wind_v_10m[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
                mWindV_hd[m,n] = sum(temp[:])
                temp = m3x3_Weights * precip[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
                mPrecip_hd[m,n] = sum(temp[:])
                
            except:
                pass

















































































    return
