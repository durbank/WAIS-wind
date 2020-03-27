# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:50:48 2020

@author: Eric Johnson
"""

def Load_ERA5():
    
    import netCDF4
    import numpy as np
    
    # Extents
    iLat_min = -75.58234
    iLat_max = -79.9461
    iLon_min = -115.86289
    iLon_max = -75.58234
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extract necessary portions of files
    
    # Extract lat/lon data
    data = netCDF4.Dataset('10m_u_component_of_wind.nc')
    vLat = data.variables['latitude']
    vLat = np.asarray(vLat[:])
    vLon = data.variables['longitude']
    vLon = np.asarray(vLon[:])
    vLon[vLon > 180] = -(360 - vLon[vLon > 180]) # Is this right? What's with the 360 deg latitudes?
    vTime = data.variables['time']
    vTime = vTime[:]
    
    # Find indices for required extents
    iLat_min_idx = np.argmin(np.abs(vLat-iLat_min))
    iLat_max_idx = np.argmin(np.abs(vLat-iLat_max))
    iLon_min_idx = np.argmin(np.abs(vLon-iLon_min))
    iLon_max_idx = np.argmin(np.abs(vLon-iLon_max))
    
    # Read NetCDF files
    data = netCDF4.Dataset('10m_u_component_of_wind.nc')
    wind_u_10m = data.variables['u10']
    data = netCDF4.Dataset('10m_v_component_of_wind.nc')
    wind_v_10m = data.variables['v10']
    data = netCDF4.Dataset('total_precipitation.nc')
    precip = data.variables['tp']
    
    # Extract required portions of each variable
    wind_u_10m = wind_u_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    wind_v_10m = wind_v_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    precip = precip[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx] 
    vLat = vLat[iLat_min_idx:iLat_max_idx]
    vLon = vLon[iLon_min_idx:iLon_max_idx]
    
    
    return wind_u_10m, wind_v_10m, precip, vLat, vLon, vTime