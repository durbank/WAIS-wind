# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:50:48 2020

@author: Eric Johnson
"""

def Load_ERA5():
    
    import netCDF4
    import numpy as np
    from pathlib import Path
    import sys
    
    # Add path to pre-downloaded ERA5 data
    ROOT_DIR = Path(__file__).parent.parent
    ERA5_dir = ROOT_DIR.joinpath('data/ERA5')
    sys.path.append(ERA5_dir)
    
    # Extents
    iLat_min = -75.58234
    iLat_max = -79.9461
    iLon_min = -115.86289
    iLon_max = -75.58234
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extract necessary portions of files
    
    # Extract lat/lon data
    
    data = netCDF4.Dataset(ERA5_dir.joinpath('10m_u_component_of_wind.nc'))
    vLat = data.variables['latitude']
    vLat = np.asarray(vLat[:])
    vLon = data.variables['longitude']
    vLon = np.asarray(vLon[:])
    vLon[vLon > 180] = -(360 - vLon[vLon > 180]) 
    vTime = data.variables['time']
    vTime = vTime[:]
    
    # Find indices for required extents
    iLat_min_idx = np.argmin(np.abs(vLat-iLat_min))-1 # +/-1 to add buffer around ERA5 data for downscaling
    iLat_max_idx = np.argmin(np.abs(vLat-iLat_max))+1
    iLon_min_idx = np.argmin(np.abs(vLon-iLon_min))-1
    iLon_max_idx = np.argmin(np.abs(vLon-iLon_max))+1
    
    # Read NetCDF files
    data = netCDF4.Dataset(ERA5_dir.joinpath('10m_u_component_of_wind.nc'))
    wind_u_10m = data.variables['u10']
    data = netCDF4.Dataset(ERA5_dir.joinpath('10m_v_component_of_wind.nc'))
    wind_v_10m = data.variables['v10']
    data = netCDF4.Dataset(ERA5_dir.joinpath('total_precipitation.nc'))
    precip = data.variables['tp']
    data = netCDF4.Dataset(ERA5_dir.joinpath('ERA5_Orography.nc'))
    elevations = data.variables['z']
    elevations = elevations[:,:,:]
    elevations = np.kron(np.squeeze(elevations),np.ones([5,5]))

    
    # Extract required portions of each variable
    wind_u_10m = wind_u_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    wind_v_10m = wind_v_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    precip = precip[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx] 
    elevations = elevations[iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx] / 9.80665
    vLat = vLat[iLat_min_idx:iLat_max_idx]
    vLon = vLon[iLon_min_idx:iLon_max_idx]
    
    
    return wind_u_10m, wind_v_10m, precip, vLat, vLon, vTime, elevations