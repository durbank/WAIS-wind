# Script to define various custom functions used in WAIS-wind project

# Import modules required for various function defintions
import numpy as np
import pandas as pd
import geopandas as gpd
import richdem as rd
from pathlib import Path
import rasterio as rio
import requests
import shutil
import os
import cdsapi
import netCDF4
import sys
import utm
from geopy.distance import distance
from numpy.matlib import repmat

# Set project root directory
ROOT_DIR = Path(__file__).parent.parent


# Function to import and concatenate PAIPR .csv files
def import_PAIPR(input_dir):
    """
    Function to import PAIPR-derived accumulation data.
    This concatenates all files within the given directory into a single pandas dataframe.

    Parameters:
    input_dir (pathlib.PosixPath): Absolute path of directory containing .csv files of PAIPR results (gamma-distributions fitted to accumulation curves). This directory can contain multiple files, which will be concatenated to a single dataframe.
    """

    data = pd.DataFrame()
    for file in input_dir.glob("*.csv"):
        data_f = pd.read_csv(file)
        data = data.append(data_f)
    return data


# Function to format imported PAIPR data
def format_PAIPR(data_raw, start_yr, end_yr):
    """

    """
    # Remove time series with data missing from period
    # of interest (and clip to period of interest)
    traces = data_raw.groupby(['Lat', 'Lon', 'elev'])
    data = data_raw.assign(trace_ID = traces.ngroup())
    traces = data.groupby('trace_ID')
    data = traces.filter(
        lambda x: min(x['Year']) <= start_yr 
        and max(x['Year']) >= end_yr)
    data = data.query(f"Year >= {start_yr} & Year < {end_yr}")

    # Ensure each trace has only one time series 
    # (if not, take the mean of all time series)
    data = data.groupby(['trace_ID', 'Year']).mean()

    # Generate descriptive statistics based on imported 
    # gamma-fitted parameters
    alpha = data['gamma_shape']
    alpha.loc[alpha<1] = 1
    beta = 1/data['gamma_scale']
    mode_accum = (alpha-1)/beta
    var_accum = alpha/beta**2

    # New df (in long format) with accum data assigned
    data_long = (
        data.filter(['trace_ID', 'Lat', 'Lon', 'elev', 'Year'])
        .assign(accum = mode_accum, std = np.sqrt(var_accum))
        .reset_index()
    )
    return data_long


# Function to download REMA data, unzip, and save to disk
def get_REMA(tile_idx, output_dir):
    """

    """

    for idx, row in tile_idx.iterrows():
        f_dir = output_dir.joinpath(row.tile)

        if not f_dir.exists():
            f_dir.mkdir(parents=True)
            zip_path = f_dir.joinpath('tmp.tar.gz')
            r = requests.get(row.fileurl, stream=True)
            print(f"Downloading tile {f_dir.name}")
            with open(zip_path, 'wb') as zFile:
                for chunk in r.iter_content(
                        chunk_size=1024*1024):
                    if chunk:
                        zFile.write(chunk)
            print(f"Unzipping tile {f_dir.name}")
            shutil.unpack_archive(zip_path, f_dir)
            os.remove(zip_path)
        else:
            print(f"REMA tile {f_dir.name} already exists locally, moving to next download")
    print("All requested files downloaded")

# Function to calculate surface model variable rasters
def calc_topo(dem_path):
    """
    Calculates slope and aspect from given DEM and saves output.
    The function checks to see whether a slope/aspect file has already been created so as to avoid needless processing.
    
    Parameters:
    dem_path (pathlib.PosixPath): The relative or absolute path to an input DEM file.

    Dependencies: 
    richdem module
    GDAL binaries
    pathlib module
    """
    slope_path = Path(
        str(dem_path).replace("dem", "slope"))
    aspect_path = Path(
        str(dem_path).replace("dem", "aspect"))

    if ((not slope_path.is_file()) or 
            (not aspect_path.is_file())):
        
        dem = rd.LoadGDAL(str(dem_path))

    if not slope_path.is_file():
        slope = rd.TerrainAttribute(
            dem, attrib='slope_riserun')
        rd.SaveGDAL(str(slope_path), slope)
    
    if not aspect_path.is_file():
        aspect = rd.TerrainAttribute(dem, attrib='aspect')
        rd.SaveGDAL(str(aspect_path), aspect)
        
    
# Function to extract raster values at given points
def topo_vals(tile_dir, trace_locs):
    """
    Extracts elevation, slope, and aspect values at given locations.
    
    Parameters:
    tile_dir (pathlib.PosixPath): The relative or absolute path to a directory containing REMA tile DSM, slope and aspect geotiffs.
    trace_locs (geopandas.geodataframe.GeoDataFrame): A geodataframe containing the locations at which to extract raster data. These data should have the same coordinate reference system as the raster data, with the geometries stored in a column named "geometry".

    Dependencies: Requires the rasterio (as rio) module and, by extension, GDAL binaries.
    Requires the geopandas module.
    """
    trace_locs = (
        trace_locs.assign(elev=None)
        .assign(slope=None).assign(aspect=None))

    coords = (
        [(x,y) for x, y in zip(
            trace_locs.geometry.x, trace_locs.geometry.y)]
    )
    
    # Extract elevation values for all points within tile
    tile_path = [
        file for file in tile_dir.glob("*dem.tif")][0]
    src = rio.open(tile_path)
    tile_vals = np.asarray(
        [x[0] for x in src.sample(coords, masked=True)])
    tile_mask = ~np.isnan(tile_vals)
    trace_locs.elev[tile_mask] = tile_vals[tile_mask]
    src.close()

    # Extract slope values for all points within tile
    tile_path = [
        file for file in tile_dir.glob("*slope.tif")][0]
    src = rio.open(tile_path)
    tile_vals = np.asarray(
        [x[0] for x in src.sample(coords, masked=True)])
    tile_mask = ~np.isnan(tile_vals)
    trace_locs.slope[tile_mask] = tile_vals[tile_mask]
    src.close()

    # Extract aspect values for all points within tile
    tile_path = [
        file for file in tile_dir.glob("*aspect.tif")][0]
    src = rio.open(tile_path)
    tile_vals = np.asarray(
        [x[0] for x in src.sample(coords, masked=True)])
    tile_mask = ~np.isnan(tile_vals)
    trace_locs.aspect[tile_mask] = tile_vals[tile_mask]
    src.close()

    return trace_locs



def Download_ERA5():
    # -*- coding: utf-8 -*-

    """
    Created on Mon Mar  2 14:30:08 2020

    @author: Eric Johnson
    """
    # xarray for netcdfs
    # ERA5 Login Credentials:
    # username: u0929154@utah.edu
    # password: DataScience2020
    
    # Extent: Longitude (-115.86289, -94.9989) Latitude (-79.9461, -75.58234) 
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This section modified from https://stackoverflow.com/questions/44210656/how-to-check-if-a-module-is-installed-in-python-and-if-not-install-it-within-t/44210735
    # Check to insure all necessary packages are installed, install missing packages
    # import sys
    # import subprocess
    # import pkg_resources
    
    # required = {'cdsapi', 'netCDF4'}
    # installed = {pkg.key for pkg in pkg_resources.working_set}
    # missing = required - installed
    
    # if missing:
    #     python = sys.executable
    #     subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
    ###############################################################################

    ## Download ERA5 variables into files
    
    # url = 'https://cds.climate.copernicus.eu/api/v2'
    # key = '38293:414bf99c-171e-48e6-b13d-75661652a8de'
    # login_file = '.cdsapirc'
    
    
    # 10-m u-component of wind
    c = cdsapi.Client()
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : '10m_u_component_of_wind',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('10m_u_component_of_wind.nc')
    
    
    # 10-m v component of wind
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : '10m_v_component_of_wind',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('10m_v_component_of_wind.nc')
    
    
    # Total precipitation
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : 'total_precipitation',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('total_precipitation.nc')

    return


def Extract_ERA5():

    # -*- coding: utf-8 -*-

    """
    Created on Mon Mar  2 14:30:08 2020

    @author: Eric Johnson
    """
    # xarray for netcdfs

    # ERA5 Login Credentials:
    # username: u0929154@utah.edu
    # password: DataScience2020
    
    # Extent: Longitude (-115.86289, -94.9989) Latitude (-79.9461, -75.58234) 
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This section modified from https://stackoverflow.com/questions/44210656/how-to-check-if-a-module-is-installed-in-python-and-if-not-install-it-within-t/44210735
    # Check to insure all necessary packages are installed, install missing packages
    # import sys
    # import subprocess
    # import pkg_resources
    
    # required = {'cdsapi', 'netCDF4'}
    # installed = {pkg.key for pkg in pkg_resources.working_set}
    # missing = required - installed
    
    # if missing:
    #     python = sys.executable
    #     subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
    ###############################################################################
    
    ## Initial conditions
    
    # Extents
    iLat_min = -75.58234
    iLat_max = -79.9461
    iLon_min = -115.86289
    iLon_max = -75.58234
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Download ERA5 variables into files
    
    # url = 'https://cds.climate.copernicus.eu/api/v2'
    # key = '38293:414bf99c-171e-48e6-b13d-75661652a8de'
    # login_file = '.cdsapirc'
    
    
    # 10-m u-component of wind
    c = cdsapi.Client()
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : '10m_u_component_of_wind',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('10m_u_component_of_wind.nc')
    
    
    # 10-m v component of wind
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : '10m_v_component_of_wind',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('10m_v_component_of_wind.nc')
    
    
    # Total precipitation
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : 'total_precipitation',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('total_precipitation.nc')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extract necessary portions of files
    
    # Extract lat/lon data
    data = netCDF4.Dataset('10m_u_component_of_wind.nc')
    vLat = data.variables['latitude']
    vLat = np.asarray(vLat[:])
    vLon = data.variables['longitude']
    vLon = np.asarray(vLon[:])
    vLon[vLon > 180] = -(360 - vLon[vLon > 180]) # Is this right? What's with the 360 deg latitudes?
    
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
    precip = data.variables['prcp']
    
    # Extract required portions of each variable
    wind_u_10m = wind_u_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    wind_v_10m = wind_v_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    precip = precip[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]

    return wind_u_10m, wind_v_10m, precip


def Load_ERA5():
    
    # -*- coding: utf-8 -*-
    """
    Created on Thu Mar 26 11:50:48 2020

    @author: Eric Johnson
    """
    
    
    # Add path to pre-downloaded ERA5 data
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


def Downscale_ERA5(wind_u_10m, wind_v_10m, precip, elevations, vLat, vLon, vTime, accum_long_df):
    
    # -*- coding: utf-8 -*-
    """
    Created on Mon Mar 23 12:05:26 2020

    Downscale lower resolution ERA5 data to 8 m resolution using REMA DEM

    @author: Eric Johnson
    """

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