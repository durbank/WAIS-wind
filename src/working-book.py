# %%[markdown]
# # Working Notebook
# 
# This is a notebook to keep an up-to-date working version of the project.

# ## Importing/defining project modules, functions, and variables
#
# A number of the functions and variables defined in our project are currently housed in a secondary script called `myModule.py`.
# We include this list of functions below for convienence, but usually import these definitions using 
# ```{python}
# from myModule import *
# ````
# and call them later on in the notebook.

# %%
# # Import required modules
# from myModule import *

# Required imports
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import richdem as rd
from pathlib import Path
import rasterio as rio
import requests
import shutil
import os

# Set project root directory
ROOT_DIR = Path(os.getcwd()).parent
# ROOT_DIR = Path(__file__).parent.parent

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

# %% [markdown]
# ## Data loading, cleaning, and formatting
#
# The first step is to import our annual accumulation data.
# These data are stored in .csv files that include the Lat/Lon position of each accumulation time series, the year, and parameters representing a gamma-distribution fitted to distribution of accumulation estimates for each year.
# We import and concatenate these files using the `import_PAIPR` function.

# %%
# Import PAIPR-generated data
PAIPR_dir = ROOT_DIR.joinpath('data/gamma_20111109')
data_0 = import_PAIPR(PAIPR_dir)

# %%[markdown]

# We next perform some data cleaning and formatting on the imported accumulation data, using `format_PAIPR`.
# This function groups data into individual time series based on location and clips the results to only include data within a given time frame (in this case 1979-2009).
# It also calculates mean accumulation and accumulation standard deviation for each year, based on the fitted gamma distribution parameters.
# The end result is a `pandas` dataframe containing 30-year time series of accumulation estimates (and uncertainty) for ~40,000 different locations in West Antarctica.
#  
# %%
# 
accum_long = format_PAIPR(data_0, start_yr=1979, end_yr=2009)
traces = accum_long.groupby('trace_ID')

# %%[markdown]
#
# We also need to accurately represent the spatial locations of our data.
# This will be particularly important when we wish to combine it with other spatial datasets to ensure proper alignment of results.
# We use the `geopandas` module to accomplish this, converting our data to various GeoDataFrames.

# %%
# New accum and std dfs in wide format
accum = accum_long.pivot(
    index='trace_ID', columns='Year', values='accum')
accum_std = accum_long.pivot(
    index='trace_ID', columns='Year', values='std')

# Create df for mean annual accumulation
accum_mu = traces.aggregate(np.mean).drop('Year', axis=1)
accum_mu = gpd.GeoDataFrame(
    accum_mu, geometry=gpd.points_from_xy(
        accum_mu.Lon, accum_mu.Lat)
).drop(['Lat', 'Lon', 'elev'], axis=1)
accum_mu.crs = "EPSG:4326"

# Create a gdf from accum df
accum_gdf = gpd.GeoDataFrame(
    accum, geometry=gpd.points_from_xy(
        traces.aggregate(np.mean)['Lon'], 
        traces.aggregate(np.mean)['Lat'])
)
accum_gdf.crs = "EPSG:4326"

# %%[markdown]

# We also import a shapefile containing the boundaries of Antarctica, for use in plotting.
# We also convert our GeoDataFrames into the same coordinate reference system as our Antarctic map data, a [Polar Stereographic projection](https://nsidc.org/data/polar-stereo/ps_grids.html).
# This is particularly important for mapping and analysis in Antarctica, as traditional meanings of North/South/East/West break down near the poles.

# %%
# Import Antarctic outline shapefile
ant_path = ROOT_DIR.joinpath(
    'data/Ant_basemap/Coastline_medium_res_polygon.shp')
ant_outline = gpd.read_file(ant_path)

# Convert accum crs to same as Antarctic outline
accum_gdf = accum_gdf.to_crs(ant_outline.crs)
accum_mu = accum_mu.to_crs(ant_outline.crs)

# %%[markdown]

# As part of our initial data exploration, we plot an inset map marking the location of our data relative to Antarctica.
#  We further plot the mean accumulation across the data range (in millimeters water equivalent).

# %%
# Plot inset map
base = ant_outline.plot(color='grey', edgecolor='black')
accum_gdf.sample(n=1000).plot(ax=base, color='red')
plt.xlabel('Easting')
plt.ylabel('Northing')

# Plot mean accumulation spatial map
accum_mu.sample(n=2500).plot(
    column='accum', legend=True, legend_kwds={
        'label': 'Mean annual accumulation (mm/yr)'})
plt.xlabel('Easting')
plt.ylabel('Northing')

#%%[markdown]
# ## Extraction of surface topography features
#
# In order to investigate the relationships between annual accumulation and surface topography, we need topography data for our study region.
# We utilize the Reference Elevation Model of Antarctica (REMA), which consists of an 8-m resolution digital surface model for Antarctica. 
# This uses the REMA tile index shapefile (found at [this website](http://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/)) to determine which REMA tiles cover our dataset.
# These tiles are fairly large (the full dataset for our study region is ~40 Gb), so the `get_REMA` function checks for whether any of the required tiles were previously downloaded, and only downloads missing tiles.
# To save time and to avoid unnecessarily large file uploads for this project, we leave this code commented out and simply load a dataset containing the extracted data of interest at a later point.

# %%

# # Define path to index shapefile
# index_file = ROOT_DIR.joinpath("data/REMA", 
#     "REMA_Tile_Index_Rel1.1", "REMA_Tile_Index_Rel1.1.shp")

# # Import shapefile of DEM tile locations
# dem_index = gpd.read_file(index_file)

# # Keep only DEMs that contain accum traces
# dem_index = (
#     gpd.sjoin(dem_index, accum_mu, op='contains')
#     .iloc[:,0:dem_index.shape[1]]).drop_duplicates()

# # Find and download missing REMA DSM tiles
# REMA_outdir = ROOT_DIR.joinpath("data/REMA/tiles_8m_v1.1")
# tiles_list = pd.DataFrame(
#     dem_index.drop(columns='geometry'))
# # For development purposes, only using first 2 tiles
# # get_REMA(tiles_list, REMA_outdir) 
# get_REMA(tiles_list[:2], REMA_outdir) 

# %%[markdown]

# We next calculate the slope and aspect of the surface topography using the downloaded REMA data and our `calc_topo` function.
# This function similarly checks for whether slope/aspect rasters were previously calculated and, in order to avoid needless processing, only performs the calculations on missing raster data.
# These lines are once again commented out, as running them would require uploaded 10's of GBs of data.

# %%

# # Generate list of paths to downloaded DEMs for topo
# #  calculations
# REMA_dir = ROOT_DIR.joinpath("data/REMA/tiles_8m_v1.1")
# dem_list = [path for path in REMA_dir.glob('**/*dem.tif')]

# # Calculate slope and aspect for each DEM
# # For development purposes, only processing first in list
# # [calc_topo(dem) for dem in dem_list]
# [calc_topo(dem) for dem in dem_list[:1]]

# %%[markdown]

# We next extract the calculated elevation, slope, and aspect values at each accumulation time series location (using our `topo_vals` fucntion) and save the results to disk (as a GeoJSON file) for later import and merging with the dataframe containing the accumulation data.
# Similar to the above step, to avoid large file sizes we performed these calculations previously, and instead load the final product.

# %%
# # Extract elevation, slope, and aspect values for each trace 
# # location
# raster_vals = accum_mu.drop(['accum', 'std'], axis=1)
# tif_dirs = [path.parent for path in dem_list[:1]]
# for path in tif_dirs:
#     raster_vals = topo_vals(path, raster_vals)

# # Save raster_vals to GeoJSON
# file_out = ROOT_DIR.joinpath('data/REMA/raster_vals.geojson')
# raster_vals.to_file(file_out, driver='GeoJSON')

# %%[markdown]

# We now load the pre-processed topography results generated using the REMA raster data and merge them into the accumulation geodataframes.
# 
 # *NOTE:* We currently have only run the above topo functions for a small subset of REMA tiles. The current version of 'raster_vals' is therefore missing most of the data, but the method is ready to scale to all REMA tiles.

# %%

# Import topography values from rasters
raster_vals = gpd.read_file(
    ROOT_DIR.joinpath('data/REMA/raster_vals.geojson'))

# Merge topo values into accumulation geodataframe
accum_mu = accum_mu.merge(
    raster_vals.drop('geometry', axis=1), 
    how='left', on='trace_ID')

# %%[markdown]

# We then download the ERA5 reanalysis gridded climate data using the online API. 
# To save time and avoid unnecessarily large files (and also because it requires login credentials to download), we comment out the code and simply load the data at a later step.

# %%

# def Download_ERA5():
    
#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     ## Imports
    
#     import cdsapi
#     import netCDF4
#     import numpy as np
    
#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     ## Download ERA5 variables into files
    
#     # 10-m u-component of wind
#     c = cdsapi.Client()
#     r = c.retrieve(
#         'reanalysis-era5-land-monthly-means', {
#                 'format'      : 'netcdf',
#                 'product_type': 'monthly_averaged_reanalysis',
#                 'variable'    : '10m_u_component_of_wind',
#                 'stream'      : 'moda',
#                 'year'        : [
#                     '1980','1981','1982','1983','1984','1985','1986','1987','1988',
#                     '1989','1990','1991','1992','1993','1994','1995','1996','1997',
#                     '1998','1999','2000','2001','2002','2003','2004','2005','2006',
#                     '2007','2008','2009','2010','2011','2012','2013','2014','2015'
#                     ],
#                 'month'       : [
#                     '01','02','03','04','05','06','07','08','09','10','11','12'
#                     ],
#                 'time'        : '00:00'
#         })
#     r.download('10m_u_component_of_wind.nc')
    
    
#     # 10-m v component of wind
#     r = c.retrieve(
#         'reanalysis-era5-land-monthly-means', {
#                 'format'      : 'netcdf',
#                 'product_type': 'monthly_averaged_reanalysis',
#                 'variable'    : '10m_v_component_of_wind',
#                 'stream'      : 'moda',
#                 'year'        : [
#                     '1980','1981','1982','1983','1984','1985','1986','1987','1988',
#                     '1989','1990','1991','1992','1993','1994','1995','1996','1997',
#                     '1998','1999','2000','2001','2002','2003','2004','2005','2006',
#                     '2007','2008','2009','2010','2011','2012','2013','2014','2015'
#                     ],
#                 'month'       : [
#                     '01','02','03','04','05','06','07','08','09','10','11','12'
#                     ],
#                 'time'        : '00:00'
#         })
#     r.download('10m_v_component_of_wind.nc')
    
    
#     # Total precipitation
#     r = c.retrieve(
#         'reanalysis-era5-land-monthly-means', {
#                 'format'      : 'netcdf',
#                 'product_type': 'monthly_averaged_reanalysis',
#                 'variable'    : 'total_precipitation',
#                 'stream'      : 'moda',
#                 'year'        : [
#                     '1980','1981','1982','1983','1984','1985','1986','1987','1988',
#                     '1989','1990','1991','1992','1993','1994','1995','1996','1997',
#                     '1998','1999','2000','2001','2002','2003','2004','2005','2006',
#                     '2007','2008','2009','2010','2011','2012','2013','2014','2015'
#                     ],
#                 'month'       : [
#                     '01','02','03','04','05','06','07','08','09','10','11','12'
#                     ],
#                 'time'        : '00:00'
#         })
#     r.download('total_precipitation.nc')
    

#     return

# %%[markdown]

# We then load the ERA5 data and clip it down to only the required extents for this project.

# %%
def Load_ERA5():
    
    import netCDF4
    
    # Extents
    iLat_min = -75.58234
    iLat_max = -79.9461
    iLon_min = -115.86289
    iLon_max = -75.58234
    
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
    
    # ERA5 Elevations
    data = netCDF4.Dataset('ERA5_Orography.nc')
    elv_Lat = data.variables['latitude']
    elv_Lat = elv_Lat[:]
    elv_Lon = data.variables['longitude']
    elv_Lon = elv_Lon[:]
    elv_orog = data.variables['z'] # Surface geopotential
    elv_orog = elv_orog[:,:,:]
    g = 9.80665 # Gravity
    mElev = elv_orog/g # Surface geopotential height = surface geopotential / gravity
    
    return wind_u_10m, wind_v_10m, precip, vLat, vLon, vTime, mElev

# %%[markdown]

# Next, the ERA5 data are downscaled from their native low resolution to a higher resolution using a simple linear interpolation method. 
# This portion of the code is still under development, and thus does not currently run. 
# It is being included here only to show where the current progress of the project is at, but remains commented out.

# %%

# def Downscale_ERA5(wind_u_10m, wind_v_10m, precip, vLat_ERA, vLon_ERA, vTime):

#     import utm
#     import numpy as np

#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     # Load ERA5 data
    
#     # wind_u, wind_v, precip
#     # iLat_min, iLat_max, iLong_min, iLong_max
#     # Elevation_ERA5

#     mElevations_ERA5 = [[1,2,3], [4,5,6], [7,8,9]] # placeholder
    
    
#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     # Load REMA DEM
#     # mDEM
#     # iLat_min, iLat_max, iLong_min, iLong_max
#     # iResolution

#     mLat_REMA = 0 # placeholder
#     mLon_REMA = 0 # placeholder

#     vLat_REMA = mREMA_Lat[:,1] # placeholder
#     vLon_REMA = mREMA_Lon[1,:] # placeholder
    
#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     # Mask out areas with radar trace data
    
#     mMask_ERA  = 0 # placeholder
#     mMask_REMA = 0 # placeholder

#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     # Convert lat/lon to gridded UTM
    
#     # ERA5
#     mEasting_ERA = np.zeros([len(vLat_ERA),len(vLon_ERA)])
#     mNorthing_ERA = np.zeros([len(vLat_ERA),len(vLon_ERA)])
        
#     for lat_idx in range(len(vLat_ERA)):
#         lat = vLat_ERA(lat_idx)
#         for lon_idx in range(len(vLon_ERA)):
#             lon = vLon_ERA(lon_idx)
#             easting, northing, zone_number, zone_letter = utm.from_latlon(lat, lon)
#             mEasting_ERA[lat_idx,lon_idx]  = easting
#             mNorthing_ERA[lat_idx,lon_idx] = northing
            
#     # REMA
#     mEasting_REMA = np.zeros([len(vLat_REMA),len(vLon_REMA)])
#     mNorthing_REMA = np.zeros([len(vLat_REMA),len(vLon_REMA)])
        
#     for lat_idx in range(len(vLat_REMA)):
#         lat = vLat_REMA(lat_idx)
#         for lon_idx in range(len(vLon_REMA)):
#             lon = vLon_REMA(lon_idx)
#             easting, northing, zone_number, zone_letter = utm.from_latlon(lat, lon)
#             mEasting_REMA[lat_idx,lon_idx]  = easting
#             mNorthing_REMA[lat_idx,lon_idx] = northing

#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     # Downscale ERA at each point onto high resolution (DEM) grid
            
#     mWindU_hd  = np.zeros(np.shape(mEasting_REMA))
#     mWindV_hd  = np.zeros(np.shape(mEasting_REMA))
#     mPrecip_hd = np.zeros(np.shape(mEasting_REMA))

#     for m in range(np.shape(mEasting_REMA)[0]):
#         for n in range(np.shape(mEasting_REMA)[1]):
#             try: # try/except necessary to avoid edges of ERA5 data matrix
                
#                 # Find indices of nearest ERA grid point
#                 iERA_e_idx = (np.abs(mEasting_ERA  -  mEasting_REMA(m,n))).argmin() # does this find index or value? want index
#                 iERA_n_idx = (np.abs(mNorthing_ERA - mNorthing_REMA(m,n))).argmin() # does this find index or value? want index
#                 # Calculate distances to nearest ERA5 grid points
#                 m3x3_ERA_Easting  =  mEasting_ERA[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
#                 m3x3_ERA_Northing = mNorthing_ERA[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
#                 m3x3_Distances = np.sqrt( (m3x3_ERA_Easting - mEasting_REMA[m,n])**2 + (m3x3_ERA_Northing - mNorthing_REMA[m,n])**2 )
#                 m3x3_Weights = m3x3_Distances / sum(m3x3_Distances[:])
#                 # Calculate weighted mean (based on distance) of surrounding ERA5 values
#                 temp = m3x3_Weights * wind_u_10m[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
#                 mWindU_hd[m,n] = sum(temp[:])
#                 temp = m3x3_Weights * wind_v_10m[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
#                 mWindV_hd[m,n] = sum(temp[:])
#                 temp = m3x3_Weights * precip[iERA_e_idx-1:iERA_e_idx+1,iERA_n_idx-1:iERA_n_idx+1]
#                 mPrecip_hd[m,n] = sum(temp[:])
                
#             except:
#                 pass
            
#     return

# %%
