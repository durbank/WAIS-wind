# A module script to define custom functions, variables, etc.
# This should basically replace the `REMA_import` and 
# `PAIPR_import` scripts (anything else will likely be 
# recycled into the various notebooks)


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